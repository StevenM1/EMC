dexGaussianG <- function(rt,pars) 
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

pexGaussianG <- function(rt,pars) 
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}
 
 
dRDEXG <- function(rt,pars) 
{
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dRDM(rt[ok],pars[ok,,drop=FALSE])
  out
}


pRDEXG <- function(rt,pars) 
{
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pRDM(rt[ok],pars[ok,,drop=FALSE])
  out
}



dexGaussianS <- function(rt,pars)
{

  dexGaussS <- function(rt,pars)
  {
    rt <- rt - pars[,"SSD"]
    isexp <- pars[,"sigmaS"] < 1e-4 # shifted exponential
    rt[isexp] <- dexp(rt[isexp]-pars[isexp,"muS"],1/pars[isexp,"tauS"]) 
    isnorm <- !isexp & pars[,"tauS"] < 0.05 * pars[,"sigmaS"] # normal
    rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"muS"], sd = pars[isnorm,"sigmaS"])
    isexg <- !(isexp | isnorm) 
    if (any(isexg)) {
      s2 <- pars[isexg,"sigmaS"]^2
      z <- rt[isexg] - pars[isexg,"muS"] - (s2/pars[isexg,"tauS"])
      rt[isexg] <- exp(
        log(pnorm(z/pars[isexg,"sigmaS"])) -
        log(pars[isexg,"tauS"]) - 
        (z + (s2/(2 *  pars[isexg,"tauS"])))/pars[isexg,"tauS"]
      )
    }
    rt
  }
  
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussS(rt[ok],pars[ok,,drop=FALSE])
  out
}


pexGaussianS <- function(rt,pars)
{
  
  pexGaussS <- function(rt,pars)
  {
    rt <- rt - pars[,"SSD"]
    isexp <- pars[,"sigmaS"] < 1e-4 # shifted exponential
    rt[isexp] <- pexp(rt[isexp]-pars[isexp,"muS"],1/pars[isexp,"tauS"]) 
    isnorm <- !isexp & pars[,"tauS"] < 0.05 * pars[,"sigmaS"] # normal
    rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"muS"], sd = pars[isnorm,"sigmaS"])
    isexg <- !(isexp | isnorm) 
    if (any(isexg)) {
      s2 <- pars[isexg,"sigmaS"]^2
      z <- rt[isexg] - pars[isexg,"muS"] - (s2/pars[isexg,"tauS"])
      rt[isexg] <- 
        pnorm((rt[isexg] - pars[isexg,"muS"])/pars[isexg,"sigmaS"]) - 
        exp(log(pnorm(z/pars[isexg,"sigmaS"])) + 
          ((pars[isexg,"muS"] + (s2/pars[isexg,"tauS"]))^2 - (pars[isexg,"muS"]^2) - 
             2 * rt[isexg] * (s2/pars[isexg,"tauS"]))/(2 * s2))
    }
    rt
  }

  
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussS(rt[ok],pars[ok,])
  out
}


update_ssd <- function(isstop,idx,idx1,ssd,stairstep,stairmin,stairmax) 
{
  if (isstop) {
    if (ssd[idx]+ stairstep < stairmax) 
      ssd[idx1] <- ssd[idx] + stairstep else ssd[idx1] <- ssd[idx]
    } else {
      if (ssd[idx] - stairstep > stairmin)
          ssd[idx1] <- ssd[idx] - stairstep else ssd[idx1] <- ssd[idx]
    }
  ssd
}

# lR=data$lR; staircase0=.2; stairstep=.05; stairmin=0; stairmax=Inf
# p_types=c("mu","sigma","tau","muS","sigmaS","tauS","tf","gf")

rSSexGaussian <- function(lR,pars,staircase0,stairstep,stairmin,stairmax,
                          p_types=c("mu","sigma","tau","muS","sigmaS","tauS","tf","gf"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. staircase0 is first SSD in staircase,stairstep is up or
  # down step, stairmin and stairmax are bounds on SSD
{
  
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (!any(levels(lR)=="stop")) stop("lR must have a \"stop\" level")
  # Go failures
  nacc <- length(levels(lR))
  ntrials <- dim(pars)[1]/nacc
  isstop <- lR=="stop"
  isgf <- pars[isstop,"gf"] > runif(ntrials)
  R <- factor(rep("stop",ntrials),levels=levels(lR)) 
  rt <- rep(NA,ntrials)
  # Not go failures
  nsample <- sum(!isgf)
  dt <- matrix(Inf,nrow=nacc,ncol=nsample)
  isgo <- rep(!isgf,each=nacc) & !isstop
  ngo <- sum(isgo)
  dt[-1,] <- rnorm(ngo,mean=pars[isgo,"mu"],sd=pars[isgo,"sigma"]) +
             rexp(ngo,rate=1/pars[isgo,"tau"]) 
  SSD <- pars[isstop,"SSD"]
  finiteSSD <- !is.infinite(SSD)
  gostopruns <-  c(finiteSSD & (pars[isstop,"tf"] < runif(ntrials)))[!isgf] # tfs
  if (any(!gostopruns)) { # Run go and trigger failure trials
    r <- apply(dt[,!gostopruns,drop=FALSE],2,which.min)
    R[!isgf][!gostopruns] <- levels(lR)[r]
    pick <- cbind(r,c(1:dim(dt)[2])[!gostopruns]) # Matrix to pick winner
    rt[!isgf][!gostopruns] <- dt[pick]
  }
  if (any(finiteSSD)) { # Run stop trials (and update staircase for gf and tf)
    nstop <- sum(gostopruns)
    dt[1,gostopruns] <- rexp(nstop,rate=1/pars[isstop,"tauS"][!isgf][gostopruns]) +
        rnorm(nstop,mean=pars[isstop,"muS"][!isgf][gostopruns],sd=pars[isstop,"sigmaS"][!isgf][gostopruns])
    isfixed <- gostopruns & !is.na(SSD[!isgf])
    if (any(isfixed)) { # Run fixed trials
      dt[1,isfixed] <- dt[1,isfixed] + SSD[!isgf][isfixed]
      r <- apply(dt[,isfixed],2,which.min)
      R[!isgf][isfixed] <- levels(lR)[r]
      pick <- cbind(r,c(1:dim(dt)[2])[isfixed]) # Matrix to pick winner
      rt[!isgf][isfixed] <- dt[pick]
    }
    isstaircase <- is.na(SSD)
    if (any(isstaircase)) {
      sindex <- c(1:ntrials)[isstaircase]
      stairgf <- sindex %in% c(1:ntrials)[isstaircase & isgf]
      dtindex <- c(1:sum(!isgf))[isstaircase[!isgf]]
      SSD[sindex[1]] <- staircase0
      dti <- 1
      for (i in 1:length(sindex)) {
        if (stairgf[i]) SSD <- update_ssd(TRUE,
          sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax) else 
        {
          dt[1,dtindex[dti]] <- dt[1,dtindex[dti]] + SSD[sindex[i]]
          pick <- which.min(dt[,dtindex[dti]])
          R[sindex[i]] <- levels(lR)[pick]
          rt[sindex[i]] <- dt[pick,dtindex[dti]]
          if ( i != length(sindex) ) {
            SSD <- update_ssd(R[sindex[i]]=="stop",
              sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax)
            dti <- dti+1
          }
        }
      }
    }
  }
  rt[R=="stop"] <- NA
  cbind.data.frame(R=R,rt=rt,SSD=SSD)
}


# lR=data$lR;staircase0=.2;stairstep=.05;stairmin=0;stairmax=Inf
# p_types=c("v","B","A","t0","muS","sigmaS","tauS","tf","gf")
rRDEX <- function(lR,pars,staircase0,stairstep,stairmin,stairmax,
                  p_types=c("v","B","A","t0","muS","sigmaS","tauS","tf","gf")) 
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. staircase0 is first SSD in staircase,stairstep is up or
  # down step, stairmin and stairmax are bounds on SSD
{
  
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (!any(levels(lR)=="stop")) stop("lR must have a \"stop\" level")
  # Go failures
  nacc <- length(levels(lR))
  ntrials <- dim(pars)[1]/nacc
  isstop <- lR=="stop"
  isgf <- pars[isstop,"gf"] > runif(ntrials)
  R <- factor(rep("stop",ntrials),levels=levels(lR)) 
  rt <- rep(NA,ntrials)
  # Not go failures
  nsample <- sum(!isgf)
  dt <- matrix(Inf,nrow=nacc,ncol=nsample)
  isgo <- rep(!isgf,each=nacc) & !isstop
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[isgo,c("A","B","v")] <- pars[isgo,c("A","B","v")]/pars[isgo,"s"]
  ngo <- sum(isgo)
  dt[-1,] <- pars[isgo,"t0"] + rWald(ngo,
      B=pars[isgo,"B"],v=pars[isgo,"v"],A=pars[isgo,"A"])
  SSD <- pars[isstop,"SSD"]
  finiteSSD <- !is.infinite(SSD)
  gostopruns <-  c(finiteSSD & (pars[isstop,"tf"] < runif(ntrials)))[!isgf] # tfs
  if (any(!gostopruns)) { # Run go and trigger failure trials
    r <- apply(dt[,!gostopruns,drop=FALSE],2,which.min)
    R[!isgf][!gostopruns] <- levels(lR)[r]
    pick <- cbind(r,c(1:dim(dt)[2])[!gostopruns]) # Matrix to pick winner
    rt[!isgf][!gostopruns] <- dt[pick]
  }
  if (any(finiteSSD)) { # Run stop trials (and update staircase for gf and tf)
    nstop <- sum(gostopruns)
    dt[1,gostopruns] <- rexp(nstop,rate=1/pars[isstop,"tauS"][!isgf][gostopruns]) +
        rnorm(nstop,mean=pars[isstop,"muS"][!isgf][gostopruns],sd=pars[isstop,"sigmaS"][!isgf][gostopruns])
    isfixed <- gostopruns & !is.na(SSD[!isgf])
    if (any(isfixed)) { # Run fixed trials
      dt[1,isfixed] <- dt[1,isfixed] + SSD[!isgf][isfixed]
      r <- apply(dt[,isfixed],2,which.min)
      R[!isgf][isfixed] <- levels(lR)[r]
      pick <- cbind(r,c(1:dim(dt)[2])[isfixed]) # Matrix to pick winner
      rt[!isgf][isfixed] <- dt[pick]
    }
    isstaircase <- is.na(SSD)
    if (any(isstaircase)) {
      sindex <- c(1:ntrials)[isstaircase]
      stairgf <- sindex %in% c(1:ntrials)[isstaircase & isgf]
      dtindex <- c(1:sum(!isgf))[isstaircase[!isgf]]
      SSD[sindex[1]] <- staircase0
      dti <- 1
      for (i in 1:length(sindex)) {
        if (stairgf[i]) SSD <- update_ssd(TRUE,
          sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax) else 
        {
          dt[1,dtindex[dti]] <- dt[1,dtindex[dti]] + SSD[sindex[i]]
          pick <- which.min(dt[,dtindex[dti]])
          R[sindex[i]] <- levels(lR)[pick]
          rt[sindex[i]] <- dt[pick,dtindex[dti]]
          if ( i != length(sindex) ) {
            SSD <- update_ssd(R[sindex[i]]=="stop",
              sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax)
            dti <- dti+1
          }
        }
      }
    }
  }
  rt[R=="stop"] <- NA
  cbind.data.frame(R=R,rt=rt,SSD=SSD)
}


my.integrate <- function(...,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (class(out)=="try-error") 0 else 
  {
    if (out$subdivisions==1) 
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (class(out)=="try-error") 0 else
      {
        if (out$subdivisions==1) 0 else out$value   
      }
    } else out$value
  }
}

  
pstopEXG <- function(parstop,n_acc,
  gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS")) 
{

  pEXG <- function (q, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) 
    # ex-Gaussian cumulative density
    # Modified from gamlss.dist to make cdf in tau > 0.05 * sigma case robust,
    # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
    # robust against small sigma cases.
  {
      if (sigma <= 0) return(rep(NA,length(q)))
      if (tau <= 0) return(rep(NA,length(q)))
    
      # if (sigma < 0.05*tau) 
      if (sigma < 1e-4) 
        return(pexp(q-mu,1/tau,log=log.p,lower.tail=lower.tail)) # shfited exponential
    
      ly <- length(q)
      sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)
      tau <- rep(tau, length = ly)
      index <- seq(along = q)
      z <- q - mu - ((sigma^2)/tau)
      cdf <- ifelse(is.finite(q), 
        ifelse(tau > 0.05 * sigma, 
           pnorm((q - mu)/sigma) - exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/tau))^2 - 
             (mu^2) -  2 * q * ((sigma^2)/tau))/(2 * sigma^2)), 
           pnorm(q, mean = mu, sd = sigma)),
          ifelse(q<0,0,1)   
      )
      if (lower.tail == TRUE) 
        cdf <- cdf
      else cdf <- 1 - cdf
      if (log.p == FALSE) 
        cdf <- cdf
      else cdf <- log(cdf)
      cdf
  }
  
  dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE) 
    # ex-Gaussian density
    # gamlss.dist function, but returns NA for bad sigma or tau, and
    # robust against small sigma cases.
  {
      if (sigma <= 0) return(rep(NA,length(x)))
      if (tau <= 0) return(rep(NA,length(x)))
    
      # if (sigma < 0.05*tau) 
      if (sigma < 1e-4) 
        return(dexp(x-mu,1/tau,log=log)) # shfited exponential
    
      ly <- length(x)
      sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)
      tau <- rep(tau, length = ly)
      z <- x - mu - ((sigma^2)/tau)
      logfy <- ifelse(tau > 0.05 * sigma, 
        -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)), 
        dnorm(x, mean = mu, sd = sigma, log = TRUE))
      if (log == FALSE) 
        fy <- exp(logfy)
      else fy <- logfy
      fy
  }

  dEXGrace <- function(dt,mu,sigma,tau)
    # Generates defective PDF for win by first runner, dt (decison time) is 
    # a matrix with length(mu) rows, one row for each runner, and one column
    # for each decision time for which a defective density value will be 
    # returned.
  {
    dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
    if (length(mu)>1) for (i in 2:length(mu))
      dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower.tail=FALSE)
    dt[1,]
  }
  
  stopfn <- function(t,mu,sigma,tau,SSD) 
    # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
  {
    dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
    dt[1,] <- dt[1,]-SSD
    dEXGrace(dt,mu,sigma,tau)
  }

  
  sindex <- seq(1,dim(parstop)[1],by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
  ntrials <- length(SSDs)
  pgo <- array(parstop[-sindex,gpars],dim=c(n_acc-1,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- character(ntrials)
  for (i in 1:ntrials) # choice parameter order doesn't matter
    cells[i] <- paste(SSDs[i],ps[i,],sort(apply(pgo[,i,],1,paste,collapse="")),collapse="")
  uniq <- !duplicated(cells)
  ups <- numeric(sum(uniq))
  for (i in 1:length(ups)) {
    ups[i] <- my.integrate(f=stopfn,lower=-Inf,SSD=SSDs[uniq][i],
      mu=c(ps[uniq,"muS"][i],pgo[,uniq,"mu",drop=FALSE][,i,1]),
      sigma=c(ps[uniq,"sigmaS"][i],pgo[,uniq,"sigma",drop=FALSE][,i,1]),
      tau=c(ps[uniq,"tauS"][i],pgo[,uniq,"tau",drop=FALSE][,i,1]))    
  }
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


pstopEXG <- function(parstop,n_acc,
  gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS")) 
{

  pEXG <- function (q, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) 
    # ex-Gaussian cumulative density
    # Modified from gamlss.dist to make cdf in tau > 0.05 * sigma case robust,
    # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
    # robust against small sigma cases.
  {
      if (sigma <= 0) return(rep(NA,length(q)))
      if (tau <= 0) return(rep(NA,length(q)))
    
      # if (sigma < 0.05*tau) 
      if (sigma < 1e-4) 
        return(pexp(q-mu,1/tau,log=log.p,lower.tail=lower.tail)) # shfited exponential
    
      ly <- length(q)
      sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)
      tau <- rep(tau, length = ly)
      index <- seq(along = q)
      z <- q - mu - ((sigma^2)/tau)
      cdf <- ifelse(is.finite(q), 
        ifelse(tau > 0.05 * sigma, 
           pnorm((q - mu)/sigma) - exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/tau))^2 - 
             (mu^2) -  2 * q * ((sigma^2)/tau))/(2 * sigma^2)), 
           pnorm(q, mean = mu, sd = sigma)),
          ifelse(q<0,0,1)   
      )
      if (lower.tail == TRUE) 
        cdf <- cdf
      else cdf <- 1 - cdf
      if (log.p == FALSE) 
        cdf <- cdf
      else cdf <- log(cdf)
      cdf
  }
  
  dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE) 
    # ex-Gaussian density
    # gamlss.dist function, but returns NA for bad sigma or tau, and
    # robust against small sigma cases.
  {
      if (sigma <= 0) return(rep(NA,length(x)))
      if (tau <= 0) return(rep(NA,length(x)))
    
      # if (sigma < 0.05*tau) 
      if (sigma < 1e-4) 
        return(dexp(x-mu,1/tau,log=log)) # shfited exponential
    
      ly <- length(x)
      sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)
      tau <- rep(tau, length = ly)
      z <- x - mu - ((sigma^2)/tau)
      logfy <- ifelse(tau > 0.05 * sigma, 
        -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)), 
        dnorm(x, mean = mu, sd = sigma, log = TRUE))
      if (log == FALSE) 
        fy <- exp(logfy)
      else fy <- logfy
      fy
  }

  dEXGrace <- function(dt,mu,sigma,tau)
    # Generates defective PDF for win by first runner, dt (decison time) is 
    # a matrix with length(mu) rows, one row for each runner, and one column
    # for each decision time for which a defective density value will be 
    # returned.
  {
    dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
    if (length(mu)>1) for (i in 2:length(mu))
      dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower.tail=FALSE)
    dt[1,]
  }
  
  stopfn <- function(t,mu,sigma,tau,SSD) 
    # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
  {
    dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
    dt[1,] <- dt[1,]-SSD
    dEXGrace(dt,mu,sigma,tau)
  }

  
  sindex <- seq(1,dim(parstop)[1],by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
  ntrials <- length(SSDs)
  pgo <- array(parstop[-sindex,gpars],dim=c(n_acc-1,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- character(ntrials)
  for (i in 1:ntrials)
    cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],collapse="")
  uniq <- !duplicated(cells)
  ups <- numeric(sum(uniq))
  for (i in 1:length(ups)) {
    ups[i] <- my.integrate(f=stopfn,lower=-Inf,SSD=SSDs[uniq][i],
      mu=c(ps[uniq,"muS"][i],pgo[,uniq,"mu",drop=FALSE][,i,1]),
      sigma=c(ps[uniq,"sigmaS"][i],pgo[,uniq,"sigma",drop=FALSE][,i,1]),
      tau=c(ps[uniq,"tauS"][i],pgo[,uniq,"tau",drop=FALSE][,i,1]))    
  }
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


pstopRDEX <- function(parstop,n_acc,
  gpars=c("v","B","A","t0"),spars=c("muS","sigmaS","tauS")) 
{

  pWald <- function(t,v,B,A,t0)
    # cumulative density for single accumulator
  {
    pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
      # cdf of inverse gaussian at t with k +/- a/2 uniform variability
      # returns pigt.0 if a<=0
    
      pigt.0 <- function(t,k=1,l=1) {
        # cdf of inverse gaussian at t with no k variability
        # much faster than statmod's pinvgauss funciton
      
        mu <- k/l
        lambda <- k^2
      
        e <- exp(log(2*lambda) - log(mu))
        add <- sqrt(lambda/t) * (1 + t/mu)
        sub <- sqrt(lambda/t) * (1 - t/mu)
      
        p.1 <- 1 - pnorm(add)
        p.2 <- 1 - pnorm(sub)
        x <- exp(e + log(p.1)) + p.2
      
        x[t<0] <- 0
        x
      }
    
      options(warn=-1)
      if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
      if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
      if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    
      tpos <- t<=0
    
      atiny <- a<=tiny & !tpos
      a[atiny] <- 0
    
      ltiny <- (l<=tiny) & !atiny & !tpos
      notltiny <- (l>tiny) & !atiny & !tpos
      l[l<=tiny] <- 0
    
      x <- numeric(length(t))
    
      # No threshold variability
      if ( any(atiny) )
        x[atiny] <- pigt.0(t[atiny],k[atiny],l[atiny])
    
      # Threshold variability
      if ( any(!atiny) ) {
      
        if ( any(notltiny) ) { # rate non-zero
        
          log.t <- log(t[notltiny])
          sqr.t <- sqrt(t[notltiny])
        
          term.1a <- .5*log.t-.5*log(2*pi)
          term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
          term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
          term.1 <- exp(term.1a)*(term.1b-term.1c)
        
          term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) + 
            log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
          term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) + 
            log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
          term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
          
          term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
          term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
          term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
          term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
          term.4 <- term.4c*term.4a + term.4d*term.4b
        
          x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
        }
      
        if ( any(ltiny) ) {  # rate zero
          sqr.t <- sqrt(t[ltiny])
          log.t <- log(t[ltiny])
          term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
          term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
          term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
          
          term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
          term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
          term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
          x[ltiny] <- term.5 + term.6
        }
      
      }
    
      x[x<0 | is.nan(x) ] <- 0
      x
    }
  
    out <- numeric(length(t))
    t <- t-t0[1]
    ok <- v>0
    out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
    out[!ok] <- 0
    out
  
  }

  dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE) 
    # ex-Gaussian density
    # gamlss.dist function, but returns NA for bad sigma or tau, and
    # robust against small sigma cases.
  {
      if (sigma <= 0) return(rep(NA,length(x)))
      if (tau <= 0) return(rep(NA,length(x)))
    
      # if (sigma < 0.05*tau) 
      if (sigma < 1e-4) 
        return(dexp(x-mu,1/tau,log=log)) # shfited exponential
    
      ly <- length(x)
      sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)
      tau <- rep(tau, length = ly)
      z <- x - mu - ((sigma^2)/tau)
      logfy <- ifelse(tau > 0.05 * sigma, 
        -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)), 
        dnorm(x, mean = mu, sd = sigma, log = TRUE))
      if (log == FALSE) 
        fy <- exp(logfy)
      else fy <- logfy
      fy
  }

  dRDEXrace <- function(dt,mu,sigma,tau,v,B,A,t0)
    # Generates defective PDF for win by first runner, dt (decison time) is 
    # a matrix with length(mu) rows, one row for each runner, and one column
    # for each decision time for which a defective density value will be 
    # returned.
  {
    dt[1,] <- dEXG(dt[1,],mu,sigma,tau)
    for (i in 1:length(v))
      dt[1,] <- dt[1,]*(1-pWald(dt[i+1,],v[i],B[i],A[i],t0[i]))
    dt[1,]
  }
  
  stopfn <- function(t,n_acc,mu,sigma,tau,v,B,A,t0,SSD) 
    # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
  {
    dt <- matrix(rep(t+SSD,each=n_acc),nrow=n_acc)
    dt[1,] <- dt[1,]-SSD
    dRDEXrace(dt,mu,sigma,tau,v,B,A,t0)
  }

  
  if (any(dimnames(parstop)[[2]]=="s")) # rescale
     parstop[,c("A","B","v")] <- parstop[,c("A","B","v")]/parstop[,"s"]
  sindex <- seq(1,dim(parstop)[1],by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
  ntrials <- length(SSDs)
  pgo <- array(parstop[-sindex,gpars],dim=c(n_acc-1,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- character(ntrials)
  for (i in 1:ntrials)
    cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],collapse="")
  uniq <- !duplicated(cells)
  ups <- numeric(sum(uniq))
  for (i in 1:length(ups)) {
    ups[i] <- my.integrate(f=stopfn,lower=0,
      n_acc=n_acc,SSD=SSDs[uniq][i],
      mu=ps[uniq,"muS"][i],
      sigma=ps[uniq,"sigmaS"][i],
      tau=ps[uniq,"tauS"][i],
      v=pgo[,uniq,"v",drop=FALSE][,i,1],
      B=pgo[,uniq,"B",drop=FALSE][,i,1],
      A=pgo[,uniq,"A",drop=FALSE][,i,1],
      t0=pgo[,uniq,"t0",drop=FALSE][,i,1])    
  }
  
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


log_likelihood_race_ss <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood, one accumulator (lR=="stop") is ExGaussian
  # others are Wald.
{
  
  pars <- get_pars(p_vector,dadm)
  
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  
  # standard race, ignores R=NA
  isstop <- dadm[,"lR"]=="stop"
  stopwinner <- dadm$winner & isstop
  gowinner <- dadm$winner & !isstop
  stoplooser <- !dadm$winner & isstop
  golooser <- !dadm$winner & !isstop
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[gowinner] <- log(attr(dadm,"model")$dfunG(rt=dadm$rt[gowinner],pars=pars[gowinner,,drop=FALSE]))
  lds[stopwinner] <- log(attr(dadm,"model")$dfunS(rt=dadm$rt[stopwinner],pars=pars[stopwinner,,drop=FALSE]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) {
    lds[golooser] <- log(1-attr(dadm,"model")$pfunG(rt=dadm$rt[golooser],pars=pars[golooser,,drop=FALSE]))
    lds[stoplooser] <- log(1-attr(dadm,"model")$pfunS(rt=dadm$rt[stoplooser],pars=pars[stoplooser,,drop=FALSE]))
  }
  lds[is.na(lds) | !ok] <- 0
  # lds[is.na(lds)] <- 0
  lds <- lds[attr(dadm,"expand")] # decompress
  winner <- dadm$winner[attr(dadm,"expand")]
  if (n_acc>1) {
    ll <- lds[winner]
    if (n_acc==2)
      ll <- ll + lds[!winner] else
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    ll[is.na(ll)] <- 0
  } else ll <- lds
  like <- exp(ll)
  
  # stop success
  stopsucess_accs <- is.finite(dadm$SSD[attr(dadm,"expand")]) & (dadm$R[attr(dadm,"expand")] == "stop")
  stoptrial <- is.finite(dadm$SSD[attr(dadm,"expand")][winner])
  goresp <- dadm$R[attr(dadm,"expand")][winner] != "stop"
  if (any(stopsucess_accs)) {
    parstop <- pars[attr(dadm,"expand"),][stopsucess_accs,]
    like[stoptrial & !goresp] <- attr(dadm,"model")$sfun(parstop,n_acc)
  }
  
  # trigger failures  
  tf <- pars[attr(dadm,"expand"),"tf"][winner]
  dotf <- goresp & stoptrial & (tf > 0)
  if (any(dotf)) {
    tflike <- exp(apply(matrix(lds,nrow=n_acc)[-1,dotf,drop=FALSE],2,sum))
    like[dotf] <- tf[dotf]*tflike + (1-tf[dotf])*like[dotf]
    stopsucess <- !dotf & stoptrial
    like[stopsucess] <- like[stopsucess]*(1-tf[stopsucess])
  }

  # go failures
  gf <- pars[,"gf"][attr(dadm,"expand")][winner]
  if (any(gf>0)) {
    like[!goresp] <- gf[!goresp] + (1-gf[!goresp])*like[!goresp]
    like[goresp] <- like[goresp]*(1-gf[goresp]) 
  }
  
  like[!ok[winner]] <- 0
  
  sum(pmax(min_ll,log(like)))

  # # trigger failures
  # tf <- pars[attr(dadm,"expand"),"tf"][winner]
  # stoptrial <- is.finite(dadm$SSD[attr(dadm,"expand")][winner])
  # goresp <- dadm$R[attr(dadm,"expand")][winner] != "stop"
  # dotf <- goresp & stoptrial & (tf > 0)
  # if (any(dotf)) {
  #   tflike <- exp(apply(matrix(lds,nrow=n_acc)[-1,dotf,drop=FALSE],2,sum))
  #   like[dotf] <- tf[dotf]*tflike + (1-tf[dotf])*like[dotf]
  #   stopsucess <- !dotf & stoptrial
  #   like[stopsucess] <- like[stopsucess]*(1-tf[stopsucess]) 
  # }
  # 
  # 
  # # go failures
  # gf <- pars[,"gf"][attr(dadm,"expand")][winner]
  # isgf <- is.na(dadm$rt[attr(dadm,"expand")][winner])
  # like[isgf] <- gf[isgf]  
  # like[!isgf] <- like[!isgf]*(1-gf[!isgf])
  # sum(pmax(min_ll,log(like)))
}

