# Discrete Activation-suppression normal or wald race

#### Wald ----


f_wald <- function(t,a,v) 
  # a = threshold, v = rate
  (a/sqrt(2*pi*t^3))*exp(-((a-v*t)^2)/(2*t))


F_wald <- function(t,a,v) 
  # a = threshold, v = rate 
  # Uses Derenzo's approx for large negative arguments to the second phi
  # v=0;t=10;a=11;s=1;(v*t+a)/(s*sqrt(t))
{
  a <- rep(a,length.out=length(t))
  v <- rep(v,length.out=length(t))
  ssrt <- sqrt(t)
  vt <- v*t
  ok <- (vt+a)/ssrt < 5.5
  out <- numeric(length(t))
  if (any(ok)) out[ok] <- pnorm((vt[ok]-a[ok])/ssrt[ok]) + 
    exp(2*a[ok]*v[ok])*pnorm(-(vt[ok]+a[ok])/ssrt[ok]) 
  if (any(!ok)) {
    s2t <- t[!ok]
    out[!ok] <- pnorm((vt[!ok]-a[!ok])/ssrt[!ok]) + 
      (ssrt[!ok]/(sqrt(2*pi)*(vt[!ok]+a[!ok])))*exp(-((vt[!ok]-a[!ok])^2/(2*s2t))-0.94*s2t/(vt[!ok]+a[!ok])^2)
    
  }
  out
}  


#### DAS distribution functions ----

# rt=rt[i];tau=pwin[i,3:4];pwin=pwin[i,1:2];ploose=matrix(ploose[,i,],ncol=2)
# dfun=f_wald; pfun=F_wald
# dfun=attr(dadm,"model")$dfun;pfun=attr(dadm,"model")$pfun
dDAS <- function(rt,pwin,ploose,tau,dfun,pfun,soa=0,rel.tol = .Machine$double.eps^0.25) 
  # Winner of race convolved with exponential
  # pwin 3 vector of threshold, no interference (-) and interference (+) rates
  # ploose like pwin, but a matrix with as many rows as loosing accumulators
  # tau is exponential mean for suppression and encoding processes
{
  
  # Race equation
  race <- function(rt,pwin,ploose,dfun,pfun) {
    out <- dfun(rt,pwin[1],pwin[2])*(1-pfun(rt,ploose[1,1],ploose[1,2]))
    if (dim(ploose)[1]>1) for (i in 2:dim(ploose)[1]) 
      out <- out*(1-pfun(rt,ploose[i,1],ploose[i,2]))  
    out
  }
  # race(1:10,pwin,ploose,dfun,pfun)
  
  # Convolution 
  fun <- function(w,rt,pwin,ploose,tau,dfun,pfun)
    race(w,pwin,ploose,dfun,pfun)*dexp(rt-w,1/tau)

  # # Uncomment and comment out lines below to test speed relative to simple Wald race
  # race(rt,pwin,ploose,dfun,pfun) # normal ~5 times quicker, wald ~10 times faster
  
  # rates
  a <- 1/tau[1]
  b <- 1/tau[2]
  ab <- a+b
  tauAB <- 1/ab
  if (soa==0) eas <- 1 else eas <- exp(-a*soa)
  abeasc <- a+b*(1-eas)
  q <- (b/ab)*eas # probability of interference
 
  # No interference component
  if (q > .999) f_minus <- 0 else {
    c_minusAB <- try(
      integrate(fun,0,Inf,rt=rt,pwin=pwin[c(1,2)],ploose=ploose[,c(1,2),drop=FALSE],
                tau=tauAB,dfun=dfun,pfun=pfun,rel.tol=rel.tol),
    silent=TRUE)
    if (class(c_minusAB)=="try-error") return(0)
    c_minusB <- try(
      integrate(fun,0,Inf,rt=rt,pwin=pwin[c(1,2)],ploose=ploose[,c(1,2),drop=FALSE],
                tau=tau[2],dfun=dfun,pfun=pfun,rel.tol=rel.tol),
    silent=TRUE)
    if (class(c_minusB)=="try-error") return(0)
    f_minus <- (ab/abeasc)*c_minusB$value - (b*eas/abeasc)*c_minusAB$value  
  }
  # Interference component
  if (q < .001) f_plus <- 0 else {
    if (pwin[2]==pwin[3] & q < .999) f_plus <- c_minusAB else
      f_plus <- try( integrate(fun,0,Inf,rt=rt,pwin=pwin[c(1,3)],
        ploose=ploose[,c(1,3),drop=FALSE],tau=tauAB,dfun=dfun,pfun=pfun,rel.tol=rel.tol),
    silent=TRUE) 
  }
  if (class(f_plus)=="try-error") return(0)
  # combine plus and minus
  q*f_plus$value + (1-q)*f_minus
}


# a1 <- a2 <- 1; v1 <- 3; l1 <- 1; v2 <- 2; l2 <- -1; tauA <- .25; tauB <- .5
# l1=l2=0 # DASC0
# pwin <- c(a1,v1,v1+l1); ploose <- matrix(c(a2,v2,v2+l2),nrow=1) # 2 choice
# ploose <- matrix(c(a2,a2,v2,v2+l2,v2+l2),nrow=2)                # 3 choice
# tau <- c(tauA,tauB)
# dDAS(rt=1,pwin,ploose,tau,f_wald,F_wald)
# 
# # Check speed of big likelihood calculations
# # Guts of log_likelihood_DAS
# lguts <- function(dadm,pars,p_names,rel.tol = .Machine$double.eps^0.25) {
#   n_acc <- length(levels(dadm$R))
#   ntrials <- dim(dadm)[1]/n_acc
#   like <- numeric(ntrials)
#   pwin <- pars[dadm$winner,c(p_names,"tauA","tauB")]
#   ploose <- array(pars[!dadm$winner,p_names],dim=c(n_acc-1,ntrials,3))
#   if (!any(dimnames(pars)[[2]]=="t0")) rt <- dadm[dadm$winner,"rt"] else
#     rt <- dadm[dadm$winner,"rt"] - pars[dadm$winner,"t0"]
#   for (i in 1:ntrials)
#     like[i] <- dDAS(rt=rt[i],pwin=pwin[i,1:3],rel.tol=rel.tol,
#                     ploose=matrix(ploose[,i,],ncol=3),tau=pwin[i,4:5],
#                     dfun=attr(dadm,"model")$dfun,pfun=attr(dadm,"model")$pfun)
#   invisible(like)
# }
# n=1000; rt <- c(1:n)/n
# pars <- parsi <- cbind(rbind(pwin,ploose),tauA,tauB)        # 2 choice
# pars <- parsi <- cbind(rbind(pwin,ploose,ploose),tauA,tauB) # 3 choice
# dimnames(pars) <- list(NULL,c("a","vm","vp","tauA","tauB"))
# nacc <- dim(parsi)[1]
# for (i in 2:n) pars <- rbind(pars,parsi)
# winneri <- c(T,rep(F,nacc-1))
# dadm <- cbind.data.frame(rt=rep(rt,each=dim(parsi)[1]),winner=rep(winneri,times=n),
#                          R=factor(rep(1:nacc,times=n)))
# 
# # Wald race
# attr(dadm,"model") <- list(dfun=f_wald,pfun=F_wald)
# dimnames(pars)[[2]][1:3] <- c("a","vm","vp")
# system.time({res <- lguts(dadm,pars,p_names=c("a","vm","vp"))})
# system.time({res1 <- lguts(dadm,pars,p_names=c("a","vm","vp"),rel.tol=.001)})
# system.time({res2 <- lguts(dadm,pars,p_names=c("a","vm","vp"),rel.tol=.01)})
# system.time({res3 <- lguts(dadm,pars,p_names=c("a","vm","vp"),rel.tol=.1)})
# sum(log(res[res>0]))
# sum(log(res1[res1>0]))
# sum(log(res2[res2>0]))
# sum(log(res3[res3>0]))


#### DAS random function ----

rDAS <- function(lR,pars,rtype=c("wald")[1]) 
  # lR is an empty latent response factor lR with one level for each accumulator. 
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in 
  # contiguous rows. 
{

  rWald <- function(n,a,v,t0=0,s=NULL,dorep=FALSE)
  # random function for single accumulator
  {
  
    rwaldt <- function(n,a,v,tau,tiny=1e-6) {
      # random sample of n from a Wald (or Inverse Gaussian)
      # a = criterion, v = rate, assumes s=1 
    
      rlevy <- function(n=1, m=0, c=1) {
        if (any(c<0)) stop("c must be positive")
        c/qnorm(1-runif(n)/2)^2+m
      }
    
      flag <- v>tiny
      x <- rep(NA,times=n)
      x[!flag] <- rlevy(sum(!flag),0,v[!flag]^2)
      mu <- a/v
      lambda <- a^2
      y <- rnorm(sum(flag))^2
      mu.0 <- mu[flag]
      lambda.0 <- lambda[flag]
      x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
        sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)
      z <- runif(length(x.0))
      test <- mu.0/(mu.0+x.0)
      x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
      x[flag] <- x.0
      x[x<0] <- max(x)
      x
    }
  
    if (dorep) {
      a <- rep(a,length.out=n)
      v <- rep(v,length.out=n)
      tau <- rep(tau,length.out=n)
      t0 <- rep(t0,length.out=n)
      if (!is.null(s)) s <- rep(s,length.out=n)
    }
  
    if (!is.null(s)) { # rescale
      a <- a/s
      v <- v/s
    }

    # Kluge to return Inf for negative rates
    out <- numeric(n)
    ok <- !v<0  
    nok <- sum(ok)
    as <- a[ok]
    out[ok] <- rwaldt(nok,a=as,v=v[ok])
    out[!ok] <- Inf
    out
  }
  
  nacc <- length(levels(lR))
  n <- dim(pars)[1]/nacc
  eA <- rexp(n,1/pars[,"tauA"])
  eB <- rexp(n,1/pars[,"tauB"])
  is.intf <- eA > eB
  is.intf_nacc <- rep(is.intf,each=nacc)
  dt <- matrix(ncol=n,nrow=length(levels(lR)))
  if (rtype=="wald") {  
    pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives 
    dt[,is.intf] <- rWald(sum(is.intf),a=pars[is.intf_nacc,"B"],
                          v=pars[is.intf_nacc,"vp"],s=pars[is.intf_nacc,"s"])
    dt[,!is.intf] <- rWald(sum(!is.intf),a=pars[!is.intf_nacc,"B"],
                           v=pars[!is.intf_nacc,"v"],s=pars[!is.intf_nacc,"s"])
  } 
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt+eB)
}


# # Test Wald
# n=2e7
# p <- c(B=1,v=3,l=1,tauA=.5,tauB=.5,s=1,t0=.25)
# p <- c(B=1,v=4,l=0,tauA=.5,tauB=.5,s=1,t0=.25)
# p <- c(B=1,v=3,l=1,tauA=.01,tauB=.5,s=1,t0=.25)
# p <- c(B=1,v=3,l=1,tauA=10,tauB=.5,s=1,t0=.25)
# pars <- matrix(rep(p,each=n),nrow=n,dimnames=list(NULL,names(p)))
# pars[(1:(n/2))*2,"v"] <- 3
# pars[(1:(n/2))*2,"l"] <- -pars[(1:(n/2))*2,"l"]
# 
# # w_s > 0
# p <- c(B=1,v=3,l=2,tauA=.5,tauB=.5,s=1,t0=.25)
# pars <- matrix(rep(p,each=n),nrow=n,dimnames=list(NULL,names(p)))
# pars[(1:(n/2))*2,"v"] <- 3
# pars[(1:(n/2))*2,"l"] <- 1
# 
# pars <- cbind(pars,vp=pars[,"v"]+pars[,"l"])
# Rrt <- rDAS(factor(1:2),pars)
# pc <- mean(Rrt[,"R"]=="1")
# d1 <- density(Rrt[Rrt[,"R"]==1,"rt"]); d1$y <- d1$y*pc
# d2 <- density(Rrt[Rrt[,"R"]==2,"rt"]); d2$y <- d2$y*(1-pc)
# dadm <- cbind.data.frame(rt=rep(d1$x,each=2),R=factor(rep(1:2,times=length(d1$x))),
#                          winner=rep(c(TRUE,FALSE),times=length(d1$x)))
# attr(dadm,"expand") <- 1:length(d1$x)
# pmat <- pars[1:dim(dadm)[1],]
# attr(dadm,"model") <- list(dfun=f_wald,pfun=F_wald)
# p_names=c("B","v","vp")
# rel.tol=.Machine$double.eps^0.25
# D1 <- lguts(dadm,pmat,p_names,rel.tol=rel.tol)
# dadm$winner <- !dadm$winner
# D2 <- lguts(dadm,pmat,p_names,rel.tol=rel.tol)
# 
# 
# # Plot results
# plot(d1); lines(d2)
# lines(d1$x,D1,col="red")
# lines(d2$x,D2,col="red")

