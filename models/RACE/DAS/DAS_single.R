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

# rt=rt[i];tau=pwin[i,3];pwin=pwin[i,1:2];ploose=matrix(ploose[,i,],ncol=2)
# dfun=attr(dadm,"model")$dfun;pfun=attr(dadm,"model")$pfun
dDAS <- function(rt,pwin,ploose,tau,dfun,pfun) 
  # Winner of race convolved with exponential
  # pwin 2 vector of first and second parameters of dfun and pfun
  # ploose is a matrix with as many rows as loosing accumulators
  # tau is exponential mean
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
  
  # For full model normal ~2.5 times quicker than Wald
  out <- try(
    integrate(fun,0,Inf,rt=rt,pwin=pwin,ploose=ploose,tau=tau,dfun=dfun,pfun=pfun),
  silent=TRUE)
  if (class(out)=="try-error") 0 else out$value
}


  # # Race equation
  # race <- function(rt,pwin,ploose,dfun,pfun) {
  #   out <- dfun(rt,pwin[1],pwin[2])*(1-pfun(rt,ploose[1,1],ploose[1,2]))
  #   if (dim(ploose)[1]>1) for (i in 2:dim(ploose)[1]) 
  #     out <- out*(1-pfun(rt,ploose[i,1],ploose[i,2]))  
  #   out
  # }
# a1 <- a2 <- 1; v1 <- 4; v2 <- 2; tau <- .5
# pwin <- c(a1,v1); ploose <- matrix(c(a2,v2),nrow=1) # 2 choice
# ploose <- matrix(c(a2,a2,v2,v2),nrow=2)             # 3 choice
# # CHECK integrates to p_win
# dDASvec <- function(rt,pwin,ploose,tau,dfun,pfun) {
#   out <- numeric(length(rt))
#   for (i in 1:length(rt))
#     out[i] <- dDAS(rt[i],pwin,ploose,tau,dfun,pfun)
#   out
# }
# # Check Wald race
# dDAS(rt=1,pwin,ploose,tau,f_wald,F_wald)
# integrate(dDASvec,0,Inf,pwin=pwin,ploose=ploose,tau=tau,dfun=f_wald,pfun=F_wald)
# integrate(race,0,Inf,pwin=pwin,ploose=ploose,dfun=f_wald,pfun=F_wald)
# # Check normal race
# dDAS(rt=1,pwin,ploose,tau,dnorm,pnorm)
# integrate(dDASvec,0,Inf,pwin=pwin,ploose=ploose,tau=tau,dfun=dnorm,pfun=pnorm)
# integrate(race,0,Inf,pwin=pwin,ploose=ploose,dfun=dnorm,pfun=pnorm)


# # Check speed of big likelihood calculations
# n=1000; rt <- c(1:n)/n
# pars <- parsi <- cbind(rbind(pwin,ploose),c(tau,tau)) # 2 choice
# pars <- parsi <- cbind(rbind(pwin,ploose,ploose),c(tau,tau,tau)) # 3 choice
# dimnames(pars) <- list(NULL,c("a","v","tau"))
# nacc <- dim(parsi)[1]
# for (i in 2:n) pars <- rbind(pars,parsi)
# winneri <- c(T,rep(F,nacc-1))
# dadm <- cbind.data.frame(rt=rep(rt,each=dim(parsi)[1]),winner=rep(winneri,times=n),
#                          R=factor(rep(1:nacc,times=n)))
# 
# # Guts of log_likelihood_DAS
# lguts <- function(dadm,pars,p_names) {
#   n_acc <- length(levels(dadm$R))
#   ntrials <- dim(dadm)[1]/n_acc
#   like <- numeric(ntrials)
#   pwin <- pars[dadm$winner,c(p_names,"tau")]
#   ploose <- array(pars[!dadm$winner,p_names],dim=c(n_acc-1,ntrials,2))
#   if (!any(dimnames(pars)[[2]]=="t0")) rt <- dadm[dadm$winner,"rt"] else
#     rt <- dadm[dadm$winner,"rt"] - pars[dadm$winner,"t0"]
#   for (i in 1:ntrials)
#     like[i] <- dDAS(rt=rt[i],pwin=pwin[i,1:2],
#                     ploose=matrix(ploose[,i,],ncol=2),tau=pwin[i,3],
#                     dfun=attr(dadm,"model")$dfun,pfun=attr(dadm,"model")$pfun)
#   invisible(like)
# }
# 
# # Wald race
# attr(dadm,"model") <- list(dfun=f_wald,pfun=F_wald)
# dimnames(pars)[[2]][1:2] <- c("a","v")
# system.time({lguts(dadm,pars,p_names=c("a","v"))})
# # Normal race
# attr(dadm,"model") <- list(dfun=pnorm,pfun=dnorm)
# dimnames(pars)[[2]][1:2] <- c("mean","sd")
# system.time({lguts(dadm,pars,p_names=c("mean","sd"))})
# 


#### DAS random function ----

rDAS <- function(lR,pars,rtype=c("normal","wald")[1]) 
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
  
  if (rtype=="wald") 
    p_types=c("v","B","s","tau","t0") else
    p_types=c("mean","sd","tau","t0") # normal
  
  if (!all(p_types %in% dimnames(pars)[[2]])) 
    stop("pars must have columns ",paste(p_types,collapse = " "))
  
  if (rtype=="wald") {  
    pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives 
    dt <- matrix(rWald(dim(pars)[1],a=pars[,"B"],v=pars[,"v"],
                       s=pars[,"s"]),nrow=length(levels(lR)))
  } else  dt <- matrix(rnorm(dim(pars)[1],mean=pars[,"mean"],
                             sd=pars[,"sd"]),nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt+rexp(length(rt),1/pars[,"tau"]))
}


# # Test Wald
# n=1e7
# p <- c(B=1,v=5,tau=.5,s=1,t0=.25)
# pars <- matrix(rep(p,each=n),nrow=n,dimnames=list(NULL,names(p)))
# pars[(1:(n/2))*2,"v"] <- 3
# Rrt <- rDAS(factor(1:2),pars,rtype="wald")
# pc <- mean(Rrt[,"R"]=="1")
# d1 <- density(Rrt[Rrt[,"R"]==1,"rt"]); d1$y <- d1$y*pc
# d2 <- density(Rrt[Rrt[,"R"]==2,"rt"]); d2$y <- d2$y*(1-pc)
# dadm <- cbind.data.frame(rt=rep(d1$x,each=2),R=factor(rep(1:2,times=length(d1$x))),
#                          winner=rep(c(TRUE,FALSE),times=length(d1$x)))
# attr(dadm,"expand") <- 1:length(d1$x)
# pmat <- pars[1:dim(dadm)[1],]
# attr(dadm,"model") <- list(dfun=f_wald,pfun=F_wald)
# p_names=c("B","v")
# D1 <- lguts(dadm,pmat,p_names)
# dadm$winner <- !dadm$winner
# D2 <- lguts(dadm,pmat,p_names)
# 
# # Test Normal
# n=1e7
# p <- c(mean=1,sd=.2,tau=.5,s=1,t0=.25)
# pars <- matrix(rep(p,each=n),nrow=n,dimnames=list(NULL,names(p)))
# pars[(1:(n/2))*2,"mean"] <- 1.2
# Rrt <- rDAS(factor(1:2),pars,rtype="normal")
# pc <- mean(Rrt[,"R"]=="1")
# d1 <- density(Rrt[Rrt[,"R"]==1,"rt"]); d1$y <- d1$y*pc
# d2 <- density(Rrt[Rrt[,"R"]==2,"rt"]); d2$y <- d2$y*(1-pc)
# dadm <- cbind.data.frame(rt=rep(d1$x,each=2),R=factor(rep(1:2,times=length(d1$x))),
#                          winner=rep(c(TRUE,FALSE),times=length(d1$x)))
# attr(dadm,"expand") <- 1:length(d1$x)
# pmat <- pars[1:dim(dadm)[1],]
# attr(dadm,"model") <- list(dfun=dnorm,pfun=pnorm)
# p_names=c("mean","sd")
# D1 <- lguts(dadm,pmat,p_names)
# dadm$winner <- !dadm$winner
# D2 <- lguts(dadm,pmat,p_names)
# 
# 
# # Plot results
# plot(d1); lines(d2)
# lines(d1$x,D1,col="red")
# lines(d2$x,D2,col="red")

