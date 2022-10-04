dexGaussianS <- function(rt,pars)
{
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

pexGaussianS <- function(rt,pars)
{
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

rSSexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau","muS","sigmaS","tauS"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
  #
  # test
  # pars=cbind(mu=c(.5,.6),sigma=c(.1,.1),tau=c(.2,.2)); lR=factor(c(1))
{
  if (!any(dimnames(pars)[[2]]=="SSD")) pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (!any(levels(lR)=="stop")) stop("lR must have a \"stop\" level")
  isstop <- lR=="stop"
  dt <- matrix(nrow=length(levels(lR)),ncol=dim(pars)[1]/length(levels(lR)))
  nstop <- sum(isstop)
  dt[isstop] <- rnorm(nstop,mean=pars[isstop,"muS"],sd=pars[isstop,"sigmaS"]) +
                rexp(nstop,rate=1/pars[isstop,"tauS"])
  ngo <- sum(!isstop)
  dt[!isstop] <- rnorm(ngo,mean=pars[!isstop,"mu"],sd=pars[!isstop,"sigma"]) +
                 rexp(ngo,rate=1/pars[!isstop,"tau"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}
