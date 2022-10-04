rRDEX <- function(lR,pars,p_types=c("v","B","A","t0","muS","sigmaS","tauS")) 
  # lR is an latent response factor with one level for each accumulator. 
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in 
  # contiguous rows. "s" parameter will be used but can be omitted
  #
  # test
  # pars=cbind(B=c(NA,1,1),v=c(NA,1,1),A=c(NA,0,0),t0=c(NA,.2,.2),
  #            mu=rep(.5,3),sigma=rep(.05,3),tau=rep(.1,3))
  # lR=factor(c("stop","left","right"))
  # pars <- rbind(pars,pars); lR <- c(lR,lR)
{
  if (!any(dimnames(pars)[[2]]=="SSD")) pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
  if (!all(p_types %in% dimnames(pars)[[2]])) 
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (!any(levels(lR)=="stop")) stop("lR must have a \"stop\" level")
  isstop <- lR=="stop"
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[!isstop,c("A","B","v")] <- pars[!isstop,c("A","B","v")]/pars[!isstop,"s"]
  dt <- matrix(nrow=length(levels(lR)),ncol=dim(pars)[1]/length(levels(lR)))
  nstop <- sum(isstop)
  dt[isstop] <- rnorm(nstop,mean=pars[isstop,"muS"],sd=pars[isstop,"sigmaS"]) +
                 rexp(nstop,rate=1/pars[isstop,"tauS"])
  dt[!isstop] <- pars[!isstop,"t0"] + rWald(sum(!isstop),
    B=pars[!isstop,"B"],v=pars[!isstop,"v"],A=pars[!isstop,"A"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}



