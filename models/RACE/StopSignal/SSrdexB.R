# RDEX stop signal model
# GO = RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Stop = exGaussian

source("models/RACE/StopSignal/ss.R")
source("models/RACE/RDM/rdm.R")

SSrdexB <- list(
  type="RACE",
  p_types=c("v","B","A","t0","s","muS","sigmaS","tauS","tf","gf"),
  # Transform to natural scale
  Ntransform=function(x) {
    # transform parameters back to real line 
    isprobit <- dimnames(x)[[2]] %in% c("tf","gf")
    x[,!isprobit] <- exp(x[,!isprobit])
    x[,isprobit] <- pnorm(x[,isprobit])
    x
  },
  # p_vector transform 
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
                                 pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) & 
      (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) & 
      ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)
    pars
  },
  # Density function (PDF) for single go racer
  dfunG=function(rt,pars) dRDEXG(rt,pars),
  # Probability function (CDF) for single go racer
  pfunG=function(rt,pars) pRDEXG(rt,pars),
  # Density function (PDF) for single stop racer
  dfunS=function(rt,pars) dexGaussianS(rt,pars),
  # Probability function (CDF) for single stop racer
  pfunS=function(rt,pars) pexGaussianS(rt,pars),
  # Stop probability integral
  sfun=function(pars,n_acc) pstopRDEX(pars,n_acc),
  # Random function for SS race 
  rfun=function(lR,pars) rRDEX(lR,pars,
    staircase0=.2,stairstep=.05,stairmin=0,stairmax=Inf),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race_ss(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

