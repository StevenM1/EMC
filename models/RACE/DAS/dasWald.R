# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)

source("models/RACE/DAS/DAS.R")

dasWald <- list(
  type="RACE",
  p_types=c("v","B","tau","l","t0","s"),
  # Transform to natural scale
  Ntransform=function(x) {
    pars <- exp(x)
    pars$vp <- pars$v + pars$l
    pars$vp[]
  },
  # p_vector transform 
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rDAS(lR,pars,rtype="wald"),
  # Density function (PDF) of choice race runner
  dfun=f_wald,
  # Density function (PDF) of choice race runner
  pfun=F_wald,
  # Density function (PDF) winner convolved with exponential
  dDAS=dDAS,
  # Race likelihood looping dDAS over trials
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_DAS(p_vector=p_vector, dadm = dadm, min_ll = min_ll,
                       pnames=c("B","v")) # must be in B, v order
)

