# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)

source("models/RACE/RDM/rdm.R")

rdmB <- list(
  type="RACE",
  p_types=c("v","B","A","t0","s"),
  # Transform to natural scale
  Ntransform=function(x) {
    exp(x)
  },
  # mapped parameter transform
  Mtransform = function(pars,dadm=NULL) 
    # transform parameters back to real line 
    # pars is a matrix output by map_p_vector  
  {
    pars
  },
  # p_vector transform 
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rRDM(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dRDM(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pRDM(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

