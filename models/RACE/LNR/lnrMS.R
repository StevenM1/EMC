# LNR mu and sigma parameterization

source("models/RACE/LNR/lnr.R")

lnrMS <- list(
  type="RACE",
  p_types=c("m","s","t0"),
  Ntransform=function(x) {
    # Transform to natural scale
   x[,dimnames(x)[[2]] != "m"] <- exp(x[,dimnames(x)[[2]] != "m"])
    x
  },
  # p_vector transform scaling parameter by s=1 assumed in lnr.R
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rLNR(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dLNR(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pLNR(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

