# Normal, mu and log(sigma) parameterization

source("models/RACE/exGaussian/exG.R")

exGaussian <- list(
  type="RACE",
  p_types=c("mu","sigma","tau"),
  Ntransform=function(x) {
  # Transform to natural scale
    exp(x)
  },
  # p_vector transform
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rexGaussian(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dexGaussian(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pexGaussian(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

