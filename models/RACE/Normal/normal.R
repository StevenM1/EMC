# Normal, mu and log(sigma) parameterization

source("models/RACE/Normal/norm.R")

normal <- list(
  type="RACE",
  p_types=c("mean","sd"),
  # Transform to natural scale
  Ntransform=function(x) {
    x[,dimnames(x)[[2]] == "sd"] <- exp(x[,dimnames(x)[[2]] == "sd"])
    x
  },
  # p_vector transform
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rNORMAL(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dNORMAL(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pNORMAL(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)
