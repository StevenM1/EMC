# Not working

source("models/RACE/DAS/DAS.R")

dasNormal <- list(
  type="RACE",
  p_types=c("mean","sd","tau","t0"),
  # Transform to natural scale
  Ntransform=function(x) {
    exp(x)
  },
  # p_vector transform 
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
  # Random function for racing accumulators
  rfun=function(lR,pars) rDAS(lR,pars,rtype="normal"),
  # Density function (PDF) of choice race runner
  dfun=dnorm,
  # Density function (PDF) of choice race runner
  pfun=pnorm,
  # Density function (PDF) winner convolved with exponential
  dDAS=dDAS,
  # Race likelihood looping dDAS over trials
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_DAS(p_vector=p_vector, dadm = dadm, min_ll = min_ll,
                       pnames=c("mean","sd")) # must me in mean, sd order
)

