# lba_B parameterization with sv=1 scaling

source("models/RACE/LBA/lba.R")

MlbaB <- list(
  type="RACE",
  p_types=c("v","sv","B","A","t0"),
  Ntransform=function(x) {
  # Transform to natural scale
    x[,dimnames(x)[[2]] != "v"] <- exp(x[,dimnames(x)[[2]] != "v"])
    x
  },
  # p_vector transform, sets sv as a scaling parameter
  transform = function(p) p,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
    pars
  },
  # Random function for racing accumulator
  rfun=function(lR,pars) rLBA(lR,pars,posdrift=TRUE),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

