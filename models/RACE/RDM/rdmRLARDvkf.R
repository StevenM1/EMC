# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)
# Added v0 parameter for initial rates in RL version.

source("models/RACE/RDM/rdm.R")

rdmRLARDvkf <- list(
  type="RACE",
  p_types=c("v0","v","B","A","t0","s", "wd", "ws",
            "alpha", "q0", "w0", "volatility0"),
  
  Ntransform=function(x) {
    # Transform to natural scale
    probit_scale <- c('alpha', 'q0')  
    # actually, q0 can take more flexible norms outside of the binomial reward function?
    x[,!dimnames(x)[[2]] %in% probit_scale] <- exp(x[,!dimnames(x)[[2]] %in% probit_scale])
    x[,dimnames(x)[[2]] %in% probit_scale] <- pnorm(x[,dimnames(x)[[2]] %in% probit_scale])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    parsList <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                         levels(dadm$subjects))
    for (i in levels(dadm$subjects))
      parsList[[i]] <- update_pars(i,pars,dadm)
    pars <- do.call(rbind,parsList)
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) 
    pars
  },
  # p_vector transform 
  transform = function(x) x,
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

