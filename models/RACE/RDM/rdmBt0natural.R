# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)

source("models/RACE/RDM/rdm.R")

rdmBt0natural <- list(
  type="RACE",
  p_types=c("v","B","A","t0","s"),
  # Transform to natural scale
  Ntransform=function(x) {
    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams != "t0"] <- exp(x[nams != "t0"]) 
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams != "t0"] <- exp(x[,nams != "t0"])
    }
    x
  },
  # p_vector transform 
  transform = function(x) x,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
    pars
  },
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

