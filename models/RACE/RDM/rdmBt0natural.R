# RDM_B parameterization with s=1 scaling (B = b-A done in rdm.R)

source("models/RACE/RDM/rdm.R")

rdmBt0natural <- list(
  type="RACE",
  p_types=c("v","B","A","t0"),
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
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except t0 back to real line 
    # pars is a matrix output by map_p_vector  
  {
    pars <- rdmBt0natural$Ntransform(pars)
    attr(pars,"ok") <- pars[,"t0"]>0
    pars

  },
  # p_vector transform scaling parameter by s=1 assumed in rdm.R
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
