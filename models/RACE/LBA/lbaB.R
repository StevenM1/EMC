# lba_B parameterization 

source("models/RACE/LBA/lba.R")

lbaB <- list(
  type="RACE",
  # p_vector transform, sets sv as a scaling parameter
  p_types=c("v","sv","B","A","t0"),
  transform = function(p) p,
  # Transform to natural scale
  Ntransform=function(x) {
    get_p_types <- function(nams)
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams != "v"] <- exp(x[nams != "v"])
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams != "v"] <- exp(x[,nams != "v"])
    }
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    if (!is.matrix(pars)) {
      pars <- c(pars,b=pars["B"] + pars["A"])
      attr(pars,"ok") <- (pars["t0"] > .05) & ((pars["A"] > 1e-6) | pars[,"A"] == 0)
    } else {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
    }
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
    log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
)

