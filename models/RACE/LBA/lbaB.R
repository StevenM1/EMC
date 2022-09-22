# lba_B parameterization 

source("models/RACE/LBA/lba.R")

lbaB <- list(
  type="RACE",
  p_types=c("v","sv","B","A","t0"),
  Ntransform=function(x) {
  # Transform to natural scale
   
    get_p_types <- function(nams)
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams != "v"] <- exp(x[nams != "v"])
      x <- c(x,b=x["B"] + x["A"])
      attr(x,"ok") <- (x["t0"] > .05) & ((x["A"] > 1e-6) | x[,"A"] == 0)
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams != "v"] <- exp(x[,nams != "v"])
      x <- cbind(x,b=x[,"B"] + x[,"A"])
      attr(x,"ok") <- (x[,"t0"] > .05) & ((x[,"A"] > 1e-6) | x[,"A"] == 0)
    }
    x
  },
  # p_vector transform, sets sv as a scaling parameter
  transform = function(p) p,
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) pars,
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

