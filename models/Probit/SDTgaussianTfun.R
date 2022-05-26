# Normal, natural mu and log(sigma), increasing threshold (first natural,
# others on log scale) parameterization

source("models/Probit/probit.R")

probitTfun <- list(
  type="SDT", # Discrete choice based on continuous latent, no RT
  p_types=c("mean","sd","threshold","slope"),
  # Transform to natural scale
  Ntransform=function(x) {

    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      x[nams == "sd"] <- exp(x[nams == "sd"])
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      x[,nams == "sd"] <- exp(x[,nams == "sd"])
    }
    x
  },
  # mapped parameter transform
  Mtransform = function(pars) 
    # transform parameters except v back to real line 
    # pars is a matrix output by map_p_vector  
  {
    probit$Ntransform(pars)
  },
  # p_vector transform
  transform = function(x) {
    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    nams <- get_p_types(names(x))
    x[nams == "threshold"][-1] <- 
      diff(c(x[nams == "threshold"][1] + x["slope"]*(x[nams == "threshold"][-1]),1e100))
    x
  },
  # Random function for discrete choices
  rfun=function(lR,pars) rPROBIT(lR,pars),
  # probability of choice given by rt (actually response)
  pfun=function(rt,pars) pPROBIT(rt,pars),
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood, lb is lower bound threshold for first response 
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_sdt(p_vector=p_vector, dadm = dadm, min_ll = min_ll, lb=-Inf)
)

