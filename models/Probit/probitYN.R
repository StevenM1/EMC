# Normal, mu and log(sigma) parameterization

source("models/Probit/probit.R")

probit <- list(
  type="SDT", # Discrete choice based on continuous latent, no RT
  p_types=c("mean","sd","threshold"),
  # Transform to natural scale
  Ntransform=function(x) {

    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    if (!is.matrix(x)) {
      nams <- get_p_types(names(x))
      # x[nams %in% c("sd","threshold")] <- exp(x[nams %in% c("sd","threshold")])
      x[nams == "sd"] <- exp(x[nams == "sd"])
      # is_threshold <- nams == "threshold"
      # if (sum(is_threshold)>1) # transform all but intercept
      #   x[is_threshold][-1] <- exp(x[is_threshold][-1])
    } else {
      nams <- get_p_types(dimnames(x)[[2]])
      # x[,nams %in% c("sd","threshold")] <- exp(x[,nams %in% c("sd","threshold")])
      x[,nams == "sd"] <- exp(x[,nams == "sd"])
      # is_threshold <- nams == "threshold"
      # if (sum(is_threshold)>1) # transform all but intercept
      #   x[,is_threshold][,-1] <- exp(x[,is_threshold][,-1])
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
    x[nams == "threshold"][-1] <- exp(x[nams == "threshold"][-1])
    x
  },
  # Random function for discrete choices
  rfun=function(lR,pars) rPROBIT(lR,pars),
  # probability of choice given by rt
  pfun=function(rt,pars) pPROBIT(rt,pars),
  # Likelihood, lb is lower bound threshold for first response 
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_dc(p_vector=p_vector, dadm = dadm, min_ll = min_ll, lb=-Inf)
)

