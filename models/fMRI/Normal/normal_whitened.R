# Normal, mu and log(sigma) parameterization

source("models/fMRI/Normal/norm.R")

normal <- list(
  type="MRI",
  p_types=c("sd", "rho"), #This is a bit hacky for now
  # Transform to natural scale
  Ntransform=function(x) {
    
    get_p_types <- function(nams) 
      unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    
    nams <- get_p_types(dimnames(x)[[2]])
    x[,nams == "sd"] <- exp(x[,nams == "sd"])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm=NULL) 
  {
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rNORMAL(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dNORMAL(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pNORMAL(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
    y <- dadm[,!colnames(dadm) %in% c("subjects", "trials")]
    X <- attr(dadm, "model")$design_matrix[[dadm$subjects[1]]]
    order <- attr(dadm, "model")$toeplitz[[dadm$subjects[1]]]
    is_sd <- grepl("sd", names(p_vector))
    is_rho <- grepl("rho", names(p_vector))
    betas <- p_vector[!is_sd & !is_rho]
    sigma <- exp(p_vector[is_sd])
    rho <- p_vector[is_rho]
    V <- rho ** order * (sigma**2)
    return(max(mvtnorm::dmvnorm(y, mean=X %*% betas, sigma = V, log = T, checkSymmetry = F), min_ll*length(y)))
  
  }
)

