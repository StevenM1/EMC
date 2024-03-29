# Normal, natural mu and log(sigma), increasing threshold (first natural,
# others on log scale) parameterization

source("models/SDT/Probit/probit.R")

probit <- list(
  type="SDT", # Discrete choice based on continuous latent, no RT
  p_types=c("mean","sd","threshold"),
  # Transform to natural scale
  Ntransform=function(x) {
    is_sd <- grepl("sd",dimnames(x)[[2]]) 
    x[,is_sd] <- exp(x[,is_sd]) 
    x
  },
  # p_vector transform
  transform = function(x) {
    if (!is.matrix(x)) {
      increasing <- grepl("threshold",names(x)) & grepl(":lR",names(x)) | grepl("threshold_lR",names(x)) 
      x[increasing] <- exp(x[increasing])
      x
    } else {
      increasing <- grepl("threshold",dimnames(x)[[2]]) & grepl(":lR",dimnames(x)[[2]]) | 
        grepl("threshold_lR",dimnames(x)[[2]]) 
      x[,increasing] <- exp(x[,increasing])
      x
    }
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    pars
  },
  # Random function for discrete choices
  rfun=function(lR,pars) rPROBIT(lR,pars),
  # probability of choice between lower and upper thresholds (lt & ut)
  pfun=function(lt,ut,pars) pPROBIT(lt,ut,pars),
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood, lb is lower bound threshold for first response 
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)) 
    log_likelihood_sdt(p_vector=p_vector, dadm = dadm, min_ll = min_ll, lb=-Inf)
)

