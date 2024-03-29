prior_samples_alpha <- function(theta_mu,theta_var,n=1e3)
  # samples from prior for alpha implied by hyper model 
  # (mixture over hyper posterior samples)
{
  t(sapply(sample.int(dim(theta_mu)[2],n,TRUE),function(x){
    mvtnorm::rmvnorm(1,theta_mu[x,],theta_var[,,x])
  }))
}


prior_samples <- function(samps,type=c("mu","variance","covariance","correlation","sigma")[1],n=1e4)
  # Samples from prior for standard at hyper levels for three types, 
  # or for single at subject level (in which case type="mu")
{
  if (type=="mu") {
    out <- mvtnorm::rmvnorm(n, mean = samps$prior$theta_mu_mean,
                      sigma = samps$prior$theta_mu_var)
    dimnames(out)[[2]] <- samps$par_names
    return(out) 
   } else {
    var <- array(NA_real_, dim = c(samps$n_pars, samps$n_pars, n),
                 dimnames = list(samps$par_names, samps$par_names, NULL))
    # Can't specify n in riwish?
    hyper <- attributes(samps)
    for(i in 1:n){
      a_half <- 1 / rgamma(n = samps$n_pars,shape = 1/2,
                           rate = 1/(hyper$A_half^2))
      var[,,i] <- riwish(hyper$v_half + samps$n_pars - 1, 2 * hyper$v_half * diag(1 / a_half))
    }
    if (type=="sigma") return(var)
    if (type=="variance") return(t(apply(var,3,diag)))
    if (type=="correlation")
      var <- array(apply(var,3,cov2cor),dim=dim(var),dimnames=dimnames(var))
    lt <- lower.tri(var[,,1])
    return(t(apply(var,3,function(x){x[lt]})))
  }
}


get_prior_samples <- function(samples,selection,filter,thin,subfilter,n_prior) 
  # get matrix of prior samples for different parameter types  
{
  if (class(samples)=="pmwgs") samps <- samples else samps <- samples[[1]]
  if (selection=="alpha") {
    if (samps$source !="samplers/pmwg/variants/single.R") {
      if (class(samples)=="pmwgs") {
        theta_mu <- as_Mcmc(samples,selection="mu",filter=filter,
                         thin=thin,subfilter=subfilter)
        theta_var <- get_sigma(samples,filter=filter,
                         thin=thin,subfilter=subfilter)
      } else {
        theta_mu <- do.call(rbind,as_mcmc.list(samples,selection="mu",
          filter=filter,thin=thin,subfilter=subfilter))
        theta_var <- get_sigma_chains(samples,filter=filter,
                         thin=thin,subfilter=subfilter) 
      }
      psamples <- prior_samples_alpha(theta_mu,theta_var,n_prior)
      colnames(psamples) <- colnames(theta_mu)
      return(psamples)
    } else return(prior_samples(samps,"mu",n_prior)) 
  } else {
    if (samps$source!="samplers/pmwg/variants/standard.R") {
      warnings("Prior plot not yet implemented for this type")
      return(NULL)  
    } else return(prior_samples(samps,selection,n_prior))
  }
}

