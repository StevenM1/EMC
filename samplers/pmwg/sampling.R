
source("samplers/pmwg/messaging.R")

pmwgs <- function(dadm, pars = NULL, ll_func = NULL, prior = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(ll_func)) ll_func <- dadm$ll_func
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list
  
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  samples <- variant_funs$sample_store(dadm, pars, ...)
  sampler <- list(
    data = dadm_list,
    par_names = c(pars, names(dadm$subject_covariates)),
    subjects = subjects,
    n_pars = length(pars) + length(dadm$subject_covariates),
    subject_covariates = dadm$subject_covariates,
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- variant_funs$add_info(sampler, prior, ...)
  return(sampler)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  startpoints <- variant_funs$get_startpoints(pmwgs, start_mu, start_var)
  if(n_cores > 1){
    proposals <- mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals, parameters = startpoints, n_particles = particles, pmwgs = pmwgs, mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals, parameters = startpoints, n_particles = particles, pmwgs = pmwgs)
  }
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
  
  # Sample the mixture variables' initial values.
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars, verbose), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                epsilon = epsilon, j = 1, n_pars = pmwgs$n_pars)
  pmwgs$init <- TRUE
  return(pmwgs)
}


start_proposals <- function(s, parameters, n_particles, pmwgs){
  #Draw the first start point
  group_pars <- variant_funs$get_group_level(parameters, s)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,dadm = pmwgs$data[[which(pmwgs$subjects == s)]])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL, 
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0), 
                          likelihood_func = NULL, epsilon = NULL, subjects,
                          components, prev_ll, shared_ll_idx) 
{
  unq_components <- unique(components)
  group_pars <- variant_funs$get_group_level(parameters, s)
  proposal_out <- numeric(length(group_pars$mu))
  ll <- numeric(length(unq_components))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- parameters$alpha[,s]
  eff_mu_sub <- eff_mu[,s]
  if(is.null(eff_mu_sub)){
    eff_mu_sub <- subj_mu
  } 
  num_particles <- num_particles[s]
  for(i in unq_components){
    if(is.null(eff_mu[,s])){
      eff_var_curr <- eff_var * epsilon[s,i]^2
    } else{
      eff_var_curr <- eff_var
    }
    var_subj <- group_var *  epsilon[s,i]^2
    idx <- components == i
    # Draw new proposals for this component
    particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
    cumuNumbers <- cumsum(particle_numbers) + 1 # Include the particle from b4
    pop_particles <- particle_draws(particle_numbers[1], group_mu[idx], group_var[idx, idx])
    ind_particles <- particle_draws(particle_numbers[2], subj_mu[idx], var_subj[idx, idx])
    if(mix_proportion[3] == 0){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mu_sub[idx], eff_var_curr[idx, idx , s])
    }
    # Rejoin new proposals with current MCMC values for other components
    proposals <- matrix(rep(subj_mu, num_particles + 1), nrow = num_particles + 1, byrow = T)
    colnames(proposals) <- names(subj_mu)
    proposals[2:(num_particles+1),idx] <- rbind(pop_particles, ind_particles, eff_particles)
    
    # should there be any blocking within a component, include those parameters as well.
    lw_calc_idx <- shared_ll_idx == shared_ll_idx[idx][1]
    # Calculate likelihoods and other IS quantities
    lw <- apply(proposals[,lw_calc_idx], 1,  likelihood_func, dadm = data[[which(subjects == s)]])
    lw_total <- lw + prev_ll[s] - lw[1] # Bit inefficient but safes code, makes sure lls from other components are included
    lp <- mvtnorm::dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
    prop_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu, sigma = var_subj)
    if (mix_proportion[3] == 0) {
      eff_density <- 0
    }
    else {
      eff_density <- mvtnorm::dmvnorm(x = proposals, mean = eff_mu_sub, sigma = eff_var_curr[,,s])
    }
    lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density) + (mix_proportion[3] * eff_density))
    # Calculate weights and center
    l <- lw_total + lp - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = num_particles+1, size = 1, prob = weights)
    origin <- min(which(idx_ll <= cumuNumbers))
    ll[i] <- lw_total[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
  } # Note that origin only contains last components origin, just used by devs anyway
  return(list(proposal = proposal_out, ll = mean(ll), origin = origin))
}



run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      verbose = TRUE,
                      force_prev_epsilon = TRUE,
                      n_cores = 1,
                      epsilon = NULL,
                      pstar = NULL,
                      mix = NULL,
                      preburn = 50,
                      pdist_update_n = 50,
                      epsilon_upper_bound = 15,
                      n_cores_conditional = 1,
                      thin = NULL,
                      thin_eff_only = FALSE) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  components <- attr(pmwgs$data, "components")
  shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")
  # Display stage to screen
  if(verbose){
    msgs <- list(
      burn = "Phase 1: Burn in\n",
      adapt = "Phase 2: Adaptation\n",
      sample = "Phase 3: Sampling\n"
    )
    cat(msgs[[stage]])
  }
  
  alphaStar=-qnorm(pstar/2) #Idk about this one
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  
  epsilon <- fix_epsilon(pmwgs, epsilon, force_prev_epsilon, components)
  if(length(particles == 1)){
    particles <- rep(particles, pmwgs$n_subjects)
  }
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  eff_mu <- attr(pmwgs, "eff_mu")
  eff_var <- attr(pmwgs, "eff_var")
  mix <- set_mix(eff_var, mix, verbose)
  if (verbose) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  unq_components <- unique(components)
  # Main iteration loop
  for (i in 1:iter) {
    if (verbose) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i
    
    # Gibbs step
    pars <- variant_funs$gibbs_step(pmwgs, rbind(pmwgs$samples$alpha[,,j-1], pmwgs$subject_covariates))
    # Particle step
    proposals=mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu, 
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, components, 
                       prev_ll = pmwgs$samples$subj_ll[,j-1], shared_ll_idx, mc.cores =n_cores)
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))
    
    #Fill samples
    pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = pars,
                                  proposals = proposals, epsilon = rowMeans(epsilon), j = j, n_pars = pmwgs$n_pars)
    
    # Update epsilon
    if(!is.null(pstar)){
      if(j > n0){
        for(component in unq_components){
          idx <- components == component
          acc <-  pmwgs$samples$alpha[max(which(idx)),,j] != pmwgs$samples$alpha[max(which(idx)),,(j-1)]
          epsilon[,component] <-pmin(update.epsilon(epsilon[,component]^2, acc, pstar, j, sum(idx), alphaStar), epsilon_upper_bound)
        }
      }
    }

  }
  if (verbose) close(pb)
  attr(pmwgs, "epsilon") <- epsilon
  if(!is.null(thin) & !thin_eff_only){
    pmwgs <- thin_objects(pmwgs, thin, i)
  }
  return(pmwgs)
}

test_sampler_adapted <- function(pmwgs, n_unique, i, n_cores_conditional, conditionals_func, verbose = T, thin, thin_eff_only) {
  n_pars <- length(pmwgs$par_names)
  if (i < n_unique) {
    return("continue")
  }
  test_samples <- extract_samples(pmwgs, stage = "adapt", thin, i, thin_eff_only)
  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  if (all(n_unique_sub > n_unique)) {
    if(verbose){
      message("Enough unique values detected: ", n_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      mclapply(X = 1:pmwgs$n_subjects,FUN = conditionals_func,samples = test_samples, 
               n_pars, mc.cores = n_cores_conditional)
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        warning("A problem was encountered creating proposal distribution")
        warning("Increasing required unique values and continuing adaptation")
      }
      return("increase")
    }
    else {
      if(verbose) message("Successfully adapted after ", i, "iterations - stopping early")
      return("success")
    }
  }
  return("continue")
}

particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  } else if(length(mu) == 1){
    return(as.matrix(rnorm(n, mu, sqrt(covar))))
  }
  return(mvtnorm::rmvnorm(n, mu, covar))
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c=((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta=log(sqrt(epsilon2))
  Theta=Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}

set_epsilon <- function(n_pars, verbose = T) {
  if (n_pars > 15) {
    epsilon <- 0.1
  } else if (n_pars > 10) {
    epsilon <- 0.3
  } else {
    epsilon <- 0.5
  }
  if(verbose) message(sprintf("Epsilon has been set to %.1f based on number of parameters",epsilon))
  return(epsilon)
}

set_mix <- function(eff_var, mix, verbose) {
  if (!is.null(eff_var)) {
    mix <- c(0.15, 0.15, 0.7)
  } else {
    mix <- c(0.5, 0.5, 0.0)
  }
  if(verbose) message(sprintf("mix has been set to c(%s) based on the stage being run",  paste(mix, collapse = ", ")))
  return(mix)
}

numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  return(numbers)
}

extract_samples <- function(sampler, stage = c("adapt", "sample"), thin, i, thin_eff_only) {
  samples <- sampler$samples
  stage_filter <- which(samples$stage %in% stage)
  if(!is.null(thin)){
    if(thin_eff_only){
      # Assume there's no thinning done on the actual samples, just for the eff
      thin <- seq(thin, samples$idx, by = thin)
      full_filter <- intersect(thin, stage_filter)
    } else{
      # There has been thinning on prev samples, do it for current ones as well
      old_idx <- 1:(samples$idx - i - 1)
      old_filter <- intersect(old_idx, stage_filter)
      if(i > thin){ # make sure we have enough
        new_filter <- samples$idx - i + seq(thin,i,by=thin)
        full_filter <- c(old_filter, new_filter)
      } else{
        full_filter <- old_filter
      }
    }
  } else{
    # No thinning
    sampled_filter <- which(seq_along(samples$stage) <= samples$idx)
    full_filter <- intersect(stage_filter, sampled_filter)
  }
  out <- variant_funs$filtered_samples(sampler, full_filter)
  return(out)
}

sample_store_base <- function(data, par_names, iters = 1, stage = "init") {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  samples <- list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

fill_samples_base <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  # Fill samples both group level and random effects
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
  return(samples)
}

fill_samples_RE <- function(samples, proposals, epsilon, j = 1, n_pars){
  # Only for random effects, separated because group level sometimes differs.
  samples$alpha[, , j] <- proposals[1:n_pars,]
  samples$subj_ll[, j] <- proposals[n_pars + 1,]
  samples$origin[,j] <- proposals[n_pars + 2,]
  samples$idx <- j
  samples$epsilon[,j] <- epsilon
  return(samples)
}

extend_obj <- function(obj, n_extend){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  if(n_dimensions == 2){
    if(isSymmetric(round(obj, 5))) return(obj) #Don't extend priors and theta_mu_var_inv
  }  
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
  return(extended)
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. Big shout out to the rapply function
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  sampler$samples <- base::rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  return(sampler)
}

condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE) 
{
  if (missing(dependent.ind)) 
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given)) 
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind, 
                                                                          dependent.ind])))
  if (length(given.ind) == 0) 
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind, 
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind)) 
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma)) 
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08)) 
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}

trim_na <- function(sampler) {
  idx <- 1:sampler$samples$idx
  sampler$samples$stage <- sampler$samples$stage[idx]
  sampler$samples <- base::rapply(sampler$samples, f = function(x) filter_obj(x, idx), how = "replace")
  return(sampler)
}

fix_epsilon <- function(pmwgs, epsilon, force_prev_epsilon, components){
  if(is.null(epsilon) | force_prev_epsilon){
    epsilon <- attr(pmwgs, "epsilon")
    if(is.null(epsilon)){
      epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
    }
  }
  if(is.matrix(epsilon)){
    if(ncol(epsilon) != max(components)){
      epsilon <- matrix(rowMeans(epsilon), nrow = pmwgs$n_subjects, ncol = max(components))
    }
  } else if (length(epsilon) == 1 | is.vector(epsilon)){
    epsilon <- matrix(epsilon, nrow = pmwgs$n_subjects, ncol = max(components))
  }
  return(epsilon)
}

thin_objects <- function(sampler, thin, i){
  old_idx <- 1:(sampler$samples$idx - i - 1)
  new_idx <- sampler$samples$idx - i + seq(thin,i,by=thin)
  nmc_thin <- c(old_idx, new_idx)
  sampler$samples$stage <- sampler$samples$stage[nmc_thin]
  sampler$samples$idx <- length(sampler$samples$stage)
  sampler$samples <- base::rapply(sampler$samples, f = function(x) filter_obj(x, nmc_thin), how = "replace")
  return(sampler)
}