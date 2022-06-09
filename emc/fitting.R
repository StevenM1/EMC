#### Fitting automation
require(parallel)
require(abind)

run_stages <- function(sampler,iter=c(300,0,0),
                       verbose=FALSE,verbose_run_stage=FALSE,
                      max_adapt_trys=2,particles=NA,particle_factor=100, p_accept= NULL, n_cores=1,
                      epsilon = NULL, start_mu = NULL, start_var = NULL, mix=NULL,  
                      pdist_update_n=50,min_unique=200,epsilon_upper_bound=15,
                      eff_var = NULL, eff_mu = NULL) 
  # Initializes (if needed) then runs burn, adapt and sample if iter is not 
  # NA where iter[1] = burn, iter[2] = adapt, iter[3] = sample
  # Adapt stage is run repeatedly up to max_adapt_trys times. 
  # Number of particles is set at particle_factor*sqrt(number of parameters) if
  # particles is NA
{
  
  if (is.na(particles)) 
      particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
      sampler <- init(sampler, n_cores = n_cores, particles = particles, 
        epsilon = epsilon, start_mu = start_mu, start_var = start_var) 
  }
  if (all(iter==0)) return(sampler)
  if ( iter[1] != 0 ) {
    if (verbose) message("Running burn stage")
      sampler <- run_stage(sampler, stage = "burn",iter = iter[1], particles = particles, 
        n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
        min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (all(iter[2:3]==0)) return(sampler)
  }
  if (iter[2] != 0) {
    if (verbose) message("Running adapt stage")
    sampler <- run_stage(sampler, stage = "adapt",iter = iter[2], particles = particles, 
                         n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
                         min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (iter[3]==0) return(sampler)
  }
  if (verbose) message("Running sample stage")
  sampler <- run_stage(sampler, stage = "sample", iter = iter[3], epsilon = epsilon,
        pdist_update_n=pdist_update_n,particles = particles, n_cores = n_cores, 
        pstar = p_accept, verbose = verbose_run_stage, mix=mix, eff_mu = eff_mu,
        eff_var = eff_var,
        min_unique=min_unique,epsilon_upper_bound=epsilon_upper_bound)  
  sampler
}



# iter=c(500,0,0);
#   verbose=TRUE;verbose_run_stage=FALSE;mix=NULL;
#   max_adapt_trys=10;particles=NA;particle_factor=100; p_accept= 0.7;
#   cores_per_chain=1;cores_for_chains=NULL;min_unique=200;epsilon_upper_bound=2;
#   epsilon = NULL; start_mu = NULL; start_var = NULL;pdist_update_n=50
# verbose=TRUE; iter=c(3,0,0)
run_chains <- function(samplers,iter=c(300,0,0),
  verbose=TRUE,verbose_run_stage=FALSE,mix=NULL,
  max_adapt_trys=10,particles=NA,particle_factor=100, p_accept= 0.7, 
  cores_per_chain=1,cores_for_chains=NULL,min_unique=200,epsilon_upper_bound=15,
  epsilon = NULL, start_mu = NULL, start_var = NULL,pdist_update_n=50) 
  # applies run stages over chains preserving list attributes
{
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  run_try <- 0
  repeat {
    samplers_new <- mclapply(samplers,run_stages,iter=iter,
      verbose=verbose,verbose_run_stage=verbose_run_stage,
      max_adapt_trys=max_adapt_trys,particles=particles,particle_factor=particle_factor, 
      p_accept=p_accept, min_unique=min_unique, mix=mix,
      epsilon = epsilon, start_mu = start_mu, start_var = start_var,
      pdist_update_n=pdist_update_n,epsilon_upper_bound=epsilon_upper_bound,
      n_cores=cores_per_chain, mc.cores = cores_for_chains)
    if (class(samplers_new)=="try-error" || 
      any(lapply(samplers_new,class)=="try-error")) {
      save(samplers,iter,particles,particle_factor,p_accept,pdist_update_n,
           epsilon_upper_bound,min_unique,file="fail_run_stage.RData")
      run_try <- run_try + 1
      if (verbose) message("run_stage try error", run_try)
      if (run_try > 10)
        stop("Fail after 10 run_stage try errors, see saved samplers in fail_run_stage.RData")
    } else break
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  samplers
}


run_gd <- function(samplers,iter=NA,max_trys=100,verbose=FALSE,burn=TRUE,
                   max_gd=1.1,thorough=TRUE, natural=FALSE,
                   epsilon = NULL,particle_update = 5, verbose_run_stage = F,
                   particles=NA,particle_factor=100, p_accept=NULL,
                   cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                   min_es=NULL,min_iter=NULL,max_iter=NULL) 
  # Repeatedly runs burn or sample to get subject average multivariate 
  # gelman.diag of alpha samples < max_gd (if !through) or if all univariate
  # psrf for every subject and parameter and multivariate < max_gd and
  # if min_es specified at least that effective size in worst case and if 
  # min_iter specified at least that many iterations. If max_iter specified
  # pulls out after max_iter.
  # Trys adding iter (if NA 1/3 initial length) of the length and keeps all or 
  # extra minus first n_remove (initially iter, then grows by iter/2 on each
  # expansion based on mean gd of alphas (mpsrf alone or with psrf depending on thorough)
  # until success or max_trys. Verbose prints out progress after each try
  # Cores used = cores_per_chain*cores_for_chains, latter set to number of chains by default
{
  
  enough_samples <- function(samplers,min_es,min_iter,max_iter,filter="burn") {
    if (!is.null(max_iter) && (sum(samplers[[1]]$samples$stage==filter) >= max_iter)) return(TRUE)
    if (is.null(min_iter)) ok <- TRUE else
      ok <- (sum(samplers[[1]]$samples$stage==filter)) > min_iter
    if (!is.null(min_es)) {
      es_min <- min(es_pmwg(as_mcmc.list(samplers,selection="alpha",filter=filter)))
      ok <- ok & (es_min > min_es)
      attr(ok,"es") <- es_min
    }
    ok
  }
  
  if (thorough) {
    max_name <- ", Max alpha mpsrf/psrf = "
    message("Exit on max alpha psrf/mpsrf < ",max_gd)
  } else {
    max_name <- ", Max alpha msrf = "
    message("Exit on max alpha mpsrf < ",max_gd)
  }
  n_chains <- length(samplers)
  if (is.na(cores_for_chains)) cores_for_chains <- n_chains
  if (burn) stage <- "burn" else stage <- "sample"
  n <- chain_n(samplers)
  if (burn) {
    if (sum(n[,-1])>0) 
      stop("Can only use burn=TRUE if samplers only have burn samples")
    idxs <- n[,1]
  } else {
    if (sum(n[,3])==0) 
      stop("Can only use burn=FALSE if there are sample stage samples")
    idxs <- n[,3]
  }
  if (!all(idxs[1]==idxs[-1]))
    stop("Must have same number of iterations for each chain")
  if (is.na(iter)) iter <- round(idxs[1]/3)
  n_remove <- iter
  # Iterate until criterion
  trys <- 0
  shorten <- FALSE
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  repeat {
    run_try <- 0
    repeat {
      particles <- round(particle_factor*sqrt(length(samplers[[1]]$par_names)))
      samplers_new <- mclapply(samplers,run_stage,stage=stage,iter=iter, mix=mix,
                               epsilon = epsilon, particles=particles,pstar=p_accept,verbose=verbose_run_stage,
                               n_cores=cores_per_chain,mc.cores=cores_for_chains)
      if (class(samplers_new)=="try-error" || 
          any(lapply(samplers_new,class)=="try-error")) {
        save(samplers,stage,iter,particles,file="fail_run_stage.RData")
        run_try <- run_try + 1
        if (verbose) message("run_stage try error", run_try)
        if (run_try > 10)
          stop("Fail after 10 run_stage try errors, see saved samplers in fail_run_stage.RData")
      } else {
        gd <- gd_pmwg(as_mcmc.list(samplers_new,filter=stage),!thorough,FALSE,
                      filter=stage,natural=natural)
        if (is.finite(gd[[1]])) break else
          if (verbose) message("gelman diag try error", run_try)
      }
    }
    samplers <- samplers_new
    trys <- trys+1
    if (shorten) {
      samplers_short <- lapply(samplers,remove_iterations,select=n_remove,filter=stage)
      gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter=stage),!thorough,FALSE,
                          filter=stage,natural=natural)
      if (mean(gd_short) < mean(gd)) {
        gd <- gd_short
        samplers <- samplers_short
        n_remove <- iter
      } else {
        n_remove <- round(n_remove + iter/2)
      }
    }
    enough <- enough_samples(samplers,min_es,min_iter,max_iter,filter=filter)
    if (is.null(attr(enough,"es"))) es_message <- NULL else
      es_message <- paste(", Effective samples =",round(attr(enough,"es")))
    gd_add <- (gd[,3] - max_gd)
    gd_add[gd_add > 0] <- pmin(10, sqrt(gd_add[gd_add > 0])*particle_update)
    gd_add[gd_add < 0] <- -5
    particle_factor <- pmin(pmax(25, particle_factor + gd_add), 100)
    ok_gd <- all(gd < max_gd)
    shorten <- !ok_gd
    if (trys > max_trys || (ok_gd & enough)) {
      if (verbose) {
        message("Final multivariate gelman.diag per participant")
        message("\nIterations = ",samplers[[1]]$samples$idx,", Mean mpsrf= ",
                round(mean(gd),3),max_name,round(max(gd),3))
      }
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
      return(samplers)
    }
    if (verbose) {
      chain_n(samplers)[,stage][1]
      message(trys,": Iterations (",stage,") = ",chain_n(samplers)[,stage][1],", Mean mpsrf= ",
              round(mean(gd),3),max_name,round(max(gd),3),es_message)
    }
  }
}

# burn=TRUE;ndiscard=200;nstart=300;nadapt=1000;
#     discard_start=TRUE;start_particles=NA;
#     start_mu = NULL; start_var = NULL;
#     start_mix=NULL;single_start_mix=c(.5,.5);
#     verbose=TRUE;verbose_run_stage=FALSE;
#     max_gd_trys=100;max_gd=1.1;thorough=TRUE;
#     min_es=NULL;min_iter=NULL;max_iter=NULL
#     epsilon = NULL; epsilon_upper_bound=2;
#     particles=NA;particle_factor=100; sample_particle_factor=100; p_accept=0.7;
#     pdist_update_n=50;min_unique=200; mix=NULL;
#     cores_per_chain=1;cores_for_chains=NULL
# 
# burn=TRUE;ndiscard=0;nstart=500;cores_per_chain=10;max_gd_trys=0
# 
# burn=TRUE;ndiscard=0;nstart=100;max_iter=2000; cores_per_chain=19;cores_for_chains=3
# cores_per_chain=1;cores_for_chains=1


auto_burn <- function(samplers,ndiscard=200,nstart=300,nadapt=1000,
    discard_start=TRUE,start_particles=NA,
    start_mu = NULL, start_var = NULL,
    start_mix=NULL,single_start_mix=c(.5,.5),
    verbose=TRUE,verbose_run_stage=FALSE,
    max_gd_trys=100,max_gd=1.1,
    thorough=TRUE,natural=FALSE,
    min_es=NULL,min_iter=NULL,max_iter=NULL,
    epsilon = NULL, epsilon_upper_bound=2, 
    particles=NA,particle_factor=100, sample_particle_factor=100, p_accept=0.7,
    pdist_update_n=50,min_unique=200, mix=NULL,
    cores_per_chain=1,cores_for_chains=NULL)
  # Takes a pmwgs chains list, initializes it (see run_stages), if !burn adapts
  # and runs burn or sample until gd criterion satisfied (see run_gd for details)
{
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  if (is.null(samplers[[1]]$samples$idx) || samplers[[1]]$samples$idx==1) {
    if (ndiscard==0) discard_start <- FALSE
    if (!discard_start) discard_message <- NULL else
      discard_message <- paste("(first ",ndiscard," then removed)")
    message("Getting initial ",ndiscard + nstart," samples ",discard_message)
    run_try <- 0
    if (samplers[[1]]$source=="samplers/pmwg/variants/single.R") start_mix <- single_start_mix
    run_try <- 0
    repeat {
      samplers_new <- mclapply(samplers,run_stages,iter=c(nstart + ndiscard,0,0),
                               n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                               particles=start_particles,particle_factor=start_particle_factor,epsilon=epsilon,
                               verbose=FALSE,verbose_run_stage=verbose_run_stage,
                               mc.cores=cores_for_chains)
      if (class(samplers_new)=="try-error" || 
          any(lapply(samplers_new,class)=="try-error")) {
        save(samplers,nstart,particles,particle_factor,p_accept, file="fail_run_stage.RData")
        run_try <- run_try + 1
        if (verbose) message("burn run_stage try error", run_try)
        if (run_try > 10)
          stop("Fail after 10 burn run_stage try errors, see saved samplers in fail_run_stage.RData")
      } else break
    }
    samplers <- samplers_new
    if (discard_start) samplers <- lapply(samplers,remove_iterations,select=nstart+1)
    message("Finished initial run")
    attr(samplers,"data_list") <- data_list
    attr(samplers,"design_list") <- design_list
    attr(samplers,"model_list") <- model_list
  } 
  if (max_gd_trys==0) return(samplers)
  gd <- gd_pmwg(as_mcmc.list(samplers,filter="burn"),!thorough,FALSE,
                filter="burn",natural=natural)
  gd_add <- (gd[,3] - max_gd)
  gd_add[gd_add > 0] <- pmin(10, sqrt(gd_add[gd_add > 0])*particle_update)
  gd_add[gd_add < 0] <- -5
  particle_factor <- pmin(pmax(25, start_particle_factor + gd_add), 100)
  message("Beginning iterations to achieve Rhat < ",max_gd)
  run_gd(samplers, max_trys=max_gd_trys,verbose=verbose,
         max_gd=max_gd,thorough=thorough, natural=natural, p_accept = p_accept,
         epsilon=epsilon, particles=start_particles,particle_factor=particle_factor,
         particle_update = particle_update, min_es=min_es,min_iter=min_iter, 
         max_iter=max_iter, mix=mix, iter = step_size, verbose_run_stage = verbose_run_stage,
         cores_per_chain=cores_per_chain,cores_for_chains=cores_for_chains)
}

auto_adapt <- function(samplers,max_trys=25,verbose=FALSE, epsilon = NULL, particles=NA,particle_factor=40, p_accept=.7,
                       cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                       n_cores_conditional = 1, min_es=NULL,min_unique = 200, 
                       step_size = 25, thin = NULL,
                       verbose_run_stage = F){
    if(verbose) message("Running adapt stage")
    source(samplers[[1]]$source)
    if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
    data_list <- attr(samplers,"data_list")  
    design_list <- attr(samplers,"design_list")
    model_list <- attr(samplers,"model_list")
    trys <- 0
    samplers_new <- samplers
    repeat {
      trys <- trys + 1
      samplers_new <- mclapply(samplers_new,run_stages,iter=c(0,step_size,0),
                               n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                               particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                               verbose=FALSE,verbose_run_stage=verbose_run_stage,
                               mc.cores=cores_for_chains)
      test_samples <- lapply(samplers_new, extract_samples, stage = "adapt", thin = thin, 50*trys, thin_eff_only = F)
      keys <- unique(unlist(lapply(test_samples, names)))
      test_samples <- setNames(do.call(mapply, c(abind, lapply(test_samples, '[', keys))), keys)
      test_samples$iteration <- sum(test_samples$iteration)
      # Only need information like n_pars & n_subjects, so only need to pass the first chain
      adapted <- test_adapted(samplers_new[[1]], test_samples, min_unique, n_cores_conditional, verbose_run_stage)
      if(trys > max_trys | adapted) break
    }
    samplers <- samplers_new
    attr(samplers,"data_list") <- data_list
    attr(samplers,"design_list") <- design_list
    attr(samplers,"model_list") <- model_list
    return(samplers)
}


auto_sample <- function(samplers,iter=NA,verbose=FALSE,
                       epsilon = NULL, particles=NA,particle_factor=25, p_accept=.7,
                       cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                       n_cores_conditional = 1, min_es=NULL, step_size = 50, thin = NULL,
                       verbose_run_stage = F){
  if(verbose) message("Running sample stage")
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  samplers_new <- samplers
  n_steps <- ceiling(iter/step_size)
  for(step in 1:n_steps){
    if(step == n_steps){
      step_size <- iter %% step_size
    } 
    test_samples <- lapply(samplers_new, extract_samples, stage = c("adapt", "sample"), thin = thin, 50*trys, thin_eff_only = F)
    keys <- unique(unlist(lapply(test_samples, names)))
    test_samples <- setNames(do.call(mapply, c(abind, lapply(test_samples, '[', keys))), keys)
    test_samples$iteration <- sum(test_samples$iteration)
    conditionals=mclapply(X = 1:samplers_new[[1]]$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples, samplers_new[[1]]$n_pars, mc.cores = n_cores_conditional)
    conditionals <- array(unlist(conditionals), dim = c(samplers_new[[1]]$n_pars, samplers_new[[1]]$n_pars + 1, samplers_new[[1]]$n_subjects))
    eff_mu <- conditionals[,1,] #First column is the means
    eff_var <- conditionals[,2:(samplers_new[[1]]$n_pars+1),] #Other columns are the variances
    samplers_new <- mclapply(samplers_new,run_stages,iter=c(0,0,step_size),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage, eff_mu = eff_mu,
                             eff_var = eff_var,
                             mc.cores=cores_for_chains)
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  return(samplers)
}

test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1, verbose = F){
  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  n_pars <- sampler$n_pars
  if (all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      mclapply(X = 1:sampler$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples, 
               n_pars, mc.cores = n_cores_conditional)
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        message("A problem was encountered creating proposal distribution")
        message("Increasing required unique values and continuing adaptation")
      }
      return(FALSE)
    }
    else {
      if(verbose) message("Successfully adapted after ", test_samples$iteration, "iterations - stopping adaptation")
      return(TRUE)
    }
  } else{
    return(FALSE) # Not enough unique particles found
  }
}
