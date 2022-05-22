#### pmwg object functions ----

thin_pmwg <- function(pmwg,thin=10) {
  # mcmc creation functions allow thinning, but for cases where you want to 
  # reduce the size, this will thin the sample component of a pmwg object  
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (thin >= pmwg$samples$idx) stop("Thin value to large\n")
  nmc_thin <- seq(thin,pmwg$samples$idx,by=thin)
  pmwg$samples$idx <- length(nmc_thin)
  pmwg$samples$alpha <- pmwg$samples$alpha[,,nmc_thin,drop=FALSE]
  pmwg$samples$theta_mu <- pmwg$samples$theta_mu[,nmc_thin,drop=FALSE]
  pmwg$samples$theta_var <- pmwg$samples$theta_var[,,nmc_thin,drop=FALSE]
  pmwg$samples$subj_ll <- pmwg$samples$subj_ll[,nmc_thin,drop=FALSE]
  pmwg$samples$a_half <- pmwg$samples$a_half[,nmc_thin,drop=FALSE]
  pmwg$samples$origin <- pmwg$samples$origin[,nmc_thin,drop=FALSE] 
  pmwg$samples$stage <- pmwg$samples$stage[nmc_thin]
  pmwg
}


get_stage <- function(pmwg,stage="burn") 
  # Returns object with only stage samples
{
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (!(stage %in% pmwg$samples$stage))
    stop("stage not available")
  ok <- pmwg$samples$stage == stage
  pmwg$samples$alpha <- pmwg$samples$alpha[,,ok,drop=FALSE]
  pmwg$samples$theta_mu <- pmwg$samples$theta_mu[,ok,drop=FALSE]
  pmwg$samples$theta_var <- pmwg$samples$theta_var[,,ok,drop=FALSE]
  pmwg$samples$subj_ll <- pmwg$samples$subj_ll[,ok,drop=FALSE]
  pmwg$samples$a_half <- pmwg$samples$a_half[,ok,drop=FALSE]
  pmwg$samples$origin <- pmwg$samples$origin[,ok,drop=FALSE] 
  pmwg$samples$stage <- pmwg$samples$stage[ok]
  pmwg$samples$idx <- sum(ok)
  pmwg
}


remove_iterations <- function(pmwg,select,remove=TRUE,last_select=FALSE,filter=NULL) 
  # Removes samples up to scalar remove or in remove vector
  # if remove=FALSE instead keeps these iterations
  # if last_select applies scalar select from last iteration
{
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (length(select)==1) {
    if (is.null(filter)) {
      if (!last_select)
        select <- 1:select else
        select <- (pmwg$samples$idx-select+1):pmwg$samples$idx
    } else {
      n <- table(pmwg$samples$stage)
      if (filter=="burn") {
        if (select>n["burn"]) stop("Removing more than available in burn")
        if (!last_select)
          select <- 1:select else
          select <- (n["burn"]-select+1):n["burn"]
      } else if (filter=="adapt") {
        if (select>n["adapt"]) stop("Removing more than available in adapt")
        if (!last_select)
          select <- 1:select else
          select <- (n["adapt"]-select+1):n["adapt"]
        select <- select + n["burn"]
      } else if (filter=="sample") {
        if (select>n["sample"]) stop("Removing more than available in sample")
        if (!last_select)
          select <- 1:select else
          select <- (n["sample"]-select+1):n["sample"]
        select <- select + n["burn"] + n["adapt"]
      }
    }
  }
  if (any(select>pmwg$samples$idx))
    stop("select specifies iterations not present")
  ok <- 1:pmwg$samples$idx %in% select
  if (remove) ok <- !ok
  pmwg$samples$alpha <- pmwg$samples$alpha[,,ok,drop=FALSE]
  pmwg$samples$theta_mu <- pmwg$samples$theta_mu[,ok,drop=FALSE]
  pmwg$samples$theta_var <- pmwg$samples$theta_var[,,ok,drop=FALSE]
  pmwg$samples$subj_ll <- pmwg$samples$subj_ll[,ok,drop=FALSE]
  pmwg$samples$a_half <- pmwg$samples$a_half[,ok,drop=FALSE]
  pmwg$samples$origin <- pmwg$samples$origin[,ok,drop=FALSE] 
  pmwg$samples$stage <- pmwg$samples$stage[ok]
  pmwg$samples$idx <- sum(ok)
  pmwg
}


get_sigma <- function(samps,filter="samples",thin=thin,subfilter=NULL) 
  # Get sigma matrix samples  
{
  shorten <- function(x,r,d) {
    if (!is.numeric(r)) stop("subfilter must be numeric\n")
    if (length(r)==1) r <- (r+1):dim(x)[d]
    if (min(r)<1 || max(r)>dim(x)[d]) stop("subfilter out of range\n")
    if (d==2) x[,r,drop=FALSE] else x[,,r,drop=FALSE]
  }

  sigma <- samps$samples$theta_var[,,samps$samples$stage==filter]
  if (!is.null(subfilter)) sigma <- shorten(sigma,subfilter,3)
  if (thin >= samps$samples$idx) stop("Thin value to large\n")
  nmc_thin <- seq(thin,dim(sigma)[3],by=thin)
  sigma[,,nmc_thin]  
}



### pmwg object list functions ----

chain_thin <- function(samplers,thin=5) {
  data_list <- attr(samplers,"data_list")  
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  samplers <- lapply(samplers,thin_pmwg,thin=thin)
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  samplers
}

chain_shorten <- function(samplers,n) 
  # removes first n samples (scalar n) or range at some point (vector n) 
{
  data <- attr(samplers,"data_list")
  design <- attr(samplers,"design_list")
  model <- attr(samplers,"model_list")
  samplers <- lapply(samplers,remove_iterations,select=n)
  attr(samplers,"data_list") <- data
  attr(samplers,"design_list") <- design
  attr(samplers,"model_list") <- model
  samplers
}
  
p_names <- function(samplers) 
  # parameter names of a pmwg object or list of pmwg objects  
{
  if (class(samplers)=="pmwg")
    dimnames(samplers$samples$alpha)[[1]] else
    dimnames(samplers[[1]]$samples$alpha)[[1]]
} 

get_sigma_chains <- function(samples,filter="samples",thin=thin,subfilter=NULL)
  # Get sigma matrix samples for chains and put in one array  
{
  sigma <- lapply(samples,get_sigma,filter=filter,thin=thin,subfilter=subfilter)
  abind(sigma)
}

#### mcmc extraction ----

as_Mcmc <- function(sampler,selection=c("alpha","mu","variance","covariance","correlation","LL")[1],
                    filter=stages,thin=1,subfilter=NULL)
  # replacement for pmwg package as_mccmc
  # allows for selection of subj_ll, and specifying selection as an integer
  # allows filter as single integer so can remove first filter samples
  # adds ability to thin and to subfilter (i.e., remove trials after removing 
  # a stage), range of values or single integer (expanded to end)
{
  
  shorten <- function(x,r,d) {
    if (!is.numeric(r)) stop("subfilter must be numeric\n")
    if (length(r)==1) r <- (r+1):dim(x)[d]
    if (min(r)<1 || max(r)>dim(x)[d]) stop("subfilter out of range\n")
    if (d==2) x[,r,drop=FALSE] else x[,,r,drop=FALSE]
  }
  
  if (class(sampler) != "pmwgs") stop("sampler must be an pmwgs object\n")
  stages <- c("burn", "adapt", "sample")
  if (is.numeric(filter) && length(filter==1)) filter <- filter:sampler$samples$idx
  if (is.numeric(selection)) 
    selection <- selection=c("alpha","mu","variance","covariance","correlation","LL")[selection]
  if (selection %in% c("mu","variance","covariance","correlation") &&
      all(is.na((sampler$samples$theta_mu))))
    stop("Cannot select mu, variance, covariance or correlation for non-hierarchical model\n")
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices\n")
  }
  
  if (selection == "mu") {
    mu <- sampler$samples$theta_mu[, filter,drop=FALSE]
    if (is.null(subfilter)) mu <- t(mu) else
      mu <- t(shorten(mu,subfilter,2))
    if (thin > dim(mu)[1]) stop("Thin to large\n")
    out <- coda::mcmc(mu[seq(thin,dim(mu)[1],by=thin),,drop=FALSE])
    attr(out,"selection") <- selection
    return(out)
  } else if (selection %in% c("variance","covariance","correlation")) {
    tvar <- sampler$samples$theta_var[, , filter,drop=FALSE]
    if (selection=="correlation")
      tvar <- array(apply(tvar,3,cov2cor),dim=dim(tvar),dimnames=dimnames(tvar))
    if (selection %in% c("covariance","correlation")) {
      lt <- lower.tri(tvar[,,1])
      pnams <- dimnames(tvar)[[1]]
      pnams <- outer(pnams,pnams,paste,sep=".")[lt]
      tvar <- apply(tvar,3,function(x){x[lt]})
      dimnames(tvar)[[1]] <- pnams
      if (all(tvar==0))
        stop("All covariances/correlations zero, model assumes independence.")
    } else tvar <- apply(tvar,3,diag)  
    if (is.null(subfilter)) tvar <- t(tvar) else
      tvar <- t(shorten(tvar,subfilter,2))
    if (thin > dim(tvar)[1]) stop("Thin to large\n")
    out <- coda::mcmc(tvar[seq(thin,dim(tvar)[1],by=thin),,drop=FALSE])
    attr(out,"selection") <- selection
    return(out)
  } else if (selection == "alpha") {
    alpha <- sampler$samples$alpha[, , filter,drop=FALSE]
    if (!is.null(subfilter)) alpha <- shorten(alpha,subfilter,3)
    if (thin > dim(alpha)[3]) stop("Thin to large\n")
    alpha <- alpha[,,seq(thin,dim(alpha)[3],by=thin),drop=FALSE]
    out <- stats::setNames(lapply(
      seq(dim(alpha)[2]),
      function(x) {
        coda::mcmc(t(alpha[, x, ]))
      }
    ), names(sampler$data))
    attr(out,"selection") <- selection
    return(out)
  } else if (selection == "LL") {
    LL <- sampler$samples$subj_ll[, filter,drop=FALSE]
    if (!is.null(subfilter)) LL <- shorten(LL,subfilter,2)
    if (thin > dim(LL)[2]) stop("Thin to large\n")
    out <- setNames(lapply(
      data.frame(t(LL[,seq(thin,dim(LL)[2],by=thin),drop=FALSE])),coda::mcmc),
    names(sampler$data)) 
    attr(out,"selection") <- selection
    return(out)
  }
  stop("Argument `selection` should be one of mu, variance, covariance, correlation, alpha, LL\n")
}    


as_mcmc.list <- function(samplers,
  selection=c("alpha","mu","variance","covariance","correlation","LL")[1],
  filter="burn",thin=1,subfilter=NULL,natural=FALSE) 
  # Combines as_Mcmc of samplers and returns mcmc.list
  # natural = TRUE transform mu or alpha to natural scale
{
  
  stages <- c("burn", "adapt", "sample")
  if (length(filter)>1) filter <- filter[1]
  if (!class(samplers)=="list" || !all(unlist(lapply(samplers,class))=="pmwgs"))
    stop("samplers must be a list of pmwgs objects")
  if (!all(apply(do.call(rbind,lapply(samplers,function(x){
    x$par_names})),2,function(x){all(x[1]==x[-1])}))) 
    stop("samplers must have the same parameter names")
  if (!all(apply(do.call(rbind,lapply(samplers,function(x){
    x$subjects})),2,function(x){all(x[1]==x[-1])}))) 
    stop("samplers must have the same subjects")
  
  subjects <- names(samplers[[1]]$data)

  mcmcList <- lapply(samplers,as_Mcmc,selection=selection,filter=filter,
                     thin=thin,subfilter=subfilter)
  if (natural) {
    if (selection == "alpha")
      mcmcList <- lapply(mcmcList,function(mcs){lapply(mcs,function(mc){
          attributes(samplers)$design_list[[1]]$model$Ntransform(mc)  
    })}) else if (selection=="mu") mcmcList <- 
      lapply(mcmcList,attributes(samplers)$design_list[[1]]$model$Ntransform) else
    warning("Can only transform alpha or mu to natural")
  }
  if (selection=="LL")
    iter <- unlist(lapply(mcmcList,function(x){length(x[[1]])})) else 
    if (selection == "alpha") 
      iter <- unlist(lapply(mcmcList,function(x){dim(x[[1]])[1]})) else
        iter <- unlist(lapply(mcmcList,function(x){dim(x)[1]}))
  if (!all(iter[1]==iter[-1])) 
    message("Chains have different numbers of samples, using first ",min(iter))
  iter <- min(iter)
  if (selection %in% c("alpha","LL")) {
    nChains <- length(mcmcList)
    ns <- length(samplers[[1]]$subjects)
    out <- vector(mode="list",length=ns)
    lst <- vector(mode="list",length=nChains)
    for (i in 1:ns) {
      for (j in 1:nChains) if (selection=="LL")
        lst[[j]] <- as.mcmc(mcmcList[[j]][[i]][1:iter]) else
        lst[[j]] <- as.mcmc(mcmcList[[j]][[i]][1:iter,])  
      out[[i]] <- mcmc.list(lst)
    }
    out <- setNames(out,subjects)
    attr(out,"selection") <- selection
    return(out)
  } else {
    out <- mcmc.list(mcmcList)
    attr(out,"selection") <- selection
    return(out)
  }
}
