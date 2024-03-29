# Parameter transformation and mapping

pmat <- function(p_vector,design) 
    # puts vector form of p_vector into matrix form  
{
  ss <- design$Ffactors$subjects
  matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
           dimnames=list(ss,names(p_vector)))
}
    

  
add_constants <- function(p,constants) 
  # augments parameter matrix or vector p with constant parameters (also used in data)
{
  if (is.null(constants)) return(p)
  if (is.matrix(p)) {
    nams <- c(dimnames(p)[[2]],names(constants))
    p <- cbind(p,matrix(rep(constants,each=dim(p)[1]),nrow=dim(p)[1]))
    dimnames(p)[[2]] <- nams
    p
  } else c(p,constants)
}

get_pars <- function(p_vector,dadm) {
  # Add constants, transform p_vector, map to design, transform mapped parameters 
  # to the natural scale, and create trial-dependent parameters
  attr(dadm,"model")$Ttransform(
    attr(dadm,"model")$Ntransform(
      map_p(
        attr(dadm,"model")$transform(add_constants(p_vector,attr(dadm,"constants"))),
      dadm)),
  dadm)
}


add_constants_mcmc <- function(p,constants) 
  mcmc(add_constants(p,constants))


mapped_name_list <- function(design,model,save_design=FALSE)
  # makes a list, with entries for each parameter type, of names for mapped
  # parameters or with unique design columns
{
  doMap <- function(mapi,pmat) t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))
  
  constants <- design$constants
  p_vector <- attr(design,"p_vector")
  mp <- mapped_par(p_vector,design)
  map <- attr(sampled_p_vector(design),"map")
  pmat <- model$transform(add_constants(t(as.matrix(p_vector)),constants))
  plist <- lapply(map,doMap,pmat=pmat)
  if (model$type=="SDT") {
    ht <- apply(map$threshold[,grepl("lR",dimnames(map$threshold)[[2]]),drop=FALSE],1,sum)
    plist$threshold <- plist$threshold[,ht!=max(ht),drop=FALSE]
  }
  # Give mapped variables names and remove constants
  for (i in 1:length(plist)) {
    vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
    if (is.null(vars)) dimnames(plist[[i]])[2] <- names(plist)[i] else {
      uniq <- !duplicated(apply(mp[,vars],1,paste,collapse="_"))
      if (save_design) plist[[i]] <- mp[uniq,vars[-1]] else
        dimnames(plist[[i]])[[2]] <- 
          paste(vars[1],apply(mp[uniq,vars[-1],drop=FALSE],1,paste,collapse="_"),sep="_")
    }
  }
  if (save_design) plist else lapply(plist,function(x){dimnames(x)[[2]]}) 
}


# mcmc=mcmcList[[1]]; design=attr(samplers,"design_list")[[1]];model=attr(samplers,"model_list")[[1]]
map_mcmc <- function(mcmc,design,model, include_constants = TRUE) 
  # Maps vector or matrix (usually mcmc object) of sampled parameters to native 
  # model parameterization. 
{
  doMap <- function(mapi,pmat) t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))
  
  get_p_types <- function(nams)
    unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

  if (!is.matrix(mcmc)) mcmc <- t(as.matrix(mcmc))
  map <- attr(sampled_p_vector(design),"map")
  constants <- design$constants
  mp <- mapped_par(mcmc[1,],design)
  pmat <- model$transform(add_constants(mcmc,constants))
  plist <- lapply(map,doMap,pmat=pmat)
  if (model$type=="SDT") {
    ht <- apply(map$threshold[,grepl("lR",dimnames(map$threshold)[[2]]),drop=FALSE],1,sum)
    plist$threshold <- plist$threshold[,ht!=max(ht),drop=FALSE]
  }
  # Give mapped variables names and flag constant
  isConstant <- NULL
  for (i in 1:length(plist)) {
    vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
    uniq <- !duplicated(apply(mp[,vars],1,paste,collapse="_"))
    if (is.null(vars)) dimnames(plist[[i]])[2] <- names(plist)[i] else {
      dimnames(plist[[i]])[[2]] <- 
        paste(vars[1],apply(mp[uniq,vars[-1],drop=FALSE],1,paste,collapse="_"),sep="_")
    }
    if (dim(plist[[i]])[1]!=1) isConstant <- c(isConstant,
                                               apply(plist[[i]],2,function(x){all(x[1]==x[-1])}))
  }
  pmat <- do.call(cbind,plist)
  cnams <- dimnames(pmat)[[2]]
  dimnames(pmat)[[2]] <- get_p_types(cnams) 
  out <- model$Ntransform(pmat)[,1:length(cnams)]
  dimnames(out)[[2]] <- cnams
  if(!include_constants) out <- out[,!isConstant]
  out <- as.mcmc(out)
  attr(out,"isConstant") <- isConstant
  out
}


mapped_par <- function(p_vector,design,model=NULL,
                       digits=3,remove_subjects=TRUE) 
  # Show augmented data and corresponding mapped parameter  
{
  if (is.null(model)) if (is.null(design$model)) 
    stop("Must specify model as not in design") else model <- design$model
  if (remove_subjects) design$Ffactors$subjects <- design$Ffactors$subjects[1]
  if (!is.matrix(p_vector)) p_vector <- pmat(p_vector,design) 
  dadm <- design_model(make_data(p_vector,design,model,trials=1),design,model,
                       rt_check=FALSE,compress=FALSE)
  ok <- !(names(dadm) %in% c("subjects","trials","R","rt","winner"))
  out <- cbind(dadm[,ok],round(get_pars(p_vector,dadm),digits))
  if (model$type=="SDT")  out <- out[dadm$lR!=levels(dadm$lR)[length(levels(dadm$lR))],]
  out
}


get_design_matrix <- function(samples) 
  attr(sampled_p_vector(attr(samples,"design_list")[[1]]),"map")


