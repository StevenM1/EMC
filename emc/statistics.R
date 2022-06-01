#### chain statistics ----

chain_n <- function(samplers) 
  # Length of stages for each chain
{
  do.call(rbind,lapply(samplers, function(x){
    table(factor(x$samples$stage,levels=c("burn","adapt","sample")))
  }))
}


check_adapt <- function(samplers,verbose=TRUE) 
  # Checks chains to see if adapted  
{
  ok <- TRUE
  for (i in 1:length(samplers)) {
    res <- attr(samplers[[i]],"adapted")
    if (is.null(res) || is.character(res)) {
      ok <- FALSE
      if (verbose) message("Chain ",i," not adapted") 
    } else if (verbose)
      message("Chain ",i," adapted by iteration ", res)
  }
  ok
}

es_pmwg <- function(pmwg_mcmc,selection="alpha",summary_alpha=mean,
                    filter="burn",thin=1,subfilter=NULL)
  # Effective size 
{
  if (!(class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list"))) {
    if (class(pmwg_mcmc)=="pmwgs") 
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter)                        
  }
  if (attr(pmwg_mcmc,"selection")=="LL") 
    stop("Effective size not sensible for LL\n")
  out <- do.call(rbind,lapply(pmwg_mcmc,effectiveSize))
  if (attr(pmwg_mcmc,"selection")=="alpha") {
    if (!is.null(summary_alpha)) out <- apply(out,2,summary_alpha)
    out 
  } else apply(out,2,sum) 
}
  

# return_summary=FALSE;print_summary=TRUE;digits_print=2;sort_print=TRUE;
# autoburnin=FALSE;transform=TRUE
# selection="alpha";filter="burn";thin=1;subfilter=NULL
# pmwg_mcmc=sVat0;selection="mu";filter="sample"
gd_pmwg <- function(pmwg_mcmc,return_summary=FALSE,print_summary=TRUE,
    digits_print=2,sort_print=TRUE,autoburnin=FALSE,transform=TRUE,
    selection="alpha",filter="burn",thin=1,subfilter=NULL,natural=FALSE) 
  # R hat, prints multivariate summary returns each participant unless +
  # multivariate as matrix unless !return_summary
{
  
  gelman_diag_robust <- function(mcl,autoburnin,transform) 
  {
    gd <- try(gelman.diag(mcl,autoburnin=autoburnin,transform=transform),silent=TRUE)
    if (class(gd)=="try-error") list(mpsrf=Inf,psrf=matrix(Inf)) else gd
  }
  
  if ( selection=="LL" ) stop("Rhat not appropriate for LL") 
  if (class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list")) {
    if (natural) warning("Cannot transform to natural scale unless samples list provided")
  } else {
    if (class(pmwg_mcmc)=="pmwgs") {
      if (natural) warning("Cannot transform to natural scale unless samples list provided")
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter)   
    } else
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter,natural=natural)
  } 
  if (selection=="alpha") {
    gd <- lapply(pmwg_mcmc,gelman_diag_robust,autoburnin = autoburnin, transform = transform) 
    out <- unlist(lapply(gd,function(x){x$mpsrf}))
  } else {
    gd <- gelman_diag_robust(pmwg_mcmc,autoburnin = autoburnin, transform = transform)
    out <- gd$mpsrf 
  }
  if (return_summary) return(out)
  if (sort_print) out <- sort(out)
  if (print_summary) print(round(out,digits_print))
  if (selection=="alpha") invisible(
    cbind(do.call(rbind,lapply(gd,function(x){x[[1]][,1]})),
          mpsrf=unlist(lapply(gd,function(x){x[[2]]})))) else
    invisible(c(gd$psrf[,1],mpsrf=gd$mpsrf))
}


iat_pmwg <- function(pmwg_mcmc,
    print_summary=TRUE,digits_print=2,sort_print=TRUE,summary_alpha=mean,
    selection="alpha",filter="burn",thin=1,subfilter=NULL) 
  # Integrated autocorrelation time, prints multivariate summary returns each participant unless +
  # multivariate as matrix unless !return_summary
{
  
  IAT <- function (x,verbose=FALSE) 
  # From LaplacesDemon
  {
    dt <- x
    n <- length(x)
    mu <- mean(dt)
    s2 <- var(dt)
    maxlag <- max(3, floor(n/2))
    Ga <- rep(0, 2)
    Ga[1] <- s2
    lg <- 1
    Ga[1] <- Ga[1] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    m <- 1
    lg <- 2 * m
    Ga[2] <- sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    lg <- 2 * m + 1
    Ga[2] <- Ga[2] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    IAT <- Ga[1]/s2
    while ((Ga[2] > 0) & (Ga[2] < Ga[1])) {
      m <- m + 1
      if (2 * m + 1 > maxlag) {
        if (verbose) cat("Not enough data, maxlag=", maxlag, "\n")
        break
      }
      Ga[1] <- Ga[2]
      lg <- 2 * m
      Ga[2] <- sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
      lg <- 2 * m + 1
      Ga[2] <- Ga[2] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
        IAT <- IAT + Ga[1]/s2
    }
    IAT <- -1 + 2 * IAT
    return(IAT)
  }
  
  get_IAT <- function(mcs) {
    if (class(mcs) != "mcmc.list") apply(mcs,2,IAT) else
      apply(do.call(rbind,lapply(mcs,function(x){apply(x,2,IAT)})),2,mean)
  }
  
  if (!(class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list"))) {
    if (class(pmwg_mcmc)=="pmwgs") 
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter)                        
  }
  if ( selection=="LL" ) stop("IAT not appropriate for LL") else
  if (selection=="alpha") {
    out <- do.call(rbind,lapply(pmwg_mcmc,get_IAT) ) 
    if (!is.null(summary_alpha)) out <- apply(out,2,summary_alpha)
  } else out <- get_IAT(pmwg_mcmc)
  if (sort_print & selection != "alpha") out <- sort(out)
  if (print_summary) print(round(out,digits_print))
  invisible(out)
}

#### Posterior parameter tests ----

# y=NULL;natural=TRUE;x_name=NULL;y_name=NULL; c_vector=NULL
# mu=0;alternative = c("less", "greater")[1];
# probs = c(0.025,.5,.975);digits=2;p_digits=3;print_table=TRUE;
# x_filter="burn";x_selection="alpha";x_subfilter=1;
# y_filter="burn";y_selection="alpha";y_subfilter=1
# 
# x=ddmPNASa;p_name="a_Ea-n";x_selection = "mu";x_filter="sample"
# c_vector=c(0,1,-1,0,1,-1)/4
# x=andrew;p_name="t0_CIc-i";x_selection="alpha";x_filter="burn";x_name="D"
  
p_test <- function(x,y=NULL,p_name,natural=TRUE,c_vector=NULL,
                   x_name=NULL,y_name=NULL,
                   mu=0,alternative = c("less", "greater")[1],
                   probs = c(0.025,.5,.975),digits=2,p_digits=3,print_table=TRUE,
                   x_filter="burn",x_selection="alpha",x_subfilter=1,
                   y_filter="burn",y_selection="alpha",y_subfilter=1) {

  

  get_effect <- function(x,x_vector,x_name,p_name,Ntransform=NULL,c_vector=NULL) 
    # Effect, on natural scale if Ntransform supplied, and always on natural 
    # if c_vector supplied.
  {
    if (is.null(x_name)) 
      x <- do.call(rbind,x) else
      x <- do.call(rbind,x[[x_name]])
    p_root <- strsplit(p_name,"_")[[1]]
    if (!is.null(c_vector)) { 
      design <- attr(c_vector,"design")
      model <- design$model
      design$Ffactors$subjects <- design$Ffactors$subjects[1]
      dadm <- design_model(make_data(x[1,],design,model,trials=1),design,model,
                           rt_check=FALSE,compress=FALSE)
      cvals <- apply(x,1,function(x){get_pars(x,dadm)[,p_root[1]]})
      # if (!is.null(Ntransform)) {
      #   cvals <- t(cvals)
      #   dimnames(cvals)[[2]] <- rep(p_root[1],dim(cvals)[2])
      #   cvals <- t(Ntransform(cvals))
      # }
      return((c_vector %*% cvals)[1,])  
    }
    if (length(p_root)==1) { # intercept 
     if (!is.null(Ntransform)) {
      return(Ntransform(x[,p_root,drop=FALSE])[,1])    
     } else return(x[,p_root]) 
    } else p_root <- p_root[1]
    cv <- attr(x_vector,"map")[[p_root]][,p_name] # Contrast vector
    cpos <- cv>0
    if (!any(cpos)) cpos <- Ntransform(x[,p_root,drop=FALSE]) else
      cpos <- Ntransform(x[,p_root,drop=FALSE] + sum(cv[cpos])*x[,p_name,drop=FALSE])
    cneg <- cv<0
    if (!any(cneg)) cneg  <- Ntransform(x[,p_root,drop=FALSE]) else
      cneg <- Ntransform(x[,p_root,drop=FALSE] + sum(cv[cneg])*x[,p_name,drop=FALSE])
    cpos - cneg
  }


  if (class(x[[1]])!="pmwgs") stop("x must be a list of pmwgs objects") 
  design <- attr(x,"design_list")[[1]]
  if (!is.null(c_vector)) attr(c_vector,"design") <- design
  x_vector <- sampled_p_vector(design)
  if (!natural) Ntransform <- NULL else 
    Ntransform <- attr(x,"design_list")[[1]]$model$Ntransform
  x <- as_mcmc.list(x,selection=x_selection,filter=x_filter,subfilter=x_subfilter)
  if (x_selection != "alpha") x_name <- NULL else
    if (is.null(x_name)) x_name <- 1 else
    if (!(x_name %in% names(x))) stop("Subject x_name not in x")
  x <- get_effect(x,x_vector,x_name,p_name,Ntransform,c_vector)
  if (is.null(y)) {
    p <- mean(x<mu)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),c(NA,mu,NA))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(p_name,"mu")
  } else {
    if (class(y[[1]])!="pmwgs") stop("y must be a list of pmwgs objects") 
    y_vector <- sampled_p_vector(attr(y,"design_list")[[1]])
    y <- as_mcmc.list(y,selection=y_selection,filter=y_filter,subfilter=y_subfilter)
    if (y_selection != "alpha") y_name <- NULL else
      if (is.null(y_name)) y_name <- 1 else
        if (!(y_name %in% names(y))) stop("Subject y_name not in y")
    y <- get_effect(y,y_vector,y_name,p_name,Ntransform)
    if (length(x)>length(y)) x <- x[1:length(y)] else y <- y[1:length(x)]
    d <- x-y
    p <- mean(d<0)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),quantile(y,probs),quantile(d,probs))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(paste(p_name,c(x_name,y_name),sep="_"),
                            paste(x_name,y_name,sep="-"))
  }
  if (print_table) {
    ptab <- tab
    ptab <- round(ptab,digits)
    attr(ptab,alternative) <- round(attr(ptab,alternative),p_digits)
    print(ptab)
  }
  invisible(tab)
}


# y=NULL;p_name;mapped=FALSE;c_vector=NULL;
#                    x_name=NULL;y_name=NULL;
#                    mu=0;alternative = c("less", "greater")[1];
#                    probs = c(0.025,.5,.975);digits=2;p_digits=3;print_table=TRUE;
#                    x_filter="sample";x_selection="alpha";x_subfilter=0;
#                    y_filter="sample";y_selection="alpha";y_subfilter=0
# 
# x=samples; p_name="mean_FWwords:Sold";x_selection = "mu"

p_test <- function(x,y=NULL,p_name,mapped=FALSE,c_vector=NULL,
                   x_name=NULL,y_name=NULL,
                   mu=0,alternative = c("less", "greater")[1],
                   probs = c(0.025,.5,.975),digits=2,p_digits=3,print_table=TRUE,
                   x_filter="sample",x_selection="alpha",x_subfilter=0,
                   y_filter="sample",y_selection="alpha",y_subfilter=0) {

  

  get_effect <- function(x,x_vector,x_name,p_name,c_vector=NULL) 
    # Effect, must always be on mapped scale if c_vector supplied.
  {
    if (is.null(x_name)) 
      x <- do.call(rbind,x) else
      x <- do.call(rbind,x[[x_name]])
    p_root <- strsplit(p_name,"_")[[1]]
    if (!is.null(c_vector)) { 
      design <- attr(c_vector,"design")
      model <- design$model
      design$Ffactors$subjects <- design$Ffactors$subjects[1]
      dadm <- design_model(make_data(x[1,],design,model,trials=1),design,model,
                           rt_check=FALSE,compress=FALSE)
      cvals <- apply(x,1,function(x){get_pars(x,dadm)[,p_root[1]]})
      return((c_vector %*% cvals)[1,])  
    }
    if (length(p_root)==1) { # intercept 
      return(x[,p_root]) 
    } else p_root <- p_root[1]
    cv <- attr(x_vector,"map")[[p_root]][,p_name] # Contrast vector
    cpos <- cv>0
    if (!any(cpos)) cpos <- x[,p_root,drop=FALSE] else
      cpos <- x[,p_root,drop=FALSE] + sum(cv[cpos])*x[,p_name,drop=FALSE]
    cneg <- cv<0
    if (!any(cneg)) cneg  <- x[,p_root,drop=FALSE] else
      cneg <- x[,p_root,drop=FALSE] + sum(cv[cneg])*x[,p_name,drop=FALSE]
    cpos - cneg
  }


  if (mapped & !(selection %in% c("mu","alpha")))
    stop("Can only analyze mapped mu or alpha parameters")
  if (class(x[[1]])!="pmwgs") stop("x must be a list of pmwgs objects") 
  design <- attr(x,"design_list")[[1]]
  if (!is.null(c_vector)) attr(c_vector,"design") <- design
  x_vector <- sampled_p_vector(design)
  x <- as_mcmc.list(x,selection=x_selection,filter=x_filter,
                    subfilter=x_subfilter,mapped=mapped,add_constants=TRUE)
  if (x_selection != "alpha") x_name <- NULL else
    if (is.null(x_name)) x_name <- names(x)[1] 
  if (!is.null(x_name)) {
    if (!(x_name %in% names(x))) stop("Subject x_name not in x")
    message("Testing x subject ",x_name)
  }
  x <- get_effect(x,x_vector,x_name,p_name,c_vector)
  if (is.null(y)) {
    p <- mean(x<mu)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),c(NA,mu,NA))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(p_name,"mu")
  } else {
    if (class(y[[1]])!="pmwgs") stop("y must be a list of pmwgs objects") 
    y_vector <- sampled_p_vector(attr(y,"design_list")[[1]])
    y <- as_mcmc.list(y,selection=y_selection,filter=y_filter,
                      subfilter=y_subfilter,mapped=mapped)
    if (y_selection != "alpha") y_name <- NULL else
      if (is.null(y_name)) y_name <- names(y)[1]
    if (!is.null(y_name)) {
      if (!(y_name %in% names(y))) stop("Subject y_name not in y")
      message("Testing y subject ",y_name)  
    }
    y <- get_effect(y,y_vector,y_name,p_name)
    if (length(x)>length(y)) x <- x[1:length(y)] else y <- y[1:length(x)]
    d <- x-y
    p <- mean(d<0)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),quantile(y,probs),quantile(d,probs))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(paste(p_name,c(x_name,y_name),sep="_"),
                            paste(x_name,y_name,sep="-"))
  }
  if (print_table) {
    ptab <- tab
    ptab <- round(ptab,digits)
    attr(ptab,alternative) <- round(attr(ptab,alternative),p_digits)
    print(ptab)
  }
  invisible(tab)
}


# filter="burn";subfilter=0;use_best_fit=FALSE;print_summary=FALSE;digits=0
# filter="sample"
pmwg_IC <- function(samplers,filter="burn",subfilter=0,use_best_fit=FALSE,
                    print_summary=FALSE,digits=0,subject=NULL)
  # Gets DIC, BPIC, effective parameters, mean deviance, and deviance of mean  
{
  
  
  # Mean log-likelihood for each subject
  if (class(samplers)=="pmwgs") 
      ll <- as_Mcmc(samplers,selection="LL",filter=filter,subfilter=subfilter) else
      ll <- as_mcmc.list(samplers,selection="LL",filter=filter,subfilter=subfilter)
  minDs <- -2*unlist(lapply(ll,function(x){max(unlist(x))}))
  mean_lls <- unlist(lapply(ll,function(x){mean(unlist(x))}))
  
  if (class(samplers)=="pmwgs") 
      alpha <- as_Mcmc(samplers,selection="alpha",filter=filter,subfilter=subfilter) else
      alpha <- as_mcmc.list(samplers,selection="alpha",filter=filter,subfilter=subfilter)                        
  mean_pars <- lapply(alpha,function(x){apply(do.call(rbind,x),2,mean)})
  # log-likelihood for each subject using their mean parameter vector
  ll_func <- attr(samplers,"design_list")[[1]]$model$log_likelihood
  data <- samplers[[1]]$data
  mean_pars_lls <- setNames(numeric(length(mean_pars)),names(mean_pars))
  for (sub in names(mean_pars))
    mean_pars_lls[sub] <- ll_func(mean_pars[[sub]],dadm = data[[sub]])
  Dmeans <- -2*mean_pars_lls
  if (use_best_fit) minDs <- pmin(minDs,Dmeans)
  
  if (!is.null(subject)) {
    Dmeans <- Dmeans[subject]
    mean_lls <- mean_lls[subject]
    minDs <- minDs[subject]
  }
  
  # mean deviance(-2*ll of all data) 
  mD <- sum(-2 * mean_lls)
  # Deviance of mean
  Dmean <- sum(Dmeans)
  # mimimum Deviance
  minD <- sum(minDs)
  
  # Use deviance of mean as best fit or use actual best fit
  if (!use_best_fit) Dm <- Dmean else Dm <- minD
   
  # effective number of parameters
  pD <- mD - Dm
  # DIC = mean deviance + effective number of parameters
  DIC <- mD + pD
  # BPIC = mean deviance + 2*effective number of parameters 
  # Note this is the "easy" BPIC, instead of the complex 2007 one
  BPIC <- mD + 2*pD
  out <- c(DIC = DIC, BPIC = BPIC, EffectiveN = pD,meanD=mD,Dmean=Dmean,minD=minD)
  names(out) <- c("DIC","BPIC","EffectiveN","meanD","Dmean","minD")
  if (print_summary) print(round(out,digits))
  invisible(out)
}



compare_ICs <- function(sList,filter="burn",subfilter=0,use_best_fit=FALSE,
                    print_summary=FALSE,digits=0,digits_p=3,subject=NULL) {
  
  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }
  
  ICs <- data.frame(do.call(rbind,
    lapply(sList,pmwg_IC,filter=filter,subfilter=subfilter,
           use_best_fit=use_best_fit,subject=subject)))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])
  if (print_summary) {
    tmp <- out
    tmp$wDIC <- round(tmp$wDIC,digits_p)
    tmp$wBPIC <- round(tmp$wBPIC,digits_p)
    tmp[,-c(2,4)] <- round(tmp[,-c(2,4)])
    print(tmp)
  }
  invisible(out)
}
