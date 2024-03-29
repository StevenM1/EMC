  
augment = function(s,da,design)
    # Adds attributes to augmented data
    # learn: empty array for Q values with dim = choice alternative (low,high) x 
    #   stimulus x trials (max across stimuli)
    # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
    #   stimulus matrix (rows for each choice alternative contiguous)  
{
  getIndex <- function(typei,cname,da,maxn) {
    out <- which(da[,cname]==typei)
    c(out,rep(NA,maxn-length(out)))
  }
  
  if (!is.null(design$adapt$stimulus)) {
    targets <- design$adapt$stimulus$targets
    par <- design$adapt$stimulus$output_name
    maxn <- max(sapply(dimnames(targets)[[1]],function(x){table(da[da$subjects==s,x])}))
    # da index x stimulus
    out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
                  da=da[da$subjects==s,],maxn=maxn)
    stimulus <- list(index=out)
    # accumulator x stimulus x trials
    stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
                             dimnames=list(rownames(targets),targets[1,],NULL))
    stimulus$targets <- targets
    stimulus$par <- par
  } # add other types here 
  list(stimulus=stimulus)
}

# ARCHITECTURE OF FOLLOWING FUNCTION WILL NEED UPDATING FOR MULTIPLE ADAPT TYPES  
# s="10";npars=pars;da=data; rfun=model$rfun;return_learning=FALSE;mapped_p=FALSE
update_pars = function(s,npars,da,rfun=NULL,return_learning=FALSE,mapped_p=FALSE,
                       return_all=FALSE) 
  # for subject s
  # Either return da filling in responses and RT if da has no responses (in
  # in which case rfun must be supplied), or return npars filling in adapted 
  # parameters or if return_learning returns learning (e.g., Q values)
{
  
  adapt <- attr(da,"adapt")$design
  index <- attr(da,"adapt")[[s]]$stimulus$index 
  learn <- attr(da,"adapt")[[s]]$stimulus$learn 
  npars <- npars[da$subjects==s,]
  da <- da[da$subjects==s,]  
  add_response <- any(is.na(da$R))
  nAcc <- dim(adapt$stimulus$targets)[1]
  namAcc <- dimnames(adapt$stimulus$targets)[[1]]
  nStim <- dim(adapt$stimulus$targets)[2]
  namStim <- colnames(adapt$stimulus$targets)
  # Maximum number of Q updates
  maxQupdates <- dim(index)[1]/nAcc
  # fill an array of parameters, dim = accumulator x trial x 
  # stimulus x parameter type ("v0"    "B"     "t0"    "alpha" "w"     "q0"    "A" )
  parArr <- aperm(
    array(apply(index,2,function(x){npars[x,]}),
          dim=c(nAcc,maxQupdates,dim(npars)[2],nStim),
          dimnames=list(namAcc,NULL,dimnames(npars)[[2]],namStim)),
  c(1,2,4,3)) # reorder to make look up quick in loop
  # fill prob reward array: trials x stimulus x stimulus component
  # pReward <- array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
  #   da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc),
  #   dimnames=list(namAcc,NULL,namStim,namAcc))[1,,,]
  pReward <- array(array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
    da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc))[1,,,],
    dim=c(maxQupdates,nStim,nAcc),dimnames=list(NULL,namStim,namAcc))
  
  # Extract Q values and update
  if (add_response) {
    da_reward <- da_rt <- rep(NA,dim(da)[1])         # Rewards and rts
    da_R <- factor(da_reward,levels=namAcc)          # Response factor
    Ri <- array(index,dim=c(nAcc,maxQupdates,nStim)) # Row indices
  } else {
    # maxQupdates rows, nStim columns
    Rmat <- array(da[index,"R"], dim=c(nAcc,maxQupdates,nStim))[1,,] 
    reward_mat <- array(da[index,"reward"], dim=c(nAcc,maxQupdates,nStim))[1,,] 
  }
  for (i in 1:maxQupdates)  { 
    ok <- !is.na(parArr[1,i,,1]) # stimuli that need updating
    if (i==1) learn[,ok,i] <- parArr[,1,,adapt$stimulus$init_par] # Initialize   
    # pick out fixed pars
    pars <- setNames(data.frame(matrix(parArr[,i,,][,ok,adapt$fixed_pars,drop=FALSE],
      ncol=length(adapt$fixed_pars))),adapt$fixed_pars)
    # calculate and add output_par
    pars[[adapt$stimulus$output_name]] <- adapt$stimulus$output_fun(
      output_pars=matrix(parArr[,i,ok,adapt$stimulus$output_par_names],
                         ncol=length(adapt$stimulus$output_par_names)),
      Q=as.vector(learn[,,i][,ok]))
    if (add_response) { # Simulate trial
      Rrt <- rfun(factor(levels=dimnames(adapt$stimulus$targets)[[1]]),pars)
      Rfac <- Rrt[,"R"]
      reward <- rbinom(length(Rfac),1,pReward[i,ok,][cbind(1:length(Rfac),as.numeric(Rfac))])
      # harvest new trial info
      da_rt[Ri[,i,][,ok,drop=FALSE]] <- rep(Rrt[,"rt"],each=nAcc)
      da_reward[Ri[,i,][,ok,drop=FALSE]] <- rep(reward,each=nAcc)
      da_R[Ri[,i,][,ok,drop=FALSE]] <- rep(Rfac,each=nAcc)
    } else { # Extract trial information
      Rfac <- factor(Rmat[i,ok],levels=namAcc)  
      reward <- reward_mat[i,ok]
    }
    if (i<maxQupdates) { # Update Q value
      learn[,ok,i+1] <- learn[,ok,i] # Copy last for all components
      imat <- cbind(as.numeric(Rfac),1:sum(ok)) # Components to update
      learn[,ok,i+1][imat] <- adapt$stimulus$adapt_fun(Qlast=learn[,,i][,ok,drop=FALSE][imat],
        adapt_par=parArr[,i,,adapt$stimulus$adapt_par][,ok,drop=FALSE][imat],reward)
    } 
    if (mapped_p | !add_response | return_all) # Output will be new so update parrArr 
      parArr[,i,,adapt$stimulus$output_name][!is.na(parArr[,i,,adapt$stimulus$output_name])] <- 
        pars[,adapt$stimulus$output_name]
  }
  if (add_response) {
    da$R <- da_R
    da$reward <- da_reward
    da$rt <- da_rt
    attr(da,"adapt")[[s]]$stimulus$learn <- learn
  }
  if (mapped_p | !add_response | return_all) {
    stim_output <- parArr[,,,adapt$stimulus$output_name][!is.na(index)] 
    # Protect against numerical problems, may screw up dfun for some models
    stim_output[is.na(stim_output)|is.nan(stim_output)] <- -Inf 
    npars[index[!is.na(index)],adapt$stimulus$output_name] <- stim_output
  }
  if (return_all) return(list(learn=learn,pars=npars,data=da)) 
  if (return_learning) return(learn)
  if (mapped_p) return(list(data=da,pars=npars))
  if (add_response) return(da)
  npars
}

# data=dadm
adapt_data <- function(data,design,model,pars,
  add_response=FALSE,return_learning=FALSE,mapped_p=FALSE,return_all=FALSE)
  # runs learning, and either returns learning, or data with responses and rt 
  # added (either if data has NA in R or if add_Response), or adapted parameters 
  # or data + parameters (mapped_p)   
{
  if (return_all & mapped_p) {
    warning("return_all and mapped_p incompatible, going with the latter")
    return_all <- FALSE
  }
  if (add_response) data$R <- NA # Force new data generation ignoring old data
  if ( all(is.na(data$R)) ) add_response <- TRUE # make new if none supplied
  # add augmentation
  if (is.null(attr(data,"adapt"))) {
    attr(data,"adapt") <- setNames(
      lapply(levels(data$subjects),augment,da=data,design=design),
    levels(data$subjects))
    attr(data,"adapt")$design <- design$adapt  
  }
  daList <- setNames(
    lapply(levels(data$subjects),update_pars,npars=pars,da=data,return_all=return_all,
           rfun=model$rfun,return_learning=return_learning,mapped_p=mapped_p),
    levels(data$subjects))
  if (return_learning) return(daList)
  adapt <- attr(data,"adapt")
  if (mapped_p | return_all) {
    data <- do.call(rbind,lapply(daList,function(x)x$data))
    pars <- do.call(rbind,lapply(daList,function(x)x$pars))
    if (mapped_p) return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
    for (i in names(daList)) adapt[[i]] <- attr(daList[[i]]$data,"adapt")[[i]]
    data <- cbind(data,pars)
    attr(data,"adapt") <- adapt
    return(list(learn=lapply(daList,function(x)x$learn),data=data))
  }
  for (i in names(daList)) adapt[[i]] <- attr(daList[[i]],"adapt")[[i]]
  data <- do.call(rbind,daList)
  attr(data,"adapt") <- adapt
  data
}
