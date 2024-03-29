#### Likelihoods ----


#### Race likelihoods ----

log_likelihood_race <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
{
  
  pars <- get_pars(p_vector,dadm)
  if (is.null(attr(pars,"ok"))) 
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")$dfun(rt=dadm$rt[dadm$winner],
                                                  pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-attr(dadm,"model")$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- min_ll    #SM changed this from 0!
  lds <- lds[attr(dadm,"expand")] # decompress
  if (n_acc>1) {
    winner <- dadm$winner[attr(dadm,"expand")]
    ll <- lds[winner]
    if (n_acc==2)
      ll <- ll + lds[!winner] else
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    ll[is.na(ll)] <- min_ll          # 0  SM changed this from 0! exp(0) = pretty good
#    print(sum(pmax(min_ll, ll)))
    sum(pmax(min_ll,ll))
  } else sum(pmax(min_ll,lds))
}




log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10)) 
  # Race model summed log likelihood   
{
  
  f <- function(t,p) 
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.
  {
    if (!is.matrix(p)) p <- matrix(p,nrow=1,dimnames=list(NULL,names(p)))
    out <- attr(dadm,"model")$dfun(t,
                                   matrix(rep(p[1,],each=length(t)),nrow=length(t),
                                          dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-attr(dadm,"model")$pfun(t,
                                            matrix(rep(p[i,],each=length(t)),nrow=length(t),
                                                   dimnames=list(NULL,dimnames(p)[[2]]))))
    out
  }
  
  logTrunc_race <- function(LT,UT,dadm,pars) 
    # log truncation probability for each rt
  {
    if ( (LT > 0) | (is.finite(UT)) ) 
      # Get truncation adjustment. same length as winner
    {
      upars <- pars[attr(dadm,"unique_nort"),]
      looser <- matrix(!dadm[attr(dadm,"unique_nort"),"winner"],nrow=n_acc)
      np <- dim(upars)[1]/n_acc
      upars <- array(upars,dim=c(n_acc,np,dim(upars)[2]),
                     dimnames=list(NULL,NULL,dimnames(upars)[[2]]))
      lps <- numeric(np)
      for (i in 1:np) {
        tmp <- try(integrate(f,lower=LT,upper=UT,p=upars[,i,][order(looser[,i]),]),silent=TRUE)
        if (class(tmp)=="try-error" || tmp$value == 0 ) lps[i] <- Inf else
          lps[i] <- log(tmp$value)
      }
      lps[is.na(lps)] <- Inf
      return(lps[attr(dadm,"expand_nort")])
    } else return(0)
  }
  
  # Censoring
  logCens_race <- function(LT,UT,dadm,pars,n_trials,n_acc,upper=TRUE) 
    # b = censor bound, bT = truncation bound  
  {
    if (upper) 
      b <- attr(dadm,"UC")[1] else 
        b <- attr(dadm,"LC")[1]
      out <- numeric(n_trials)
      if (is.null(b)) return(out)
      if (upper) ok <- dadm$rt==Inf else ok <- dadm$rt==-Inf
      if (!any(ok)) return(out)
      d <- dadm[ok,]
      p <- pars[ok,]
      d[,"winner"][is.na(d[,"R"])] <- NA
      looser <- matrix(!d[,"winner"],nrow=n_acc)
      np <- dim(p)[1]/n_acc
      p <- array(p,dim=c(n_acc,np,dim(p)[2]),
                 dimnames=list(NULL,NULL,dimnames(p)[[2]]))
      if (upper) {
        lo <- b
        up <- UT
      } else {
        lo <- LT
        up <- b
      }
      lps <- numeric(np)
      for (i in 1:np) {
        if ( all(is.na(looser[,i])) ) {
          l <- !logical(length(looser[,i]))
          num <- den <- vector(mode="list",length=np)
          for (j in 1:length(looser[,i])) {
            temp <- l; temp[j] <- FALSE
            num[[j]] <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(temp),]),silent=TRUE)
            if (class(num[[j]])=="try-error") {num <- NA; break}
            den[[j]] <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(temp),]),silent=TRUE) 
            if (class(den[[j]])=="try-error" || den[[j]]$value==0) {den <- NA; break}
          }
          if (any(is.na(num)) || any(is.na(den))) lps[i] <- -Inf else {
            num <- unlist(lapply(num,function(x)x$value))
            den <- unlist(lapply(den,function(x)x$value))
            nd <- num/den
            w <- num/sum(num)
            lps[i] <- log(sum(w*nd))
          }
        } else {
          num <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(looser[,i]),]),silent=TRUE)
          if (class(num)=="try-error") lps[i] <- -Inf else {
            den <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(looser[,i]),]),silent=TRUE)
            if (class(den)=="try-error" || den$value == 0) lps[i] <- -Inf else
              lps[i] <- log(num$value/den$value)
          }
        }
      }
      lps[is.na(lps)] <- -Inf
      if (upper) return(c(0,lps)[attr(dadm,"expand_uc")]) else
        return(c(0,lps)[attr(dadm,"expand_lc")]) 
  }    
  
  pars <- get_pars(p_vector,dadm)
  
  # Density for winner rows
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[attr(dadm,"ok_dadm_winner")] <- log(attr(dadm,"model")$dfun(
    rt=dadm$rt[attr(dadm,"ok_dadm_winner")],pars=pars[attr(dadm,"ok_dadm_winner"),]))
  
  n_acc <- length(levels(dadm$R))                # number of accumulators
  n_trials <- length(attr(dadm,"expand_winner")) # number of trials
  
  # Truncation
  LT <- attr(dadm,"LT")[1] 
  UT <- attr(dadm,"UT")[1] 
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf
  # log truncation
  lt <- logTrunc_race(LT,UT,dadm,pars)
  if (is.character(lt)) return(min_ll*length(attr(dadm,"expand")))
  if (length(lt)>1) lt[!attr(dadm,"ok_trial")] <- 0
  # log upper censor
  luc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc)
  if (is.character(luc)) return(min_ll*length(attr(dadm,"expand")))
  # log lower censor
  llc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc,upper=FALSE)
  if (is.character(llc)) return(min_ll*length(attr(dadm,"expand")))
  
  # Survivors and truncation
  if (n_acc>1) {
    lds[attr(dadm,"ok_dadm_looser")] <- log(1-attr(dadm,"model")$pfun(
      rt=dadm$rt[attr(dadm,"ok_dadm_looser")],pars=pars[attr(dadm,"ok_dadm_looser"),]))
  }
  lds <- lds[attr(dadm,"expand")]  # decompress 
  lds[is.na(lds)] <- 0
  if (n_acc>1) {
    ll <- numeric(n_trials)
    ll[attr(dadm,"ok_trials")] <- lds[attr(dadm,"ok_da_winner")]
    if (n_acc==2) 
      ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] + lds[attr(dadm,"ok_da_looser")] else
        ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] + 
      apply(matrix(lds[attr(dadm,"ok_da_looser")],nrow=n_acc-1),2,sum)
    ll <- ll - lt + luc + llc
    ll[is.na(ll)] <- 0
    sum(pmax(min_ll,ll))
  } else sum(pmax(min_ll,lds - lt + luc + llc))
}



log_likelihood_ddm <- function(p_vector,dadm,min_ll=log(1e-10)) 
  # DDM summed log likelihood, with protection against numerical issues    
{
  pars <- get_pars(p_vector,dadm)
  like <- numeric(dim(dadm)[1]) 
  if (any(attr(pars,"ok"))) {
    like[attr(pars,"ok")] <- attr(dadm,"model")$dfun(dadm$rt[attr(pars,"ok")],dadm$R[attr(pars,"ok")],pars[attr(pars,"ok"),,drop=FALSE])
  }
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}



log_likelihood_DAS <- function(p_vector,dadm,min_ll=log(1e-10),
  p_names=c("B","v")) # names of race parameters, defualt is wald, (c"mean","sd") for normal
  # Discrete activation suppression: race winner time convolved with exponential
  # with mean tau
{
  pars <- get_pars(p_vector,dadm)
  n_acc <- length(levels(dadm$R))
  ntrials <- dim(dadm)[1]/n_acc
  like <- numeric(ntrials)
  pwin <- pars[dadm$winner,c(p_names,"tau")]
  ploose <- array(pars[!dadm$winner,p_names],dim=c(n_acc-1,ntrials,2))
  rt <- dadm[dadm$winner,"rt"] - pars[dadm$winner,"t0"]
  for (i in 1:ntrials)
    like[i] <- dDAS(rt=rt[i],pwin=pwin[i,1:2],
                    ploose=matrix(ploose[,i,],ncol=2),tau=pwin[i,3],
                    dfun=attr(dadm,"model")$dfun,pfun=attr(dadm,"model")$pfun)
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))  
}



#### sdt choice likelihoods ----

log_likelihood_sdt <- function(p_vector,dadm,lb=-Inf,min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{
  
  pars <- get_pars(p_vector,dadm)
  first <- dadm$lR==levels(dadm$lR)[1]
  last <- dadm$lR==levels(dadm$lR)[length(levels(dadm$lR))] 
  pars[last,"threshold"] <- Inf
  # upper threshold 
  ut <- pars[dadm$winner,"threshold"] 
  # lower threshold fixed at lb for first response
  pars[first &  dadm$winner,"threshold"] <- lb
  # otherwise threshold of response before one made
  notfirst <- !first &  dadm$winner
  pars[notfirst,"threshold"] <- pars[which(notfirst)-1,"threshold"]
  lt <- pars[dadm$winner,"threshold"]
  # fix race expand to suit SDT
  nr <- length(levels(dadm$lR))      # number of responses
  ne <- length(attr(dadm,"expand"))  # length of expand
  # Shorten expand to only one per lR set
  expand <- (attr(dadm,"expand")[(c(1:ne) %% nr)==0] + 1 ) %/% nr
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(attr(dadm,"model")$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,]))
  } else ll <- log(attr(dadm,"model")$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,]))
  ll <- ll[expand]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}
