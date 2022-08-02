require(rtdists)

rDDM <- function(lR,pars,precision=3,ok=NULL) 
  # lR is an empty latent response factor lR with one level for each boundary
  # pars is a matrix of parameter values named as in p_types
  # lower is mapped to first level of lR and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # lR <- factor(c("left","right"))
  
{
  if (is.null(ok)) ok <- rep(TRUE,length(lR))
  bad <- rep(NA,length(lR))
  out <- data.frame(response=bad,rt=bad)
  out[ok,2:1] <- rdiffusion(length(lR[ok]), a = pars[ok,"a"], v = pars[ok,"v"], t0 = pars[ok,"t0"], 
    z = pars[ok,"z"], d = pars[ok,"d"], sz = pars[ok,"sz"], sv = pars[ok,"sv"],
    st0 = pars[ok,"st0"], s = pars[ok,"s"])
  # cbind.data.frame(R=factor(out[,"response"],levels=c("lower","upper"),labels=levels(lR)),rt=out[,"rt"])
  cbind.data.frame(R=factor(out[,"response"],levels=1:2,labels=levels(lR)),rt=out[,"rt"])
}

dDDM <- function(rt,R,pars,precision=3) 
  # DDM density for response factor R with rt
  # lower is mapped to first level of R and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # R <- factor(c("left","right")); rt=c(1,1)
{
  levels(R) <- c("lower","upper")
  ddiffusion(rt,response=as.character(R),
    a = pars[,"a"], v = pars[,"v"], t0 = pars[,"t0"], z = pars[,"z"],d = pars[,"d"], 
    sz = pars[,"sz"], sv = pars[,"sv"],st0 = pars[,"st0"], s = pars[,"s"])
}


pDDM <- function(rt,R,pars,precision=3) 
  # DDM cdf for response factor R with rt
  # lower is mapped to first level of R and upper to second
{
  levels(R) <- c("lower","upper")
  pdiffusion(rt,response=as.character(R),
    a = pars[,"a"], v = pars[,"v"], t0 = pars[,"t0"], z = pars[,"z"],d = pars[,"d"], 
    sz = pars[,"sz"], sv = pars[,"sv"],st0 = pars[,"st0"], s = pars[,"s"])
}


