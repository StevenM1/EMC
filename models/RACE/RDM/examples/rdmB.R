rm(list=ls())
source("emc/emc.R")
source("models/RACE/RDM/rdmB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1),
  model=rdmB)

rdm_B <- make_samplers(dat,design_B,type="standard",rt_resolution=.02)
save(rdm_B,file="rdmPNAS_B.RData")


design_Bvt0<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~E),
  model=rdmB)

rdm_Bvt0<- make_samplers(dat,design_Bvt0,type="standard",rt_resolution=.02)
save(rdm_Bvt0,file="rdmPNAS_Bvt0.RData")

design_Bt0<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~E),
  model=rdmB)

rdm_Bt0 <- make_samplers(dat,design_Bt0,type="standard",rt_resolution=.02)
save(rdm_Bt0,file="rdmPNAS_Bt0.RData")

design_Bv<- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~1),
  model=rdmB)

rdm_Bv <- make_samplers(dat,design_Bv,type="standard",rt_resolution=.02)
save(rdm_Bv,file="rdmPNAS_Bv.RData")


print(load("rdmPNAS_B.RData"))

check_run(rdm_B)

# How many iterations
chain_n(samples)
#      burn adapt sample
# [1,]  601   110    500
# [2,]  601    83    500
# [3,]  601   105    500
     
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

# As advertised auto_burn works, this function print multivariate values
gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))

print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))


pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecoveryRDM.RData")

# run in script runRecoveryRDM
print(load("testRecoveryRDM.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))
