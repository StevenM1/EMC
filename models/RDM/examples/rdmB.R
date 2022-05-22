rm(list=ls())
library(mvtnorm)
library(coda)
source("emc/emc.R")
source("models/RDM/rdmB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

# Accuracy - Speed, Neutral - Speed
Emat <- contr.sum(3)/2
dimnames(Emat) <- list(NULL,c("as","ns"))

ns <- length(levels(dat$subjects))
design <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1),
  model=rdmB)

# Design has a p_vector attribute that can be filled in with values to simulate
# data, and which itself has a "map" attribute that defines the nature of each
# parameter
attr(design,"p_vector")
#         v     v_lMd         B     B_lRd     B_Eas     B_Ens B_lRd:Eas B_lRd:Ens         A        t0 
#         0         0         0         0         0         0         0         0         0         0 
# attr(,"map")
# attr(,"map")$v
#   v v_lMd
# 1 1   0.5
# 3 1  -0.5
# 
# attr(,"map")$B
#   B B_lRd B_Eas B_Ens B_lRd:Eas B_lRd:Ens
# 1 1  -0.5   0.5   0.0     -0.25      0.00
# 3 1   0.5   0.5   0.0      0.25      0.00
# 2 1  -0.5   0.0   0.5      0.00     -0.25
# 4 1   0.5   0.0   0.5      0.00      0.25
# 
# attr(,"map")$A
# [1] 1
# 
# attr(,"map")$t0
# [1] 1

# Full (variance-covariance) fits ----

# Combine data and design, define type of PMWG sampler and rt_resolution
# By default three chains are created (n_chains). Can also define priors_list
# (if none supplied defaults are used) and extras required by certian types
# (par_groups, n_factors, contraintMat, covariates)
samplers <- make_samplers(dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRDM.RData")

# I only look at the standard number of particles here in one of my test runs
print(load("testRDM.RData"))


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
