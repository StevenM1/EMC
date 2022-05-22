rm(list=ls())
library(mvtnorm)
library(coda)
source("emc/emc.R")
source("models/DDM/ddmTZD.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and half difference = (v - -v)/2 = v matrix
ADmat <- matrix(c(-1/4,1/4),ncol=1,dimnames=list(NULL,"d"))  

# Accuracy - Speed, Neutral - Speed
Emat <- contr.sum(3)/2
dimnames(Emat) <- list(NULL,c("as","ns"))

ns <- length(levels(dat$subjects))
design <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=ADmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~E, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0)),
  model=ddmTZD)

# Design has a p_vector attribute that can be filled in with values to simulate
# data, and which itself has a "map" attribute that defines the nature of each
# parameter
# attr(design,"p_vector")
# #     v  v_Sd     a a_Eas a_Ens    sv    t0     Z Z_Eas Z_Ens 
# #     0     0     0     0     0     0     0     0     0     0 
# # attr(,"map")
# # attr(,"map")$v
# #   v  v_Sd
# # 1 1 -0.25
# # 2 1  0.25
# # 
# # attr(,"map")$a
# #   a a_Eas a_Ens
# # 1 1   0.5   0.0
# # 2 1   0.0   0.5
# # 
# # attr(,"map")$sv
# # [1] 1
# # 
# # attr(,"map")$t0
# # [1] 1
# # 
# # attr(,"map")$st0
# # [1] 1
# # 
# # attr(,"map")$s
# # [1] 1
# # 
# # attr(,"map")$Z
# #   Z Z_Eas Z_Ens
# # 1 1   0.5   0.0
# # 2 1   0.0   0.5
# # 
# # attr(,"map")$SZ
# # [1] 1
# # 
# # attr(,"map")$DP
# # [1] 1

# p_vector <- attr(design,"p_vector")
# p_vector[1:2] <- c(0,2)
# p_vector[3:7] <- log(c(2,1,1,1,.3))
# p_vector[8:10] <- qnorm(c(0.5,.5,.5))
# dataDDM <- make_data(p_vector,design=design,model=ddmTZD,trials=1000)
# plot_defective_density(dataDDM,layout=c(1,2),xlim=c(0,3))
# dadmDDM <- design_model(data=dat,design=design,model=ddmTZD)
# head(ddmTZD$log_likelihood(p_vector,dadmDDM))
 
# Full (variance-covariance) fits ----

# Combine data and design, define type of PMWG sampler and rt_resolution
# By default three chains are created (n_chains). Can also define priors_list
# (if none supplied defaults are used) and extras required by certian types
# (par_groups, n_factors, contraintMat, covariates)
samplers <- make_samplers(dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testDDM.RData")

# I only look at the standard number of particles here in one of my test runs
print(load("testDDM.RData"))
samples=samples1

# How many iterations
chain_n(samples)

# Pick which you want to look at, below I first go through burn
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

# na1t not there after 100, due to sv
gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
# kd6t ku4t rt5t rmbt as1t bd6t kmat zk1t kd9t vf1t rt3t scat ta5t bl1t kh6t hsgt rt2t hsft na1t 
# 1.02 1.02 1.03 1.03 1.04 1.05 1.06 1.06 1.06 1.07 1.07 1.08 1.08 1.09 1.09 1.10 1.11 1.11 1.19 

print(round(gd,2))
#         v v_Sd    a a_Eas a_Ens   sv t0    Z Z_Eas Z_Ens mpsrf
# as1t 1.01 1.00 1.00  1.00  1.01 1.05  1 1.00  1.00  1.01  1.04
# bd6t 1.00 1.00 1.00  1.01  1.00 1.06  1 1.00  1.01  1.01  1.05
# bl1t 1.00 1.00 1.01  1.01  1.01 1.14  1 1.00  1.01  1.01  1.09
# hsft 1.00 1.00 1.00  1.00  1.00 1.16  1 1.00  1.00  1.00  1.11
# hsgt 1.00 1.00 1.00  1.00  1.00 1.13  1 1.00  1.00  1.00  1.10
# kd6t 1.00 1.00 1.01  1.00  1.01 1.05  1 1.01  1.00  1.01  1.02
# kd9t 1.00 1.00 1.00  1.00  1.00 1.07  1 1.00  1.01  1.01  1.06
# kh6t 1.00 1.00 1.00  1.00  1.00 1.15  1 1.00  1.01  1.00  1.09
# kmat 1.00 1.00 1.00  1.00  1.00 1.08  1 1.00  1.00  1.00  1.06
# ku4t 1.00 1.00 1.01  1.00  1.00 1.02  1 1.00  1.01  1.01  1.02
# na1t 1.00 1.00 1.00  1.00  1.00 1.29  1 1.00  1.01  1.00  1.19
# rmbt 1.00 1.00 1.00  1.00  1.00 1.05  1 1.00  1.01  1.00  1.03
# rt2t 1.01 1.00 1.00  1.00  1.00 1.15  1 1.00  1.00  1.00  1.11
# rt3t 1.00 1.00 1.00  1.00  1.01 1.09  1 1.00  1.01  1.01  1.07
# rt5t 1.01 1.00 1.00  1.00  1.00 1.06  1 1.01  1.01  1.00  1.03
# scat 1.00 1.01 1.01  1.00  1.00 1.09  1 1.01  1.01  1.00  1.08
# ta5t 1.00 1.00 1.00  1.00  1.00 1.11  1 1.00  1.01  1.01  1.08
# vf1t 1.00 1.00 1.00  1.00  1.01 1.07  1 1.00  1.01  1.01  1.07
# zk1t 1.00 1.00 1.00  1.00  1.00 1.09  1 1.00  1.00  1.01  1.06

# Likelihoods look good
plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))

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
# An attribute keeps the pars used to create
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecoveryDDM.RData")

print(load("testRecoveryDDM.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))
