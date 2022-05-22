rm(list=ls())
library(mvtnorm)
library(coda)
source("emc/emc.R")
source("models/LBA/lbaB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

dat$subjects <- factor(as.numeric(factor(dat$subjects)))

# # NB: This data has been truncated at 0.25s and 1.5s
# par(mfcol=c(2,3))
# sacc <- plot_defective_density(dat,correct_fun=function(dat)dat$S==dat$R,xlim=c(0,1.5))
# round(sort(sacc),2)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

# Accuracy - Speed, Neutral - Speed
Emat <- contr.sum(3)/2
dimnames(Emat) <- list(NULL,c("as","ns"))

# # Could also do: Av acc/neut vs. speed and acc vs. neut
# Emat <- matrix(c(1/4,1/4,-1/2,1/2,-1/2,0),nrow=3,dimnames=list(NULL,c("anS","an")))

design <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~lR*E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

design <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~lR*E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

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
# attr(,"map")$sv
# [1] 1
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
# save(samplers,file="testStandard.RData")


print(load("testStandard.RData"))


# How many iterations
chain_n(samples)
#      burn adapt sample
# [1,] 1471    99    500
# [2,] 1471    91    500
# [3,] 1471   101    500
     
# Lets look at some chains, default do no thinning (thin argument) and 
# dont remove early samples (subfilter argument)

# Pick which you want to look at, below I first go through burn
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

# As advertised auto_burn works, this function print multivariate values
gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
#   12    7   11    2    8    4    3   17   15    9    6   13   19    5   18    1   10   16   14 
# 1.01 1.01 1.01 1.02 1.02 1.02 1.02 1.02 1.03 1.03 1.03 1.03 1.04 1.04 1.05 1.05 1.06 1.08 1.10 

# Here it is for sample, with pdist_update_n=50 (the default in run_stages), maybe with more samples
# this would work out ... but pretty clear the adapt stage pulls chains apart so there will be
# time to settle back in and not clear what the gain is over just burn in this case ... will have
# to explore different cases and how this goes with fewer particles, the advertised advantage
# of the final stage ...
#   11    7    6    5   10    2    9    3   12   18    8   17   15    1   16    4   19   14   13 
# 1.06 1.07 1.07 1.08 1.09 1.11 1.13 1.14 1.14 1.15 1.16 1.18 1.19 1.20 1.20 1.23 1.27 1.38 1.56 

# Getting back to burn ...

# Print the save for a more detailed summary
print(round(gd,2))
#       v v_lMd    B B_lRd B_Eas B_Ens B_lRd:Eas B_lRd:Ens    A   t0 mpsrf
# 1  1.01  1.01 1.05  1.01  1.02  1.00      1.00      1.00 1.00 1.07  1.05
# 2  1.01  1.00 1.01  1.00  1.00  1.00      1.00      1.00 1.01 1.00  1.02
# 3  1.00  1.00 1.00  1.01  1.00  1.01      1.02      1.01 1.00 1.01  1.02
# 4  1.01  1.01 1.03  1.01  1.02  1.03      1.01      1.01 1.00 1.02  1.02
# 5  1.02  1.02 1.10  1.00  1.04  1.04      1.01      1.02 1.01 1.03  1.04
# 6  1.01  1.00 1.05  1.00  1.01  1.03      1.01      1.01 1.02 1.04  1.03
# 7  1.01  1.00 1.05  1.00  1.02  1.03      1.01      1.01 1.00 1.01  1.01
# 8  1.00  1.01 1.02  1.00  1.01  1.01      1.00      1.01 1.00 1.03  1.02
# 9  1.00  1.00 1.01  1.01  1.01  1.00      1.02      1.02 1.01 1.00  1.03
# 10 1.02  1.01 1.05  1.02  1.02  1.02      1.01      1.03 1.02 1.06  1.06
# 11 1.01  1.00 1.01  1.00  1.00  1.00      1.01      1.01 1.01 1.01  1.01
# 12 1.00  1.00 1.00  1.00  1.00  1.00      1.00      1.02 1.00 1.00  1.01
# 13 1.01  1.01 1.02  1.01  1.02  1.02      1.01      1.02 1.00 1.01  1.03
# 14 1.00  1.02 1.07  1.03  1.00  1.01      1.01      1.00 1.05 1.13  1.10
# 15 1.00  1.00 1.02  1.00  1.01  1.01      1.01      1.01 1.00 1.03  1.03
# 16 1.00  1.01 1.05  1.02  1.00  1.00      1.01      1.02 1.06 1.07  1.08
# 17 1.01  1.01 1.05  1.01  1.04  1.04      1.02      1.01 1.00 1.01  1.02
# 18 1.01  1.00 1.02  1.00  1.02  1.02      1.01      1.03 1.00 1.01  1.05
# 19 1.00  1.01 1.01  1.01  1.01  1.01      1.00      1.03 1.00 1.01  1.04

# Likelihoods look good
plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
# Some evidence in hypers that not quite there at start
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
# For that an neatness you might remove the first 71 as follows then re-run the above
samples <- chain_shorten(samples,71)
# Can see some evidence of one chain wandering off at the end, with really clear
# case for s4 alpha
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))
# Lets clean that up too and keep just 1000 burn
samples <- chain_shorten(samples,1001:1400)

# Sometimes shortening like this will hurt R hat, lets check (have not got around
# to implementing this for hyper yet).
gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
#   12    2   11    3    7    8   13    4    6   17    5   15    9   18   19   10    1   16   14 
# 1.01 1.02 1.02 1.02 1.02 1.02 1.03 1.03 1.04 1.05 1.05 1.05 1.06 1.06 1.06 1.08 1.10 1.11 1.19 

# Multivariate can often be sensitive like this, detailed view attributes it to
# s14 t0 wandering off a bit, it comes back so I think this is OK for now
round(gd,2)
#       v v_lMd    B B_lRd B_Eas B_Ens B_lRd:Eas B_lRd:Ens    A   t0 mpsrf
# 1  1.02  1.01 1.11  1.02  1.04  1.02      1.02      1.01 1.01 1.14  1.10
# 2  1.01  1.00 1.00  1.00  1.00  1.01      1.00      1.01 1.01 1.01  1.02
# 3  1.01  1.00 1.01  1.01  1.00  1.02      1.01      1.01 1.00 1.02  1.02
# 4  1.01  1.01 1.02  1.01  1.02  1.02      1.01      1.01 1.01 1.02  1.03
# 5  1.01  1.01 1.02  1.00  1.00  1.02      1.02      1.02 1.01 1.04  1.05
# 6  1.01  1.00 1.03  1.00  1.02  1.00      1.00      1.01 1.01 1.05  1.04
# 7  1.00  1.00 1.01  1.00  1.01  1.00      1.02      1.01 1.00 1.01  1.02
# 8  1.00  1.00 1.02  1.00  1.00  1.00      1.01      1.01 1.00 1.03  1.02
# 9  1.01  1.01 1.02  1.02  1.02  1.01      1.02      1.04 1.03 1.01  1.06
# 10 1.02  1.01 1.06  1.02  1.03  1.02      1.02      1.05 1.02 1.08  1.08
# 11 1.01  1.00 1.01  1.00  1.01  1.01      1.00      1.01 1.01 1.02  1.02
# 12 1.00  1.00 1.00  1.00  1.00  1.00      1.00      1.03 1.00 1.00  1.01
# 13 1.01  1.02 1.02  1.00  1.02  1.03      1.02      1.01 1.01 1.01  1.03
# 14 1.02  1.02 1.17  1.07  1.01  1.04      1.01      1.00 1.06 1.26  1.19
# 15 1.00  1.00 1.02  1.00  1.01  1.01      1.04      1.01 1.00 1.04  1.05
# 16 1.00  1.02 1.07  1.02  1.02  1.00      1.00      1.02 1.10 1.09  1.11
# 17 1.01  1.02 1.02  1.02  1.02  1.02      1.04      1.01 1.01 1.02  1.05
# 18 1.01  1.01 1.04  1.01  1.04  1.06      1.01      1.03 1.00 1.03  1.06
# 19 1.01  1.01 1.03  1.01  1.02  1.03      1.00      1.04 1.01 1.02  1.06


# Note you could also look at acf on a chain by chain basis, for example for chain
# one 
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
# Generally very good except s13 is a little high in places. Hypers can also be
# looked at, mu is great, variance and covariance a bit higher in places but not
# so bad. Leave that as an exercise

# Effective samples of, but probably would want more for posterior CIs
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
# Definitely need more for B interactions!
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))


# But does it fit? Niek could you take a look at the "hyper" bit of this
# function and give me some advice on how to simulate subject paraemters from
# the hyper?
pp <- post_predict(samples,filter="burn")
# Here is the average over subjects, thick lines and black points are data
# thin lines and grey points (small = 100 random posterior parameters,
# large = their average). Points are q_points = .1, .3, .5, .7, .9 quantiles 
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
# Messier for participants, missing line when less than length(q_points) data  
plot_fit(dat,pp,layout=c(2,3))
# Looks like the fitting is working but the model is a bit too simple.

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
# An attribute keeps the pars used to create
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecovery.RData")

print(load("testRecovery.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))



# Single -----------------------------------------------------------

samplers <- make_samplers(dat,design,type="single",rt_resolution=.02)
# save(samplers,file="testSingle.RData")


print(load("testSingle.RData"))
# kd6t not quite there after 100 trys
gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
# rt5t bl1t ku4t ta5t rmbt kmat scat as1t na1t kh6t zk1t vf1t rt3t hsgt bd6t kd9t hsft rt2t kd6t 
# 1.00 1.01 1.01 1.02 1.02 1.02 1.02 1.02 1.03 1.03 1.03 1.03 1.05 1.05 1.05 1.06 1.08 1.08 1.21 
print(round(gd,2))
plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="single",rt_resolution=.02)
# save(samplers,file="testRecoverySingle.RData")

print(load("testRecoverySingle.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))


##########################################  UP TO HERE 

# Factors different constraint --------------------------------------------
# With the factor group level, you can also set different constraints on the loadings matrix
# E.g. this one fixes the diagonals (which will be 1s in estimation), whereas in the default they're freely estimated.
n_factors <- 3
constraintMat <- matrix(1, nrow = length(p_vector), ncol = n_factors)
constraintMat[upper.tri(constraintMat, diag = T)] <- 0 
constraintMat

sampler <- pmwgs(
  dadm = dadm,
  n_factors = n_factors,
  constraintMat = constraintMat
)

samples <- mclapply(list(sampler,sampler, sampler),run_stages,
                    iter=c(200,NA,NA),cores_per_chain=10,cores_for_chains=3,p_accept = .7)
samples <- mclapply(samples,run_stages,
                    iter=c(NA,200,NA),cores_per_chain=10,cores_for_chains=3,p_accept = .7)
report_adapt(samples)
samples <- mclapply(samples,run_stages,
                    iter=c(NA,NA,500),cores_per_chain=10,cores_for_chains=3,p_accept = .7)
save(samples, file = "testFactorsConstrained.RData")

load("testFactorsConstrained.RData")

# Factor regression -------------------------------------------------------
source("pmwg/sampling_factorRegression.R")

# Add a covariate like age at the group level! Here I could have also set different constraint but I use the default
# I just used a random covariate
n_factors <- 2

sampler <- pmwgs(
  dadm = dadm,
  n_factors = n_factors,
  covariates = matrix(runif(ns), ncol = 1)
)

samples <- mclapply(list(sampler,sampler, sampler),run_stages,
                    iter=c(200,NA,NA),cores_per_chain=10,cores_for_chains=3, p_accept = .7)
samples <- mclapply(samples,run_stages,
                    iter=c(NA,200,NA),cores_per_chain=10,cores_for_chains=3,p_accept = .7)
report_adapt(samples)
samples <- mclapply(samples,run_stages,
                    iter=c(NA,NA,500),cores_per_chain=10,cores_for_chains=3,p_accept = .7)
save(samples, file = "testFactorRegression.RData")

load("testFactorRegression.RData")
