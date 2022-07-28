rm(list=ls())
source("emc/emc.R")

#### DDM ----
# Create Wiener data, add uniform movement time (MT) noise, then compare fit of 
# DDM allowing st0 and Wiener using t0 = Intercept + t0_MT * MT model 

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# t0 estimated on natural scale so can add MT
source("models/DDM/DDM/ddmTZDt0natural.R")

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))  

# Simulate the MT covariate as uniform noise, width st-=.12s.
Fcovariates <- list(MT=function(data) runif(dim(data)[1],-.06,.06))

# Standard 7 parameter Wiener diffusion model with an extra parameter for
# the regression slope of the MT effect on t0 with t0~MT formula which will
# create intercept (t0) and regression slope (t0_MT) parameters.
design_a_MT <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Fcovariates = c("MT"),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~MT, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(0)),
  model=ddmTZDt0natural)


#### Check model works with large single subject data set (6000 trials)
# Use Wiener model fit to get plausible parameters
print(load("models/DDM/DDM/examples/samples/sPNAS_a.RData")) 
p_vector <- sampled_p_vector(design_a_MT,doMap = FALSE)
p_vector[c(1:6,8)] <- apply(parameters_data_frame(sPNAS_a),2,mean)
p_vector[6] <- exp(p_vector[6]) # put t0 on natural scale
p_vector[7] <- 1 # Unit regression slope parameter t0_MT

# Note that make_data allows the covariates to be included, in this case
# adding uniform noise to RT.
data <- make_data(p_vector,design=design_a_MT,trials=1000,Fcovariates=Fcovariates)
plot_defective_density(data,layout=c(2,3),factors=c("S","E"))

dadmDDM <- design_model(data,design_a_MT)
# head(log_likelihood_ddm(p_vector,dadm=dadmDDM))

# Shows log-likelihood varying each parameter around true value while holding 
# others at true. Following is slow so look at pdf instead of running,
# vignettes/TrialCovariates/samples/profile.pdf
par(mfrow=c(2,4))
profile_pmwg(pname="v",p=p_vector,p_min=-1,p_max=1,dadm=dadmDDM)
profile_pmwg(pname="v_S",p=p_vector,p_min=1,p_max=3,dadm=dadmDDM)
profile_pmwg(pname="a",p=p_vector,p_min=.1,p_max=.6,dadm=dadmDDM)
profile_pmwg(pname="a_Ea-n",p=p_vector,p_min=-.25,p_max=1,dadm=dadmDDM)
profile_pmwg(pname="a_Ea-s",p=p_vector,p_min=-.25,p_max=1,dadm=dadmDDM)
profile_pmwg(pname="Z",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmDDM)
profile_pmwg(pname="t0",p=p_vector,p_min=.1,p_max=.3,dadm=dadmDDM)
profile_pmwg(pname="t0_MT",p=p_vector,p_min=.5,p_max=1.5,dadm=dadmDDM)
# # Here is the console output it makes
#      true.v         max      miss.v 
# -0.10036675 -0.11111111  0.01074436 
#     true.v_S          max     miss.v_S 
#  2.022537676  2.030303030 -0.007765354 
#      true.a         max      miss.a 
# 0.337603121 0.332323232 0.005279889 
#   true.a_Ea-n           max   miss.a_Ea-n 
#  0.0776771507  0.0782828283 -0.0006056775 
# true.a_Ea-s         max miss.a_Ea-s 
#   0.4706128   0.4823232  -0.0117104 
#       true.Z          max       miss.Z 
# -0.031940058 -0.035353535  0.003413477 
#     true.t0         max     miss.t0 
# 0.250864677 0.249494949 0.001369727 
#  true.t0_MT         max  miss.t0_MT 
# 1.000000000 0.994949495 0.005050505 

#### Do a simulation study with a data set like PNAS with 19 subjects and 
#### 900 trials each

# Extract individual parameters from Wiener diffusion fit 
ipars <- parameters_data_frame(sPNAS_a,selection="alpha")
# Calculate individual means
iMeans <- matrix(nrow=19,ncol=8,dimnames=list(levels(ipars$subjects),c(dimnames(ipars)[[2]][-1],"t0_MT")))
for (i in levels(ipars$subjects))
  iMeans[i,c(1:7)] <- apply(ipars[ipars$subjects==i,-1],2,mean)

# Draw zero mean normally distributed t0_MT values, sd = 0.2 for each participant
# (i.e., individual differences)
iMeans[,8] <- rnorm(19,0,.2); iMeans[,8] <- 1 + iMeans[,8] - mean(iMeans[,8])
iMeans[,6] <- exp(iMeans[,6]) # transform t0 to natural scale

# Make the data including uniform MT noise and save the simulation parameters
# with for a later parameter recovery check
data <- make_data(iMeans,design=design_a_MT,trials=150,Fcovariates=Fcovariates)
# Data look reasonable
plot_defective_density(data,layout=c(2,3),factors=c("S","E"))

# # Save values for later parameter recovery
# pars <- iMeans
# attr(data,"pars") <- pars
# save(pars,file="vignettes/TrialCovariates/samples/ddm_pars.RData")
print(load("vignettes/TrialCovariates/samples/ddm_pars.RData"))
head(pars)


# Create samples object and fit. Note there is no compression as MT covariate
# makes every trial unique.
ddm_a_MT <- make_samplers(data,design_a_MT,type="standard",rt_resolution=.02)
# save(ddm_a_MT,pars,file="ddm_a_MT.RData")
# fit by run_ddm_MT.R, ~3.3hr run time

# Now for a comparison fit DDM with uniform non-decision time noise, which is 
# exactly what the simulated data has, see if we can recover its width,
# st0 = 0.12. 
source("models/DDM/DDM/ddmTZD.R")
design_a_st0 <- make_design(
  Ffactors=list(subjects=levels(data$subjects),S=levels(data$S),E=levels(data$E)),
  Rlevels=levels(data$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),DP=qnorm(0.5),SZ=qnorm(0),sv=log(0)),
  model=ddmTZD)

# remove the irrelevant t0_MT parameter
pars <- pars[,-8]
attr(data,"pars") <- pars
ddm_a_st0 <- make_samplers(data,design_a_st0,type="standard",rt_resolution=.02)
# save(ddm_a_st0,pars,file="ddm_a_st0.RData")
# fit by run_ddm_st0.R, ~ 1.9hr run time

# Look at results 
print(load("vignettes/TrialCovariates/samples/ddm_a_st0.RData"))
print(load("vignettes/TrialCovariates/samples/ddm_a_MT.RData"))

# Looks good
check_run(ddm_a_st0)
check_run(ddm_a_MT)

# Slightly fewer effective parameters, MASSIVE improvement in fit for MT model
# because it captures trial-to-trial fluctuations!
compare_IC(list(st0=ddm_a_st0,MT=ddm_a_MT))
#        DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# st0  -8154    0  -8008     0        147  -8301  -8430  -8447
# MT  -11814    1 -11670     1        144 -11958 -12094 -12102


# # Large number so can average at trial level (st0 expected to be constant)
# pp_MT <- post_predict(ddm_a_MT,n_post=1000,n_cores=19)
# pp_st0 <- post_predict(ddm_a_st0,n_post=1000,n_cores=19)
# save(pp_MT,pp_st0,file="pp_MT.RData")
print(load("vignettes/TrialCovariates/samples/pp_MT.RData"))

# Standard (cell) fit virtually identical. Again a bit slow so look at pdf
# vignettes/TrialCovariates/samples/fit_MT.pdf and fit_st0.pdf
plot_fit(attr(ddm_a_st0,"data_list")[[1]],pp_st0,layout=c(2,3),factors=c("E","S"))
plot_fit(attr(ddm_a_MT,"data_list")[[1]],pp_MT,layout=c(2,3),factors=c("E","S"))

# Examine trial by trial fit, and correlations in each cell for each subject of
# observed and predicted trial sequence
# pdf("MT.pdf")
# rMT <- plot_trials(pp=pp_MT,data=attr(ddm_a_MT,"data_list")[[1]],Fcovariates="MT",layout=c(2,3))
# dev.off()

print(load("vignettes/TrialCovariates/samples/rMTst0.RData"))

# Clear and consistent correlations
round(unlist(lapply(rMT,function(x){mean(x)})),2)
# as1t bd6t bl1t hsft hsgt kd6t kd9t kh6t kmat ku4t na1t rmbt rt2t rt3t rt5t scat ta5t vf1t zk1t 
# 0.24 0.24 0.15 0.25 0.33 0.18 0.29 0.20 0.19 0.31 0.37 0.13 0.43 0.08 0.05 0.18 0.20 0.21 0.19 

# pdf("st0.pdf")
# rst0 <- plot_trials(pp=pp_st0,data=attr(ddm_a_st0,"data_list")[[1]],Fcovariates="MT",layout=c(2,3))
# dev.off()
# save(rMT,rst0,file="vignettes/TrialCovariates/samples/rMTst0.RData")

# No correlation
round(unlist(lapply(rst0,function(x){mean(x)})),2)
 # as1t  bd6t  bl1t  hsft  hsgt  kd6t  kd9t  kh6t  kmat  ku4t  na1t  rmbt  rt2t  rt3t  rt5t  scat  ta5t  vf1t  zk1t 
 # 0.01 -0.02 -0.02 -0.06 -0.03 -0.01  0.05  0.06 -0.05  0.00 -0.03  0.03 -0.05  0.03  0.03 -0.05 -0.03  0.04 -0.01 

# priors dominated in both cases
plot_density(ddm_a_MT,selection="mu",layout=c(2,4))
plot_density(ddm_a_st0,selection="mu",layout=c(2,4))

# Check random effects recovery
pars <- attributes((attributes(ddm_a_MT)$data_list[[1]]))$pars
tabs_MT <- plot_density(ddm_a_MT,pars=pars,layout=c(2,4))
pars_st0 <- pars
dimnames(pars_st0)[[2]][8] <- "st0"
pars_st0[,"st0"] <- pars_st0[,"st0"]*0.12
pars_st0[,"t0"] <- log(pars_st0[,"t0"]-pars_st0[,"st0"]/2)
pars_st0[,"st0"] <- log(pars_st0[,"st0"])
# Width of uniform noise (st0) well recovered
tabs_st0 <- plot_density(ddm_a_st0,pars=pars_st0,layout=c(2,4))

# Good recovery
plot_alpha_recovery(tabs_MT,layout=c(2,4))
plot_alpha_recovery(tabs_st0,layout=c(2,4))
