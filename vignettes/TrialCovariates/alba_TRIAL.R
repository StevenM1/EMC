rm(list=ls())
source("emc/emc.R")

#### ALBA ----

source("models/RACE/LBA/albaB.R") # So far only advantage coding for binary case


print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
dat$b_L_  <- 1+as.numeric(dat$S=="left")
dat$b_R <- 1+as.numeric(dat$S=="right")
  head(dat)

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  
Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

design_bLR <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Fcovariates = c("b_L","b_R"),
  Ffunctions = list(SS=function(d)d$b_L+d$b_R,
    SD=function(d)ifelse(d$lR=="left",d$b_L-d$b_R,d$b_R-d$b_L)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v_0 ~ 1, v_D~1, v_S~1, sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1),v_S=1/3),
  model=albaB)
   # v_0    v_D      B B_Ea-n B_Ea-s      A     t0 
   #   0      0      0      0      0      0      0 
# NB: set v_S to constant above so can compare with standard LBA (same number
# of parameters), set 1/3 so match = 3, mismatch = 1



source("models/RACE/LBA/lbaB.R")
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)
     # v  v_lMd      B B_Ea-n B_Ea-s      A     t0 
     # 0      0      0      0      0      0      0 

#### Check model works with one subject, ~ 6000 trials
design_bLR1 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects)[[1]],S=levels(dat$S),E=levels(dat$E)),
  Ffunctions = list(
    b_L=function(d){1+as.numeric(d$S=="left")},   # NB: Need to order so b_L and
    b_R=function(d){1+as.numeric(d$S=="right")},  #     b_R are made first so can
    SS=function(d)d$b_L+d$b_R,                    #     be used in SS and SD
    SD=function(d)ifelse(d$lR=="left",d$b_L-d$b_R,d$b_R-d$b_L)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v_0~1, v_D~1, v_S~1, sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1),v_S=1/3),
  model=albaB)

# Use LBA model fit to get plausible parameters
print(load("models/RACE/LBA/examples/samples/sPNAS_B.RData"))
mpar <- apply(parameters_data_frame(sPNAS_B),2,mean)
p_vector <- sampled_p_vector(design_bLR,doMap = FALSE)
p_vector[1:2] <- c(1,1); names(p_vector)[1:2] <- c("v_0","v_D") 
p_vector[3:7] <- mpar[3:7]

# Note that here we use the real data to specify trial numbers etc.
data <- make_data(p_vector,design=design_bLR1,trials=1000,
  Fcovariates=data.frame(b_L=rep(c(2,1),times=3000),b_R=rep(c(1,1),times=3000)))
plot_defective_density(data,layout=c(2,3),factors=c("S","E"),xlim=c(0,1.5))

dadmALBA <- design_model(data,design_bLR1)
head(dadmALBA)
# log_likelihood_ddm(p_vector,dadm=dadmALBA)

# Shows log-likelihood varying each parameter around true value while holding 
# others at true. Following is slow so look at pdf instead of running,
# vignettes/TrialCovariates/samples/profile.pdf
par(mfrow=c(2,4))
profile_pmwg(pname="v_0",p=p_vector,p_min=0,p_max=2,dadm=dadmALBA)
profile_pmwg(pname="v_D",p=p_vector,p_min=0,p_max=2,dadm=dadmALBA)
profile_pmwg(pname="B",p=p_vector,p_min=-1,p_max=0,dadm=dadmALBA)
profile_pmwg(pname="B_Ea-n",p=p_vector,p_min=-.25,p_max=.75,dadm=dadmALBA)
profile_pmwg(pname="B_Ea-s",p=p_vector,p_min=.25,p_max=1.25,dadm=dadmALBA)
profile_pmwg(pname="A",p=p_vector,p_min=-.75,p_max=0,dadm=dadmALBA)
profile_pmwg(pname="t0",p=p_vector,p_min=-2,p_max=-1.5,dadm=dadmALBA)

#### Do a simulation study with a data set like PNAS with 19 subjects and 
#### 900 trials each

# Extract individual parameters from LBA fit 
ipars <- parameters_data_frame(sPNAS_B,selection="alpha")
# Calculate individual means
iMeans <- matrix(nrow=19,ncol=7,dimnames=list(levels(ipars$subjects),dimnames(ipars)[[2]][-1]))
for (i in levels(ipars$subjects))
  iMeans[i,c(3:7)] <- apply(ipars[ipars$subjects==i,-c(1:3)],2,mean)
colnames(iMeans)[1:2] <- c("v_0","v_D")

# Draw zero mean normally distributed values, sd = 0.2 for each participant
# (i.e., individual differences)
iMeans[,1] <- rnorm(19,1,.2); 
iMeans[,2] <- rnorm(19,1,.2)

# To make life easy auto-calculate brightness
design_bLRa <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Ffunctions = list(
    b_L=function(d){1+as.numeric(d$S=="left")},   # NB: Need to order so b_L and
    b_R=function(d){1+as.numeric(d$S=="right")},  #     b_R are made first so can
    SS=function(d)d$b_L+d$b_R,                    #     be used in SS and SD
    SD=function(d)ifelse(d$lR=="left",d$b_L-d$b_R,d$b_R-d$b_L)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v_0~1, v_D~1, v_S~1, sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1),v_S=1/3),
  model=albaB)

data <- make_data(iMeans,design=design_bLRa,trials=150)
# Data look reasonable
plot_defective_density(data,layout=c(2,3),factors=c("S","E"),xlim=c(0,2),mfcol=FALSE)

# Save values for later parameter recovery
pars <- iMeans
attr(data,"pars") <- pars
save(pars,file="vignettes/TrialCovariates/samples/lba_pars.RData")


# Create samples object and fit. Note there is no compression as MT covariate
# makes every trial unique.
lba_bLR <- make_samplers(data,design_bLRa,type="standard",rt_resolution=.02)
# save(lba_bLR,pars,file="lba_bLR.RData")
# fit by run_lba_bLR.R

# Now fit standard LBA for comparison
lba_B <- make_samplers(data,design_B,type="standard",rt_resolution=.02)
# save(lba_B,pars,file="lba_B.RData")
# fit by run_lba_B.R

# Look at results 
load("vignettes/TrialCovariates/samples/lba_B.RData")
load("vignettes/TrialCovariates/samples/lba_bLR.RData")

# Looks good
check_run(lba_bLR)
check_run(lba_B)


# pplba_bLR <- post_predict(lba_bLR,n_cores=19)
# save(lba_bLR,pplba_bLR,file="lba_bLR.RData")
# pplba_B <- post_predict(lba_B,n_cores=19)
# save(lba_B,pplba_B,file="lba_B.RData")


# Identical as expected
compare_IC(list(lba=lba_B,bLR=lba_bLR))
#        DIC  wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# lba -14083 0.438 -13954  0.41        129 -14213 -14312 -14342
# bLR -14084 0.562 -13955  0.59        129 -14213 -14313 -14342

# Very good fits
plot_fit(attr(lba_bLR,"data_list")[[1]],pplba_bLR,layout=c(2,3),factors=c("E","S"))
plot_fit(attr(lba_B,"data_list")[[1]],pplba_B,layout=c(2,3),factors=c("E","S"))


# priors dominated in both cases
plot_density(lba_bLR,selection="mu",layout=c(2,4))
plot_density(lba_B,selection="mu",layout=c(2,4))

# Check random effects recovery
pars_bLR <- attributes((attributes(lba_bLR)$data_list[[1]]))$pars
tabs_bLR <- plot_density(lba_bLR,pars=pars_bLR,layout=c(2,4))
pars_B <- attributes((attributes(lba_B)$data_list[[1]]))$pars
dimnames(pars_B)[[2]][1:2] <- c("v","v_lMd")
pars_B[,"v"] <- 1 + pars_bLR[,"v_0"]
pars_B[,"v_lMd"] <- 2*pars_bLR[,"v_D"]
tabs_B <- plot_density(lba_B,pars=pars_B,layout=c(2,4))

# Good recovery
plot_alpha_recovery(tabs_bLR,layout=c(2,4))
plot_alpha_recovery(tabs_B,layout=c(2,4))
