rm(list = ls())
source("emc/emc.R")
source("models/RACE/LBA/lbaB.R")

load("~/Documents/UVA/2022/EMC/Data/PNAS_InAndOut.RData")
data$R <- as.factor(data$resp)
data$subjects <- data$subject
data$E <- as.factor(data$cond)
data <- data[,c("subjects", "R", "E", "scan", "rt")]
# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

# This is our design for out of scanner
dat1 <- data[data$scan == 0,]
design_ses1 <- make_design(
  Ffactors=list(subjects=unique(dat1$subjects), E=unique(dat1$E)),
  Rlevels=levels(dat1$R),
  Clist=list(R=ADmat, E=Emat),
  Flist=list(v~E*R,sv~R,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

# The design for in scanner is the same, but with different data.
dat2 <- data[data$scan == 1,]
design_ses2 <- make_design(
  Ffactors=list(subjects=unique(dat2$subjects), E=unique(dat2$E)),
  Rlevels=levels(dat2$R),
  Clist=list(R=ADmat, E=Emat),
  Flist=list(v~E*R,sv~R,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

# just make the sampler by adding our designs and data together in a list
samplers <- make_samplers(list(dat1, dat2), list(design_ses1, design_ses2),type="standard",rt_resolution=.02)

# we can just run sample now!
samples <- run_emc(samplers, nsample = 2500, cores_for_chains = 1, cores_per_chain = 20, verbose = T, verbose_run_stage = T)
save(samples, file = "joint_PNAS.Rdata")
load("joint_PNAS.Rdata")

# Post_predict will output a list of the two submodels/components in our joint model
ppJoint <- post_predict(samples,n_cores=1, filter = "sample")
# We run plot_fit for each component separately
plot_fit(dat1,ppJoint[[1]],layout=c(2,3))
plot_fit(dat2,ppJoint[[2]],layout=c(2,3))

# To look at the correlations, first we merge all the samples to be one chain object
samples_merged <- merge_samples(samples)
library(corrplot)
median_corrs <- cov2cor(apply(samples_merged$samples$theta_var, 1:2, median))
corrplot(median_corrs)
