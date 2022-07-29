rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/fMRI/Normal/normal.R")
source("models/RACE/LBA/lbaB.R")

# Load in the data
load("~/EMC/Data/VM2011/fmri_dm_VM2011.RData")
load("~/EMC/Data/VM2011/fmri_data_VM2011.RData")
load("~/EMC/Data/VM2011/behav_data_VM2011.RData")

behav_data <- behav_data[,c("subject", "cond", "stimulus", "response", "RT")]
colnames(behav_data) <- c('subjects', 'E', 'S', 'R', 'rt')

# Remove all extremely fast responses (<.15s)
behav_data <- behav_data[behav_data$rt>.15,]

# Ensure the following columns are of type factor
behav_data$subjects <- factor(behav_data$subjects)
behav_data$E <- as.factor(behav_data$E)
behav_data$S <- as.factor(behav_data$S)
behav_data$R <- as.factor(behav_data$R)

# subset data to only have one roi
data_lSTR <- fmri_data[, c("lSTR", "subjects")]
data_rSTR <- fmri_data[, c("rSTR", "subjects")]

# Make our design of the fmri data, we only need the design matrix and the data for this.
design_lSTR <- make_design_fmri(design_matrix = fmri_dm, data = data_lSTR, model = normal)
design_rSTR <- make_design_fmri(design_matrix = fmri_dm, data = data_rSTR, model = normal)

Emat <- matrix(c(-1,1),ncol=1)
dimnames(Emat) <- list(NULL,"acc-spd")     # estimate the contrast between ACC-SPD trials

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

# Make our design of the behavioral data
design_behav <- make_design(
  Ffactors=list(subjects=levels(behav_data$subjects),S=levels(behav_data$S),E=levels(behav_data$E)),
  Rlevels=levels(behav_data$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~E,A~1,t0~1),       # <-- only B affected by E
  constants=c(sv=log(1), B=log(0.05)),
  model=lbaB)

# combine all our models together
sampler <- make_samplers(list(data_lSTR, data_rSTR, behav_data), list(design_lSTR, design_rSTR, design_behav))
# we can just run sample now!
samples <- run_emc(sampler, nsample = 5000, cores_per_chain = 8, cores_for_chains = 3, verbose_run_stage = T)
save(samples, file = "joint_VM2011.RData")
plot_chains(samples)

# Here we can only examine the fit of the behavioral model (our 3rd model), so we extract it from the joint samples
samples_behav <- single_out_joint(samples, 3)
pp_behav <- post_predict(samples_behav,n_cores=1, filter = "sample", subfilter = 500)
plot_fit(behav_data,pp_behav,layout=c(2,2),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))     # aggregate over subjects and plot 'overall' CDFs

# To look at the correlations, first we merge all the samples to be one chain object
samples_merged <- merge_samples(samples)
library(corrplot)
median_corrs <- cov2cor(apply(samples_merged$samples$theta_var, 1:2, median))
corrplot(median_corrs)

# Lastly we examine the credible intervals of the correlations
CIs <- plot_density(samples,layout=c(2,4),selection="correlation", filter = "sample", do_plot=T)
CIs