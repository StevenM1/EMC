rm(list=ls())
source("emc/emc.R")

# Data from experiment 1 of https://elifesciences.org/articles/63055 
# 55 participants performing 52 trials on each of 4 pairs of stimuli with 
# different levels of probabilistic reinforcement for or choosing one or
# other pair member. 

print(load("Data/mileticExp1.RData"))
data <- dat[!dat$excl,] # Excluded participants
data$pp <- factor(as.character(data$pp))
data <- data[!is.na(data$rt),]
# Ignore presentation side, low/high reward pairs are 
# 0.6(0.8-0.2) = tT, 0.4(0.7-0.3) = eD, 0.3(.65 - 0.35) = hJ 0.2(0.6 - 0.4) = kK
# Make S factor in order of increasing ease (.2,.3,.4,.6)
data$S <- factor(paste(data$low_stim,data$high_stim,sep="_"),levels=c("k_K","h_J","e_D","t_T"))
data <- data[,c("pp","TrialNumber","S","reward","choiceIsHighP","rt")]
data$reward <- data$reward/100
names(data) <- c("subjects","trials","S","reward","R","rt")
data$R <- factor(data$R,labels=c("low","high"))
head(data)

# Two lowest probability pairs not much different
round(tapply(data$R=="high",data$S,mean),2)

# NB: This example collapses across stimulus side, so just designates stimuli as
#     high vs. low (reward probability).  

# Racing diffusion model with reinforcement learning
source("models/RACE/RDM/rdmRL.R")

# Special component of design for adaptive models, here with only a component
# for adapting stimulus value, but architecture allows other components.
adapt <- list( 
  stimulus=list(
    # 2AFC stimulus components whose value is to be adapted
    targets=matrix(c("K","k","J","h","D","e","T","t"),nrow=2,ncol=4,
               dimnames=list(c("high","low"),c("k_K","h_J","e_D","t_T"))),
    # adaptive equation, creating "Q" values for each stimulus component
    init_par=c("q0"),
    adapt_par=c("alpha"),
    feedback=c("reward"),
    adapt_fun=function(Qlast,adapt_par,reward) # new Q value  
      Qlast + adapt_par*(reward - Qlast),
    # Output equation, uses Q values and parameters to update rates for stimulus components
    output_par_names=c("v0","w"),
    output_name = c("v"),
    output_fun=function(output_pars,Q) # rate based on Q value
            output_pars[,1] + output_pars[,2]*Q
  ),
  fixed_pars=c("B","t0","A") # Non-adapted parameters
)

# NB: Adapt design requires reward Fcovariate and Ffunctions making one column
#     for each stimulus (e.g., two columns for 2AFC) and corresponding columns
#     names (p_ + stimulus name) with the probability of getting a reward
sLow <- function(d){factor(unlist(lapply(strsplit(as.character(d$S),"_"),function(x){x[[1]]})))}
sHigh <- function(d){factor(unlist(lapply(strsplit(as.character(d$S),"_"),function(x){x[[2]]})))}
pLow <- function(d){levels(d$S) <- 1-c(.6,.65,.7,.8); as.numeric(as.character(d$S))}
pHigh <- function(d){levels(d$S) <- c(.6,.65,.7,.8); as.numeric(as.character(d$S))}

matchfun <- function(d){d$lR=="high"} # "correct" (higher expected reward) response
design <- make_design(
  Ffactors=list(subjects=levels(data$subjects),S=levels(data$S)),
  Ffunctions=list(low=sLow,high=sHigh,p_low=pLow,p_high=pHigh),
  Fcovariates = c("reward"),
  Rlevels=levels(data$R),matchfun=matchfun,
  Flist=list(v0~1,v~1,B~1,A~1,t0~1,alpha~1,w~1,q0~1),
  constants = c(A=log(0),v=log(0)),
  adapt=adapt, # Special adapt component
  model=rdmRL)


# Add covariates to data now (as it is used by adapt so cant wait until
# make_samplers adds it)
data <- add_Ffunctions(data,design)
head(data)


#### Check likelihood is correct ----

p_vector <- sampled_p_vector(design,doMap = FALSE)
# Average parameter table from Militic for Exp 1
p_vector["v0"] <- log(1.92) 
p_vector["B"] <- log(2.16)
p_vector["t0"] <- log(.1)
p_vector["w"] <- log(3.09)
p_vector["alpha"] <- log(0.12)
p_vector["q0"] <- log(.1) # was fixed to zero in estimated model


# All identical subjects
simdata <- make_data(p_vector,design=design,trials=52)

plot_defective_density(simdata,layout=c(2,2),factors=c("S"))

# Slow to compute so view pdf
# vignettes/ReinforcementLearning/samples/profile1.pdf

dadmRL <- design_model(simdata,design)
# log_likelihood_race(p_vector,dadm=dadmRL)
par(mfrow=c(2,3))
profile_pmwg(pname="v0",p=p_vector,p_min=.3,p_max=1,dadm=dadmRL)
profile_pmwg(pname="B",p=p_vector,p_min=.5,p_max=1,dadm=dadmRL)
profile_pmwg(pname="t0",p=p_vector,p_min=-3,p_max=-2,dadm=dadmRL)
profile_pmwg(pname="alpha",p=p_vector,p_min=-2.5,p_max=-1.5,dadm=dadmRL)
profile_pmwg(pname="w",p=p_vector,p_min=.75,p_max=1.5,dadm=dadmRL)
profile_pmwg(pname="q0",p=p_vector,p_min=-2.75,p_max=-2,dadm=dadmRL)
# q0 is worst estimated, in some runs down biased, in other runs up biased. 
#
#      true.v0          max      miss.v0 
#  0.652325186  0.660606061 -0.008280875 
#      true.B         max      miss.B 
# 0.770108222 0.767676768 0.002431454 
#     true.t0         max     miss.t0 
# -2.30258509 -2.32323232  0.02064723 
#  true.alpha         max  miss.alpha 
# -2.12026354 -2.10606061 -0.01420293 
#       true.w          max       miss.w 
#  1.128171091  1.136363636 -0.008192545 
#    true.q0        max    miss.q0 
# -2.3025851 -2.1666667 -0.1359184 


#  Now try a larger value of q0
p_vector["q0"] <- log(1) # was fixed to zero in estimated model
simdata1 <- make_data(p_vector,design=design,trials=52)

plot_defective_density(simdata1,layout=c(2,2),factors=c("S"))

dadmRL <- design_model(simdata1,design)
# head(log_likelihood_race(p_vector,dadm=dadmRL))
par(mfrow=c(2,3))
profile_pmwg(pname="v0",p=p_vector,p_min=.3,p_max=1,dadm=dadmRL)
profile_pmwg(pname="B",p=p_vector,p_min=.5,p_max=1,dadm=dadmRL)
profile_pmwg(pname="t0",p=p_vector,p_min=-3,p_max=-2,dadm=dadmRL)
profile_pmwg(pname="alpha",p=p_vector,p_min=-2.5,p_max=-1.5,dadm=dadmRL)
profile_pmwg(pname="w",p=p_vector,p_min=.75,p_max=1.5,dadm=dadmRL)
profile_pmwg(pname="q0",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmRL)
# q0 is now well estimated

#     true.v0          max      miss.v0 
#  0.652325186  0.653535354 -0.001210167 
#      true.B         max      miss.B 
# 0.770108222 0.767676768 0.002431454 
#     true.t0         max     miss.t0 
# -2.30258509 -2.30303030  0.00044521 
#  true.alpha         max  miss.alpha 
# -2.12026354 -2.10606061 -0.01420293 
#        true.w           max        miss.w 
#  1.1281710909  1.1287878788 -0.0006167879 
#      true.q0          max      miss.q0 
#  0.000000000 -0.005050505  0.005050505 
 
#### Fit real data

# May as well use ms resolution as no compression is possible
miletic1_rdm <- make_samplers(data,design,type="standard",rt_resolution=.001)
# save(miletic1_rdm,file="miletic1_rdm.RData")
# run by run_miletic1_rdm.R

# Sampling looks good after 500, giving 4500 good samples
print(load("vignettes/ReinforcementLearning/samples/miletic1_rdm.RData"))
check_run(miletic1_rdm,subfilter = 500)

# Fit aggregated over participants and trials fairly good, although some misfit for t_T
plot_fit(data,ppmiletic1_rdm,layout=c(2,2),factors="S")

# Priors clearly dominated
ciMiletic1_rdm <- plot_density(miletic1_rdm,layout=c(2,3),selection="mu",subfilter=500)
# Estimates fairly similar to Miletic et al.'s estimates with q0 fixed at zero.

# Plot observed rt and response probability averaged over subjects
plot_adapt(data = data,design=design,do_plot="R",ylim=c(0,1))
plot_adapt(data = data,design=design,do_plot="rt",ylim=c(0.3,1.05))

# Get adaptive estimates and fits to allow plotting of 95% CIs
# Base estimates on posterior means for each participant
pars <- parameters_data_frame(miletic1_rdm,selection="alpha",subfilter=500,stat=mean)
# # This is slow so use pre-computed. Also watch out for RAM use!
# adapted <- plot_adapt(data = attr(miletic1_rdm,"data_list")[[1]],design=design,
#   do_plot=NULL,reps=500,n_cores=11,p_vector=pars)

# Plot adaptive estimates
plot_adapt(adapted=adapted,do_plot="learn",ylim=c(0.1,.85))
plot_adapt(adapted=adapted,do_plot="par",ylim=c(2,5))
# Plot trial fits to rt and response probability averaged over subjects
plot_adapt(adapted=adapted,do_plot="R",ylim=c(0,1))
plot_adapt(adapted=adapted,do_plot="rt",ylim=c(0.45,1.15))
# There is some misfit early in learning 

# save(adapted,miletic1_rdm,ppmiletic1_rdm,
#      file="vignettes/ReinforcementLearning/samples/miletic1_rdm.RData")

#### Parameter recovery study ----

# Simulate data from the posterior means
miletic1_rdm_simdat <- post_predict(miletic1_rdm,use_par="mean",n_post=1,subfilter=500)
# Looks like the real thing.
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="R",ylim=c(0,1))
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="rt",ylim=c(0.3,1.05))
# If we provide parameters we can see the average learning and rates
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="learn",ylim=c(0,1),
           p_vector=attributes(miletic1_rdm_simdat)$pars)
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="par",ylim=c(2,5),
           p_vector=attributes(miletic1_rdm_simdat)$pars)



# can we recover these?
miletic1_rdm_sim <- make_samplers(miletic1_rdm_simdat,design,type="standard")
# save(miletic1_rdm_sim,file="miletic1_rdm_sim.RData")
# run in run_miletic1_rdm_sim.R

# Looks fine over last 1500
print(load("vignettes/ReinforcementLearning/samples/miletic1_rdm_sim.RData"))
check_run(miletic1_rdm_sim,subfilter=500)

# Excellent stationary fit
plot_fit(miletic1_rdm_simdat,ppmiletic1_rdm_sim,layout=c(2,2),factors="S")

# Some updating but likely shrinkage influential
tabs <- plot_density(miletic1_rdm_sim,selection="alpha",layout=c(2,3),
                     pars=attributes(attr(miletic1_rdm_sim,"data_list")[[1]])$pars)
# Note mapped=TRUE does not work for adaptive models.

# OK recovery, shrinkage understandable given small sample size per participant
plot_alpha_recovery(tabs,layout=c(2,3))
# Coverage is fairly decent , maybe a bit under
plot_alpha_recovery(tabs,layout=c(2,3),do_rmse=TRUE,do_coverage=TRUE)

# Plot observed rt and response probability averaged over subjects
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="R",ylim=c(0,1))
plot_adapt(data = miletic1_rdm_simdat,design=design,do_plot="rt",ylim=c(0.45,1.15))

# Get adaptive estimates and fits, which can be done by providing the samples
# object. We can provided the posterior mean parameters as before by and
# appropriate setting of p_vector (other statistics such as "median" can be 
# used) along with a subfilter if needed. Again stored to avoid slow compute
# adapted_sim <- plot_adapt(samples = miletic1_rdm_sim,do_plot=NULL,reps=500,n_cores=11,
#   p_vector="mean",subfilter=500)

# Plot adaptive estimates
plot_adapt(adapted=adapted_sim,do_plot="learn",ylim=c(0.1,.85))
plot_adapt(adapted=adapted_sim,do_plot="par",ylim=c(2,5))
# Plot trial fits to rt and response probability averaged over subjects
plot_adapt(adapted=adapted_sim,do_plot="R",ylim=c(0,1))
plot_adapt(adapted=adapted_sim,do_plot="rt",ylim=c(0.45,1.15))
# As expected when fitting the true model the misfit has disappeared.

# A better approach to fully represent uncertainty is to sample parameters
# from the posterior (one per rep) as below.
# adapted_sim_sample <- plot_adapt(samples = miletic1_rdm_sim,do_plot=NULL,reps=500,n_cores=11,
#   p_vector="sample",subfilter=500)

# Plot adaptive estimates
plot_adapt(adapted=adapted_sim_sample,do_plot="learn",ylim=c(0.1,.85))
plot_adapt(adapted=adapted_sim_sample,do_plot="par",ylim=c(2,5))
# Plot trial fits to rt and response probability averaged over subjects
plot_adapt(adapted=adapted_sim_sample,do_plot="R",ylim=c(0,1))
plot_adapt(adapted=adapted_sim_sample,do_plot="rt",ylim=c(0.45,1.15))
# In this case there is very little change to the credible intervals.


# save(adapted_sim,adapted_sim_sample,miletic1_rdm_sim,ppmiletic1_rdm_sim,
#      file="vignettes/ReinforcementLearning/samples/miletic1_rdm_sim.RData")

