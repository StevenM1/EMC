rm(list=ls())
source("emc/emc.R")
source("models/DDM/DDM/ddmTZD.R")
# NB: The "TZD" parameterization defined relative to the "rtdists" package is:
  # natural scale
  #   v = rtdists rate v (positive favors upper)
  # log scale 
  #   t0 > 0: lower bound of non-decision time 
  #   st0 > 0: rtdists width of non-decision time distribution 
  #   a > 0: rtdists upper threshold, a
  #   sv > 0: rtdists v standard deviation sv
  #   s > 0: rtdists moment-to-moment standard deviation, s
  # probit scale
  #   0 < Z < 1: rtdists start point z = Z*a 
  #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z)) 
  #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0



print(load("Data/PNAS.RData"))
# Note that this data was censored at 0.25s and 1.5s

dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# 3 x 2 design, 19 subjects
lapply(dat,levels)
# $subjects
#  [1] "as1t" "bd6t" "bl1t" "hsft" "hsgt" "kd6t" "kd9t" "kh6t"
#  [9] "kmat" "ku4t" "na1t" "rmbt" "rt2t" "rt3t" "rt5t" "scat"
# [17] "ta5t" "vf1t" "zk1t"
# 
# $E
# [1] "accuracy" "neutral"  "speed"   
# 
# $S
# [1] "left"  "right"
# 
# $R
# [1] "left"  "right"
# 
# $rt
# NULL

# 800 trials per subject
table(dat$subjects)
# as1t bd6t bl1t hsft hsgt kd6t kd9t kh6t kmat ku4t na1t rmbt rt2t 
#  810  849  843  848  849  849  849  837  846  845  848  831  842 
# rt3t rt5t scat ta5t vf1t zk1t 
#  843  691  849  845  838  806
 
# 70-90% accuracy, .35s - .5s mean RT
plot_defective_density(dat,factors=c("E","S"),layout=c(2,3))

# Test E factor with Accuracy - Neutral &  Accuracy - Speed contrasts
Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

# Test stimulus factor with intercept and right-left factor. When applied to rate 
# parameter v_S is the traditional DDM "drift rate" parameter. The 
# intercept term is drift bias. If it is zero drift rate is the same for left
# and right, if positive rate bias favors right, when negative left, e.g.,
# intercept v = 1, v_S = 2 implies an upper rate of 3 and lower rate of -1. 
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))  


# Fit a Wiener diffusion model where between-trial variability parameters are set
# to zero in "constants" (NB: this is done on the sampled scale, so qnorm (probit)
# for DP and SZ and log for st0 and sv) with the traditional characterization of
# the speed vs. accuracy emphasis factor (E) selectively influencing threshold (a), 
# and stimulus (S) affecting rate (v). In order to make the model identifiable
# we set moment-to-moment variability to the conventional value of 1. 

design_a <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(0)),
  model=ddmTZD)

# This produces a 7 parameter model
sampled_p_vector(design_a)

# v is the rate intercept, and v_S the traditional drift rate, a is the threshold
# for accuracy, a_Ea-n accuracy-neutral and a_Ea-s accuracy - neutral. t0 is the
# mean non-decision time and Z the proportional start-point bias (unbiased Z=0.5).

# We will fit a "standard" hierarchical model that allows for population 
# correlations among all parameters (7 x (7-1) /2 = 21 correlations), along with 
# population mean (mu) and variance estimates for each parameter. 

# Here we use the default hyper-prior, reproduced here explicitly (the same
# is obtained if the prior argument is omitted from make_samplers). It sets a
# mean of zero and a variance of 1 for all parameters and assumes they are 
# independent. 
prior=list(theta_mu_mean = rep(0, 7), theta_mu_var = diag(rep(5, 7)))

# make_samplers combines the data and design and makes a list of n_chains
# pmwg objects ready for sampling. Here we use the default 3 chains. Each will
# be sampled independently (multiple chains are useful for assessing convergence).
samplers <- make_samplers(dat,design_a,type="standard",
                          prior_list=prior,n_chains=3,rt_resolution=.02)
# make_samplers also compresses the data, grouping together trials with the same 
# parametersand rt's that differ less than the value specified in the rt_resolution 
# argument. Here we use the default which assumes a seconds scale (using 
# seconds is STRONGLY recommended) appropriate to visual stimuli presented on a
# monitor with 50Hz refresh (i.e., refreshing each 0.02s). Uniform random error
# in rt measurement is introduced in such a setup unless timing is coordinated
# with the sync pulse. Response devices such as the mouse or keyboard buttons
# introduce further error, so this setting is quite conservative. Taking 
# measurement resolution can substantially speed up likelihood calculation, the
# key bottleneck in sampling, in this case by a factor of 4

save(samplers,file="sPNAS_a.RData")
#  Fitting is performed in the script sPNAS_a.R (e.g., on a linux based system
# the command line is R CMD BATCH sPNAS_a.R & to run it in background)
# We run 300 burn in samples, up to 1000 adapt samples (adapt can pull out 
# early when complete, so this is an upper limit) and 1000 final "efficeint"
# samples to be used for inference, using the follwing command
sPNAS_a <- run_chains(samplers,iter=c(300,1000,1000),cores_per_chain=10)

# By default each chain gets its own cores (this can be set with the 
# cores_for_chains argument) and the cores_per_chain argument above gives each
# chain 10, so 30 are used in total (using up most of the capacity of the 32 
# core system this was run on (16 physical cores, each with two threads)

# Once sampling is completed the script also gets posterior predictive samples
# to enable model fit checks. By default this is based on randomly selecting 
# iterations from the final (sample) stage, and provides posterior predictives 
# for the random effects. Here we use one core pre participant.
ppPNAS_a <- post_predict(sPNAS_a,n_cores=19)

# Lets load in the results and look at them.
print(load("models/DDM/DDM/examples/samples/sPNAS_a.RData")) 

##### Check convergence

print(load("sPNAS_a.RData"))
sPNAS_a <- sPNAS_a_adapt
# Adapted with 533 final samples
chain_n(sPNAS_a)
plot_chains(sPNAS_a,selection="LL",layout=c(4,5),filter="burn",subfilter=100)
plot_chains(sPNAS_a,selection="epsilon",layout=c(4,5),filter="burn",subfilter=100)

par(mfrow=c(2,7))
plot_chains(sPNAS_a,selection="alpha",layout=NULL,filter="burn",subfilter=100)
plot_chains(sPNAS_a,selection="mu",layout=c(2,4),filter="burn",subfilter=100)
plot_chains(sPNAS_a,selection="variance",layout=c(2,4),filter="burn",subfilter=100)
plot_chains(sPNAS_a,selection="correlation",layout=c(3,7),filter="burn",subfilter=100)


# Focus on the sample stage from here (the default setting of filter) 

# RANDOM EFFECTS (i.e., subject level)
# Participant likelihoods all fat flat hairy caterpillars
plot_chains(sPNAS_a,selection="LL",layout=c(4,5),filter="burn")
# Plot random effects, again they look good
par(mfrow=c(2,7)) # one row per participant
plot_chains(sPNAS_a,selection="alpha",layout=NULL)

# MIXING 
# R hat indicates good mixing
gd_pmwg(sPNAS_a,selection="alpha")
# Default shows multivariate version over parameters. An invisible return
# provides full detail, here printed and rounded, all very good.
round(gd_pmwg(sPNAS_a,selection="alpha",print_summary = FALSE),2)

# SAMPLING EFFICIENCY
# Actual number samples = 3 chains x 533 per chain = 1599. 
# The average effective number shows autocorrelation is quite low
round(es_pmwg(sPNAS_a,selection="alpha"))
# Sometimes you might want to look at worst case summary
round(es_pmwg(sPNAS_a,selection="alpha",summary_alpha=min))
# To get per subject details
round(es_pmwg(sPNAS_a,selection="alpha",summary_alpha=NULL))
# Integrated autocorrelation time provides a nice summary of efficiency, a
# value of 1 means perfect efficiency, larger values indicate the approximate 
# factor by which iterations need to be increased to get a nominal value. 
iat_pmwg(sPNAS_a,selection="alpha")


# POPULATION EFFECTS
# Similar analyses as above for random effects

# Population mean
plot_chains(sPNAS_a,selection="mu",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="mu"),2) 
round(es_pmwg(sPNAS_a,selection="mu"))
iat_pmwg(sPNAS_a,selection="mu")
# Can also print autocorrelation functions for each chain (can also be done for
# alpha)
par(mfrow=c(3,7))
plot_acfs(sPNAS_a,selection="mu",layout=NULL)

# Population variance
plot_chains(sPNAS_a,selection="variance",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="variance"),2)
round(es_pmwg(sPNAS_a,selection="variance"))
iat_pmwg(sPNAS_a,selection="variance")

# There are p*(p-1)/2 correlations, where p = number of parameters.  
plot_chains(sPNAS_a,selection="correlation",layout=c(3,7),ylim=c(-1,1))
# All are estimated quite well without strong autocorrelation.
round(gd_pmwg(sPNAS_a,selection="correlation"),2)
round(es_pmwg(sPNAS_a,selection="correlation"))
iat_pmwg(sPNAS_a,selection="correlation")


########  Fit 

# Already run and saved
# ppPNAS_a <- post_predict(sPNAS_a,n_cores=5)
# save(ppPNAS_a,sPNAS_a,file="sPNAS_a.RData")

# By default the plot shows results for all subjects, putting everyone on the
# same x scale, which can make it hard to see fit for some or most subjects
plot_fit(dat,ppPNAS_a,layout=c(2,3))

# You could specify your own limits (e.g,. the data range) and move the legend)
plot_fit(dat,ppPNAS_a,layout=c(2,3),xlim=c(.25,1.5),lpos="right")

# Or plot for a single subject, e.g., the first 
plot_fit(dat,ppPNAS_a,layout=c(2,3),subject="as1t")

# This function (note the "s" in plot_fits) does the x scaling per subject
plot_fits(dat,ppPNAS_a,layout=c(2,3),lpos="right")

# Can also show the average over subjects as subjects is like any other factor,
# so just omit it. We see that the fit is OK but has various misses. 
plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))

# We can also use plot_fit to examine specific aspects of fit by supplying a 
# function that calculates a particular statistic for the data. For example to
# look at accuracy we might use (where d is a data frame containing the observed
# data or or posterior predictives) this to get percent correct
pc <- function(d) 100*mean(d$S==d$R)
# Drilling down on accuracy we see the biggest misfit is in speed for right
# responses by > 5%.
plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))

# Conversely for mean RT the biggest misses are over-estimation of RT by 50ms
# or more in the accuracy and neutral conditions
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
round(tab,2)

# However for fast responses (10th percentile) there is global under-prediction.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.375),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
round(tab,2)

# and for slow responses (90th percentile) global over-prediction .
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
round(tab,2)

# This corresponds to global under-estimation of RT variability.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
round(tab,2)


# Errors tend to be much faster than the models for all conditions.
tab <- plot_fit(dat,ppPNAS_a,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,2)


######## Posterior parameter inference ----

# Population means (mu)
# Here we use plot_density without plotting to get a table of 95% parameter CIs
ciPNAS <- plot_density(sPNAS_a,layout=c(2,4),selection="mu",do_plot=FALSE)
round(ciPNAS,2)

# Lets look at the mapping using the mean of the population mean (mu)
# samples, which is stored as an attribute (very little difference from
# the median, ciPNAS[2,])
pMU <- attr(ciPNAS,"mean")
# Now map the mean into the fill design. Between trial variability parameters
# are fixed at either zero (), 
mapped_par(pMU,design_a)
#       S        E      v     a   t0     Z
# 1  left accuracy -2.338 1.502 0.25 0.487
# 2  left  neutral -2.338 1.384 0.25 0.487
# 3  left    speed -2.338 0.907 0.25 0.487
# 4 right accuracy  2.119 1.502 0.25 0.487
# 5 right  neutral  2.119 1.384 0.25 0.487
# 6 right    speed  2.119 0.907 0.25 0.487

# Slight bias to left (lower)
p_test(sPNAS_a,p_name="v",x_selection = "mu", x_filter="sample")

# Non-decision time tightly estimated
p_test(sPNAS_a,p_name="t0",x_selection = "mu", x_filter="sample")

# Threshold less precise
p_test(sPNAS_a,p_name="a",x_selection = "mu", x_filter="sample")

# Very slight left bias
p_test(sPNAS_a,p_name="Z",x_selection = "mu", x_filter="sample",
       mu=0.5,digits=3,alternative = "greater")


######### E Main Effect ----

# Accuracy threshold is credibly but not much greater than neutral
p_test(sPNAS_a,p_name="a_Ea-n",x_selection = "mu", x_filter="sample")

# Accuracy threshold is much greater than speed.
p_test(sPNAS_a,p_name="a_Ea-s",x_selection = "mu", x_filter="sample")


######### S Main effect ----

# Average rate is ~2.4
p_test(sPNAS_a,p_name="v_S",x_selection = "mu", x_filter="sample")


#####  FULL DDM

design_a_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

samplers <- make_samplers(dat,design_a_full,type="standard")
save(samplers,file="sPNAS_a_full.RData")



