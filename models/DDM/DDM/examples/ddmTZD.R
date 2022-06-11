rm(list=ls())
source("emc/emc.R")
source("models/DDM/DDM/ddmTZD.R")

load("Data/PNAS.RData")
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
# parameter it produces the traditional DDM "drift rate" parameter. The 
# intercept term is drift bias. If it is zero drift rate is the same for left
# and right. 
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))  


# Fit a Wiener diffusion model (i.e., no between-trial variability) with
# traditional characterization where E affects only threshold (a) 

design_a <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(1)),
  model=ddmTZD)

# 7 parameter model
sampled_p_vector(design_a)

samplers <- make_samplers(dat,design_a,type="standard",rt_resolution=.02,
  prior=list(theta_mu_mean = rep(0, 7), theta_mu_var = diag(rep(5, 7))))
# Likelihood speedup factor: 4
# save(samplers,file="sPNAS_a.RData")
#  Fit in script sPNAS_a.R 
print(load("models/DDM/DDM/examples/samples/sPNAS_a.RData")) 

##### Check convergence
# Adapted with 533 final samples
chain_n(sPNAS_a)

# Focus on the sample stage from here (the default setting of filter) 

# RANDOM EFFECTS (i.e., subject level)
# Participant likelihoods all fat flat hairy caterpillars
plot_chains(sPNAS_a,selection="LL",layout=c(4,5))
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

ppPNAS_a <- post_predict(sPNAS_a,n_cores=5)
# save(ppPNAS_a,sPNAS_a,file="sPNAS_a.RData")

# By default the plot shows results for all subjects, putting everyone on the
# same x scale, which can make it hard to see fit for some or most subjects
plot_fit(dat,ppPNAS_a,layout=c(2,3))

# You could specify your own limits (e.g,. the data range) and move the 
# legend)
plot_fit(dat,ppPNAS_a,layout=c(2,3),xlim=c(.25,1.5),lpos="right")

# Or plot for a single subject, e.g., the first 
plot_fit(dat,ppPNAS_a,layout=c(2,3),subject="as1t")

# This function (note the s on fits) does the x scaling per subject
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

ciPNAS <- plot_density(sPNAS_a,layout=c(2,4),do_plot=FALSE,selection="mu")
round(ciPNAS,2)
#           v  v_S    a a_Ea-n a_Ea-s    t0     Z
# 2.5%  -0.26 1.69 0.31   0.04   0.38 -1.42 -0.07
# 50%   -0.11 2.24 0.41   0.08   0.51 -1.39 -0.03
# 97.5%  0.05 2.76 0.50   0.12   0.62 -1.35  0.00

# Lets look at the mapping using the mean of the population mean (mu)
# samples, which is stored as an attribute (very little difference from
# the median, ciPNAS[2,])
pMU <- attr(ciPNAS,"mean")
# Just take first 6 rows (full output repeats for all 19 subjects)
# NB: The correspondence to population mean  
mapped_par(pMU,design_a)[1:6,c(1:4,6,9)]
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




