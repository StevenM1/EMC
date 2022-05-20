rm(list=ls())
library(coda)
source("emc/emc.R")
source("models/DDM/ddmTZD.R")

load("Data/PNAS.RData")
# Note that this data was censored at 0.25s and 1.5s


dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
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

table(dat$subjects)
# as1t bd6t bl1t hsft hsgt kd6t kd9t kh6t kmat ku4t na1t rmbt rt2t 
#  810  849  843  848  849  849  849  837  846  845  848  831  842 
# rt3t rt5t scat ta5t vf1t zk1t 
#  843  691  849  845  838  806
 
plot_defective_density(dat,factors=c("E","S"),layout=c(2,3))

# Accuracy - Neutral, Accuracy - Speed
Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))  

ns <- length(levels(dat$subjects))
design_a <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Clist=list(S=Vmat,E=Emat),
  Flist=list(v~S,a~E,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(1)),
  model=ddmTZD)

sampled_p_vector(design_a)
     # v    v_S      a a_Ea-n a_Ea-s     t0      Z 
     # 0      0      0      0      0      0      0 

samplers <- make_samplers(dat,design_a,type="standard",rt_resolution=.02,
  prior=list(theta_mu_mean = rep(0, 7), theta_mu_var = diag(rep(5, 7))))
# Likelihood speedup factor: 4

sPNAS_a <- auto_burn(samplers,burn=FALSE,ndiscard=100,nstart=100,min_es=500,
                     cores_per_chain=4)

##### Check convergence

chain_n(sPNAS_a)

plotChains(sPNAS_a,selection="LL",filter="sample",layout=c(4,5))
plotACFs(sPNAS_a,selection="LL",filter="sample",layout=c(4,5))
  
par(mfrow=c(2,7))
plotChains(sPNAS_a,selection="alpha",filter="sample",layout=NULL)
# By default plots for for first subject
plotACFs(sPNAS_a,selection="alpha",filter="sample",layout=NULL)
# Can select to subject by number or name, this is the same as subject=1
plotACFs(sPNAS_a,selection="alpha",filter="sample",layout=NULL,subject="as1t")

round(gd_pmwg(sPNAS_a,selection="alpha",filter="sample"))
round(es_pmwg(sPNAS_a,selection="alpha",filter="sample"))

plotChains(sPNAS_a,selection="mu",filter="sample",layout=c(2,4))
plotACFs(sPNAS_a,selection="mu",filter="sample",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="mu",filter="sample"))
round(es_pmwg(sPNAS_a,selection="mu",filter="sample"))

plotChains(sPNAS_a,selection="variance",filter="sample",layout=c(2,4))
plotACFs(sPNAS_a,selection="variance",filter="sample",layout=c(2,4))
round(gd_pmwg(sPNAS_a,selection="variance",filter="sample"))
round(es_pmwg(sPNAS_a,selection="variance",filter="sample"))

plotChains(sPNAS_a,selection="covariance",filter="sample",layout=c(3,7))
plotACFs(sPNAS_a,selection="covariance",filter="sample",layout=c(3,7))
round(gd_pmwg(sPNAS_a,selection="covariance",filter="sample"))
round(es_pmwg(sPNAS_a,selection="covariance",filter="sample"))


########  Fit 

ppPNAS_a <- post_predict(sPNAS_a,filter="sample",n_cores=5)
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

ciPNAS <- plotDensity(sPNAS_a,layout=c(2,4),do_plot=FALSE,
                      selection="mu",filter="sample")
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




