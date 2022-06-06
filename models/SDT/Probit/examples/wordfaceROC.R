rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")
# This matchfun useful for SDT models, assumes binary noise/signal factor and 
# an even number of confidence ratings
matchfun <- function(d) as.numeric(d$S) == (as.numeric(d$lR)>length(levels(d$lR))/2)

print(load("Data/wordfaceROC.RData"))
# remove RT for ROC fitting
wordfaceROC$rt <- NA

#### Test fits ----
wordfaceROC1 <- wordfaceROC[wordfaceROC$subjects=="104",]
wordfaceROC1$subjects <- factor(as.character(wordfaceROC1$subjects))

# To allow different thresholds for words and faces we must nest the lR factor
# within the FW factor (FW/lR). This ensures that thresholds are always 
# increasing within each level of FW. FW*lR will not enforce the increasing 
# constraint (but FW + lR will). If you wanted to allow different but increasing
# thresholds across combinations of factors, say A and B, use (A*B)/lR etc.
designFW1 <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=104,S=c("new","old"),FW=c("faces","words")),Rlevels=1:6, 
  matchfun=matchfun,
  constants=c(mean=0,mean_FWwords=0,sd=0,sd_FWwords=0),model=probit)

p_vector <- sampled_p_vector(designFW1)
# constant mean new face = 0 
p_vector[1:2] <- c(1,1) # old faces, old words 
# constant new face sd = 1
p_vector[3:4] <- log(c(1.25,1)) # old sd = 1.25, no item effect on sd 
p_vector[5:6] <- c(-0.5,-0.5)  # first threshold for faces, shift down by 0.5 for words
p_vector[7:14] <- log(rep(c(.5,.75),4)) # .5 threshold spacing for faces, 0.75 spacing for words

# Check mapping
mapped_par(p_vector,designFW1) 

# Simulate data
simFW1 <- make_data(p_vector,design=designFW1,trials=10000)
par(mfrow=c(2,2))
plot_roc(simFW1[simFW1$FW=="faces",],main="Faces")
plot_roc(simFW1[simFW1$FW=="faces",],zROC=TRUE,qfun=qnorm,main="Faces",lim=c(-2.5,2))
plot_roc(simFW1[simFW1$FW=="words",],main="Words")
plot_roc(simFW1[simFW1$FW=="words",],zROC=TRUE,qfun=qnorm,main="Words",lim=c(-2.5,2))

# profiles
dadmFWsim <- design_model(data=simFW1,design=designFW1)
par(mfrow=c(2,8))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.25,p_max=p_vector[i]+.25,dadm=dadmFWsim))

samplers <- make_samplers(simFW1,designFW1,type="single")
save(samplers,file="probitFWsim.RData")
# runSingleProbitFWsim.R to get 1000 samples

samplers <- make_samplers(wordfaceROC1,designFW1,type="single")
# save(samplers,file="probitFW1.RData")
# runSingleProbitFW1.R to get 1000 samples

# Look at simulation, clearly not close at 1000, added another 2000
print(load("probitFWsim.RData"))
plot_chains(samples,subfilter=200,layout=c(4,4))  # Fully converged
gd_pmwg(samples,subfilter=200)
plotACFs(samples,subfilter=200,layout=c(4,4)) # Quite autocorrelated
iat_pmwg(samples,subfilter=200) 
# Excellent recovery
tabs <- plot_density(samples,subfilter=200,layout=c(4,4),pars=p_vector)

# Look at data, quick convergence 
print(load("probitFW1.RData")) 
plot_chains(samples,subfilter=100,layout=c(4,4)) # very fast convergence
gd_pmwg(samples,subfilter=100) # 1.01
plotACFs(samples,subfilter=100,layout=c(4,4)) # Quite autocorrelated
iat_pmwg(samples,subfilter=100) 
tabs <- plot_density(samples,subfilter=100,layout=c(4,4)) # prior domination ok

#### Fit full data set ----

designFW <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=levels(wordfaceROC$subjects),S=levels(wordfaceROC$S),
                FW=levels(wordfaceROC$FW)),Rlevels=1:6, 
  matchfun=matchfun,
  constants=c(mean=0,mean_FWwords=0,sd=0,sd_FWwords=0),model=probit)

samplers <- make_samplers(wordfaceROC,designFW,type="standard")
# save(samplers,file="probitFW.RData")
# runProbitFW.R to get 2000 burn in samples


print(load("probitFW.RData")) 
# All converged 
plot_chains(samples,filter="burn",subfilter=100,layout=c(3,5),selection="mu") 
plot_chains(samples,filter="burn",subfilter=100,layout=c(3,5),selection="variance") 
plot_chains(samples,filter="burn",subfilter=100,layout=c(4,4),selection="correlation") 
plot_chains(samples,filter="burn",subfilter=100,layout=c(3,5),selection="alpha") 
check_adapt(samples)
# Chain 1 adapted by iteration 292
# Chain 2 adapted by iteration 299
# Chain 3 adapted by iteration 301

# Converged, some correlations and alphas look quite autocorrelated
plot_chains(samples,filter="sample",layout=c(3,5),selection="mu",thin=10) 
plot_chains(samples,filter="sample",layout=c(3,5),selection="variance",thin=10) 
plot_chains(samples,filter="sample",layout=c(4,4),selection="correlation",thin=10) 
plot_chains(samples,filter="sample",layout=c(3,5),selection="alpha",thin=10) 

selection="mu"; selection="variance"; selection="alpha"; layout=c(3,5)
selection="correlation"; layout=c(4,7)
gd_pmwg(samples,filter="sample",selection=selection)
iat_pmwg(samples,filter="sample",selection=selection,summary_alpha=max) 
round(es_pmwg(samples,filter="sample",selection=selection,summary_alpha=min))
tabs <- plot_density(samples,filter="sample",layout=layout,selection=selection)

#### Fit
# ppWordFace <- post_predict(samples,filter="sample",n_cores=18)
# save(ppWordFace,file="ppWordFace.RData")
load("ppWordFace.RData")
# For type=SDT plot_fit requires a factor (by default "S", argument signalFactor) 
# whose first level is noise and second level is signal in order to construct an 
# ROC. Where there this 2 level structure does not apply (e.g., different types 
# of signal) subset the factor. 

# Average fit
par(mfrow=c(2,2))
plot_fit(wordfaceROC,ppWordFace,factors=c("FW","S"))
plot_fit(wordfaceROC,ppWordFace,factors=c("FW","S"),zROC=TRUE,qfun=qnorm,lim=c(-1.5,1.25))

# Individual fits
par(mfrow=c(2,4))
plot_fit(wordfaceROC,ppWordFace,zROC=TRUE,qfun=qnorm)


### Look at parameter estimates

# We can look at parameters in two ways, 1) in terms of the parameters that were
# actually sampled (and hence always transformed to have no bounds) with names
# organized by type (for the probit model the types are "mean", "sd" and 
# "threshold") that are shown by this function:
sp_names <- p_names(samples); sp_names

# 2) in terms of the parameters that are mapped to the model parameterization,
# which may be bounded, and which correspond to the combinations of factors that
# specified in the formula for the parameter type they correspond to:
mp_names <- p_names(samples,mapped=TRUE); mp_names

# Note that in this form some of the parameters can be constants. The derivation
# of their names can be seen by looking at the design, being the unique cells
# for the factors in the formula for each type
mp_design <- p_names(samples,mapped=TRUE,design=TRUE); mp_design

# Lets look at sampled parameters first

# At hyper level the posterior (black lines) shows moderately strong domination 
# of the prior (red lines).
tab_mu <- plot_density(samples,filter="sample",selection="mu",layout=c(2,7))

# Start by looking at mean and sd parameters
round(tab_mu[,sp_names$mean],2)
round(tab_mu[,sp_names$sd],2)

# To understand these, recall that mean for new faces and words (mean & 
# mean_FWwords) are set to zero and corresponding sd (sd sd_FWwords) set to 1
# (so as sd is sampled on log scale sampled values are set to zero).

# Because we did not specify contrasts R's default treatment code was used.
# We can get a list of the contrasts for each parameter type and look at the
# mean and sd design matrices (which are the same as they share formulas).
maps <- get_map(samples)
maps$mean 
# Reading across row 3 (down column Hence, mean_Sold is the study effect for 
# faces (0.69) and mean_FWwords:Sold is the extra study effect for words 
# relative to faces (1.05).
maps$sd
# The same is true for sd (0.18 and 0.22 respectively), but recall the effect is 
# on the log scale.


tab_mu_mapped <- plot_density(samples,filter="sample",selection="mu",layout=c(2,7),mapped=TRUE)

# Need to remove cells set to constant in mp_names to look at estimates.

# For means we can see that mean_words_old = mean_Sold + mean_FWwords:Sold
round(tab_mu_mapped[,mp_names$mean[-c(1:2)]],2)

# Similarly sd_words_old = exp(sd_Sold) and sd_words_old = exp(sd_Sold + sd_FWwords:Sold)
round(tab_mu_mapped[,mp_names$sd[-c(1:2)]],2)


#### Testing parameter estimates ----

# Hence mean_FWwords:Sold > 0 :Sold tests if d'(words) > d'(faces). 
# Similarly for sd_FWwords:Sold. We can see this is the case from tab_mu above, 
# as the 95% CI is well away from zero, but we can also see the same thing by 
# running a single sample test, which by default compares to mu=0.
p_test(samples,p_name="sd_FWwords:Sold",x_selection = "mu")
# Can also compare to another values, say 1 for the mean effect.
p_test(samples,p_name="mean_FWwords:Sold",x_selection = "mu",mu=1)

# We can also plot and test individual subjects. If the subject argument is 
# omitted all subjects are plotted, whereas for the test the first subject is selected  
plot_density(samples,filter="sample",selection="alpha",layout=c(2,7),subject="104")
p_test(samples,p_name="sd_FWwords:Sold",x_selection = "alpha",x_name="104")



# Turning to thresholds we get estimates for the first face and word threshold
round(tab_mu[,5:6],2)

# To test if they credibly differ:
p_test(samples,p_name="sd_FWwords:Sold",x_selection = "mu")


# We could also look at variances. Again focusing on 
var_tab <- plot_density(samples,filter="sample",selection="variance",layout=c(2,7))
round(var_tab[,1:4],2)
# Sometimes interest focuses on the "mapped" parameters
tab_mu_mapped <- plotDensity(samples,filter="sample",selection="mu",mapped=TRUE,layout=c(2,7))
round(tab_mu_mapped[,1:4],2)


tab_alpha <- plotDensity(samples,filter="sample",selection="alpha",mapped=TRUE,layout=c(2,7))

p_test(samples,p_name="sd_FWwords:Sold",x_selection = "mu",mapped=FALSE)


get_design(samples)


#### Parameter recovery study ----
print(load("probitFW.RData")) 

# Create a single simulated data set.

# Can make up new data using the mean of the alphas (median or random also possible) 
new_dat <- post_predict(samples,use_par="mean",n_post=1)
# Or by sampling a parameter vector from the hyper
new_dat_hyper <- post_predict(samples,hyper=TRUE,n_post=1)

# An attribute keeps the pars used to create
round(attr(new_dat,"pars"),2)
round(attr(new_dat_hyper,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard")
# save(samplers,file="RecoveryProbitFixed.RData")

samplers <- make_samplers(new_dat_hyper,design,type="standard")
# save(samplers,file="RecoveryProbitRandom.RData")


print(load("RecoveryProbitFixed.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))

