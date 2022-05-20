rm(list=ls())
library(mvtnorm)
library(coda)
source("emc/emc.R")
source("models/LNR/lnrMS.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

# Accuracy - Speed, Neutral - Speed
Emat <- contr.sum(3)/2
dimnames(Emat) <- list(NULL,c("as","ns"))

design <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(m~E*lM,s~lM,t0~1),
  model=lnrMS)

# Design has a p_vector attribute that can be filled in with values to simulate
# data, and which itself has a "map" attribute that defines the nature of each
# parameter
attr(design,"p_vector")
#         m     m_Eas     m_Ens     m_lMd m_Eas:lMd m_Ens:lMd         s     s_lMd        t0 
#         0         0         0         0         0         0         0         0         0 
# attr(,"map")
# attr(,"map")$m
#   m m_Eas m_Ens m_lMd m_Eas:lMd m_Ens:lMd
# 1 1   0.5   0.0   0.5      0.25      0.00
# 3 1   0.5   0.0  -0.5     -0.25      0.00
# 2 1   0.0   0.5  -0.5      0.00     -0.25
# 4 1   0.0   0.5   0.5      0.00      0.25
# 
# attr(,"map")$s
#   s s_lMd
# 1 1   0.5
# 3 1  -0.5
# 
# attr(,"map")$t0
# [1] 1

# Full (variance-covariance) fits ----

# Combine data and design, define type of PMWG sampler and rt_resolution
# By default three chains are created (n_chains). Can also define priors_list
# (if none supplied defaults are used) and extras required by certian types
# (par_groups, n_factors, contraintMat, covariates)
samplers <- make_samplers(dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testLNR.RData")

print(load("testLNR.RData"))


# How many iterations
chain_n(samples)
#      burn adapt sample
# [1,]  406    75    500
# [2,]  406    82    500
# [3,]  406    75    500
     
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))


pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecoveryLNR.RData")

# run in script runRecovery/runRecovery100 
print(load("testRecoveryLNR.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))


##################### MISSING

source("LNR/MlnrMS.R")

# #### severe censoring above and truncation below
# 
# datM <- dat
# 
# # Larger than typical, 10% above and 10% below
# lower <- tapply(datM$rt,datM$subjects,quantile,probs=.1)
# upper <- tapply(datM$rt,datM$subjects,quantile,probs=.9)
# 
# # Truncate below
# fast <- as.numeric(as.character(factor(datM$subjects,labels=lower)))
# datM <- datM[datM$rt>fast,]
# attr(datM,"LT") <- lower
# # Size and direction of censor known but not response
# slow <- as.numeric(as.character(factor(datM$subjects,labels=upper)))
# datM$rt[datM$rt>=slow] <-  Inf
# datM$R[datM$rt>=slow] <-  NA
# attr(datM,"UC") <- upper
# 
# par(mfcol=c(2,3))
# sacc <- plot_defective_density(datM,correct_fun=function(data) data$S == data$R,xlim=c(.2,1))
# round(sort(sacc),2)
# 
# 
# designM <- make_design(
#   Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
#   Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
#   Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
#   Flist=list(m~E*lM,s~lM,t0~1),
#   model=MlnrMS)
# 
# 
# # dadmM <- design_model(data=datM,design=designM,rt_resolution=.02)
# # p_vector <- attr(designM,"p_vector")
# # p_vector[1:9] <- c(0,0,0,-1,-.5,-.5,log(1),log(1),log(0.3))
# # # Intercept and quality of meanlog effectively inverse
# # # of rate so match is LESS than mismatch
# # 
# # print(attr(dadmM,"model")$log_likelihood(p_vector,dadmM))
# # 
# # # individual subject likelihood
# # dadml <- dm_list(dadmM)
# # for (i in levels(dat$subjects))
# #  print(attr(dadmM,"model")$log_likelihood(p_vector=p_vector,dadm=dadml[[i]] ))
# 
# 
# samplers <- make_samplers(data_list=datM,design_list=designM,type="standard",rt_resolution=.02)
# # save(samplers,file="testLNRM.RData")
# 
# print(load("testLNRM.RData"))

#### Actual truncation on data set

datT <- dat

attr(datT,"LT") <- 0.25
attr(datT,"UT") <- 1.5

samplers <- make_samplers(data_list=datT,design_list=designM,type="standard",rt_resolution=.02)
save(samplers,file="testLNRT.RData")

print(load("testLNR.RData")); samplesNT <- samples
#     user   system  elapsed 
# 7834.661 2661.759 2037.509 
print(load("testLNRT.RData")); samplesT <- samples
#      user    system   elapsed 
# 91055.778  6600.869 15646.189

# How many iterations
chain_n(samples)
#      burn adapt sample
# [1,]  975     0      0
# [2,]  975     0      0
# [3,]  975     0      0
     
     
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))

round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))

pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecoveryLNR.RData")

# run in script runRecovery/runRecovery100 
print(load("testRecoveryLNR.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))


#### Actual truncation on data set + 3SD censoring 

datTC <- dat
UC <- tapply(datTC$rt,datTC$subjects,function(x){mean(x)+3*sd(x)})
slow <- datTC$subjects
levels(slow) <- UC
slow <- as.numeric(as.character(slow))
isSlow <- datTC$rt > slow
round(100*(1-tapply(isSlow,datTC$subjects,mean)),1)
# as1t bd6t bl1t hsft hsgt kd6t kd9t kh6t kmat ku4t na1t rmbt rt2t rt3t rt5t scat ta5t vf1t zk1t 
#  0.6  1.3  1.5  1.3  2.4  1.2  2.1  2.2  2.5  2.0  2.5  1.9  1.5  1.7  2.3  1.8  2.2  1.9  1.7 
 
attr(datTC,"LT") <- 0.25
attr(datTC,"UT") <- 1.5

datTC$rt[isSlow] <- Inf
attr(datTC,"UC") <- UC

# Response known
samplers <- make_samplers(data_list=datTC,design_list=designM,type="standard",rt_resolution=.02)
save(samplers,file="testLNRTC.RData")

# Response unknown
datTCR <- datTC
datTCR$R[isSlow] <- NA
samplers <- make_samplers(data_list=datTCR,design_list=designM,type="standard",rt_resolution=.02)
save(samplers,file="testLNRTCR.RData")



##################### FACTOR

samplers <- make_samplers(dat,design,type="factor",rt_resolution=.02,n_factors=1)
# save(samplers,file="testLNRfactor.RData")

print(load("testLNRfactor.RData"))

# Run progressively due to crash with sample stage, use adapt run
samples=samples1

# How many iterations
chain_n(samples)
#     burn adapt sample
# [1,]  600   300      0
# [2,]  600   304      0
# [3,]  600   302      0
    
    
filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))


pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="factor",rt_resolution=.02,
                          n_factors=1)
# save(samplers,file="testRecoveryLNRfactor.RData")

print(load("testRecoveryLNRfactor.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))



##################### diagonal

samplers <- make_samplers(dat,design,type="diagonal",rt_resolution=.02)
# save(samplers,file="testLNRdiagonal.RData")

print(load("testLNRdiagonal.RData"))
samples=samples1

# How many iterations
chain_n(samples)

filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))


pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="diagonal",rt_resolution=.02)
# save(samplers,file="testRecoveryLNRdiagonal.RData")

print(load("testRecoveryLNRdiagonal.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))


##################### BLOCKED

samplers <- make_samplers(dat,design,type="blocked",rt_resolution=.02,
                          par_groups=c(rep(1,6), rep(2, 2),3))
# save(samplers,file="testLNRblocked.RData")
# save(samplers,file="testLNRblocked1.RData")

print(load("testLNRblocked.RData"))
# print(load("testLNRblocked1.RData"))
# samples <- chain_shorten(samples,1120)



# How many iterations
chain_n(samples)

filter <- c("burn","adapt","sample")[1]
# filter <- c("burn","adapt","sample")[3]

gd <- gd_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter))
print(round(gd,2))

plotChains(as_mcmc.list(samples,selection="LL",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="mu",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="variance",filter=filter),layout=c(2,5))
plotChains(as_mcmc.list(samples,selection="covariance",filter=filter),layout=c(3,5))
plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),layout=c(2,5))

plotChains(as_mcmc.list(samples,selection="alpha",filter=filter),
           plot_acf=TRUE,acf_chain=1,layout=c(2,5))
round(es_pmwg(as_mcmc.list(samples,selection="alpha",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="mu",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="variance",filter=filter)))
round(es_pmwg(as_mcmc.list(samples,selection="covariance",filter=filter)))


pp <- post_predict(samples,filter="burn")
plot_fit(dat,pp,factors=c("E","S"),layout=c(2,3))
plot_fit(dat,pp,layout=c(2,3))

# Do recovery for data stimulated from posterior means
new_dat <- post_predict(samples,n_post=1,use_par="mean")
round(attr(new_dat,"pars"),2)

# can we recover these?
samplers <- make_samplers(new_dat,design,type="standard",rt_resolution=.02)
# save(samplers,file="testRecoveryLNRfactor.RData")

print(load("testRecoveryLNRblocked.RData"))

# After some checking of samples as above all looks good without trimming
pars <- attributes(attr(samples,"data_list")[[1]])$pars
tabs <- plotDensity(as_mcmc.list(samples,selection="alpha",filter="burn"),
                    layout=c(2,5),pars=pars)
# Some shrinkage but not bad
plotAlphaRecovery(tabs,layout=c(2,5))
