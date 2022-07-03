rm(list=ls())
source("emc/emc.R")
source("models/RACE/LNR/lnrMS.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))


design_mu <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat),
  Flist=list(m~E*lM,s~1,t0~1),
  model=lnrMS)
lnr_mu <- make_samplers(dat,design_mu,type="standard",rt_resolution=.02)
save(lnr_mu,file="lnrPNAS_mu.RData")


design_mu_slM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat),
  Flist=list(m~E*lM,s~lM,t0~1),
  model=lnrMS)
lnr_mu_slM <- make_samplers(dat,design_mu_slM,type="standard",rt_resolution=.02)
save(lnr_mu_slM,file="lnrPNAS_mu_slM.RData")

design_mu_slM_t0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat),
  Flist=list(m~E*lM,s~lM,t0~E),
  model=lnrMS)
lnr_mu_slM_t0 <- make_samplers(dat,design_mu_slM_t0,type="standard",rt_resolution=.02)
save(lnr_mu_slM_t0,file="lnrPNAS_mu_slM_t0.RData")

design_mu_sElM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat),
  Flist=list(m~E*lM,s~E*lM,t0~1),
  model=lnrMS)
lnr_mu_sElM <- make_samplers(dat,design_mu_sElM,type="standard",rt_resolution=.02)
save(lnr_mu_sElM,file="lnrPNAS_mu_sElM.RData")


design_mu_sElM_t0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat),
  Flist=list(m~E*lM,s~E*lM,t0~E),
  model=lnrMS)
lnr_mu_sElM_t0 <- make_samplers(dat,design_mu_sElM_t0,type="standard",rt_resolution=.02)
save(lnr_mu_sElM_t0,file="lnrPNAS_mu_sElM_t0.RData")



print(load("lnrPNAS_mu.RData"))
print(load("lnrPNAS_mu_slM.RData"))
print(load("lnrPNAS_mu_slM_t0.RData"))
print(load("lnrPNAS_mu_sElM.RData"))
print(load("lnrPNAS_mu_sElM_t0.RData"))

check_run(lnr_mu)
check_run(lnr_mu_slM)
check_run(lnr_mu_slM_t0)
check_run(lnr_mu_sElM,layout=c(3,5),subfilter=500)
check_run(lnr_mu_sElM_t0,layout=c(3,5),subfilter=1000)

############ Model selection

# muElMt0 wins by DIC but mulMt0 by BPIC
compare_IC(list(mu=lnr_mu,mulM=lnr_mu_slM,mulMt0=lnr_mu_slM_t0,
  muElM=lnr_mu_sElM,muElMt0=lnr_mu_sElM_t0), subfilter=list(0,0,0,500,1000))
# On an individual basis muElMt0 wins or ties for first.
ICs <- compare_ICs(list(mu=lnr_mu,mulM=lnr_mu_slM,mulMt0=lnr_mu_slM_t0,
  muElM=lnr_mu_sElM,muElMt0=lnr_mu_sElM_t0), subfilter=list(0,0,0,500,1000))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$DIC)]})))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$BPIC)]})))

# Comparing with the best DDM (16 parameter) model to the best (15 parameter)
# LBA model and best (16 parameter) RDM with muElMt0 (15 parameters) the LNR
# comes in third. 
source("models/DDM/DDM/ddmTZD.R")
print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full.RData")) 
source("models/RACE/LBA/lbaB.R")
print(load("models/RACE/LBA/examples/samples/sPNAS_Bv_sv.RData"))
source("models/RACE/RDM/rdmB.R")
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0.RData"))

compare_IC(list(DDM_avt0=sPNAS_avt0_full,LBA_Bvsv=sPNAS_Bv_sv,RDM_Bvt0=rdm_Bvt0,
                LNRmulMt0=lnr_mu_sElM_t0),subfilter=list(500,2000,1500,1000))

# For the remaining analysis add 4000 iterations to Bvt0 model

####  Fit of winning (Bvsv) model ----

# post predict did not use first 1500, so inference based on 5000*3 samples
plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pprdm_Bvt0,layout=c(2,3),lpos="right")

# Good fit, slight under-estimation of 10th percentile and over-estimation of 
# error RT in speed
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
tab <- plot_fit(dat,pprdm_Bvt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")



#####################################################

#### Model selection with Bayes factors

print(load("is2.RData"))
do.call(c,is2_sPNAS_Bv_sv)
do.call(c,is2_sPNAS_avt0_full)
do.call(c,is2_sPNAS_avt0_full_nocell)

std_error_IS2 <- std_error_IS2(IS_samples)
median(IS_samples)

source("models/DDM/DDM/ddmTZD.R")
design_avt0_full_nocell <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Flist=list(v~S*E,a~E,sv~1, t0~E, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1),DP=qnorm(0.5)),
  model=ddmTZD)

samplers <- make_samplers(dat,design_avt0_full_nocell,type="standard")
save(samplers,file="sPNAS_avt0_full_nocell.RData")

print(load("sPNAS_avt0_full_nocell.RData"))
# Nicely converged after 500
check_run(sPNAS_avt0_full_nocell,subfilter=500)

# Little difference between two parameterizations and full length Bvsv model 
compare_IC(list(avt0=sPNAS_avt0_full,avt0_nocell=sPNAS_avt0_full_nocell,Bvsv=sPNAS_Bv_sv),
           subfilter=list(500,500,2000))


