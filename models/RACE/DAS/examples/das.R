rm(list=ls())
source("emc/emc.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
# head(dat)

# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

#### Models ----
source("models/RACE/DAS/dasWald.R")
design_Wald <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,t0~1,s~1,tau~1),
  constants=c(s=log(1)),
  model=dasWald)

# das_Wald <- make_samplers(dat,design_B,type="standard",rt_resolution=.02)
# save(rdm_B,file="rexdmPNAS_B.RData")

# Use Wiener model fit to get plausible parameters
print(load("models/RACE/RDM/examples/samples/rdmPNAS_B.RData"))
p_vector <- apply(parameters_data_frame(rdm_B),2,mean)
names(p_vector)[9] <- "tau"; p_vector["tau"] <- log(.5)

# All identical subjects
simdata <- make_data(p_vector,design=design_Wald,trials=150)
plot_defective_density(simdata,layout=c(2,3),factors=c("E","S"))


dadm <- design_model(simdata,design_Wald)
# log_likelihood_DAS(p_vector,dadm=dadm,pnames=c("B","v"))
par(mfrow=c(2,5))
profile_pmwg(pname="tau",p=p_vector,p_min=-1,p_max=-.25,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="v",p=p_vector,p_min=.4,p_max=1,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="v_lMd",p=p_vector,p_min=1,p_max=2,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B_lRd",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B_Ea-n",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B_Ea-s",p=p_vector,p_min=0,p_max=1,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B_lRd:Ea-n",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="B_lRd:Ea-s",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
profile_pmwg(pname="t0",p=p_vector,p_min=-2,p_max=-1,dadm=dadm,n_point=96,n_cores=32)

# # Cant get this to work
# source("models/RACE/DAS/dasNormal.R")
# design_Normal <- make_design(
#   Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
#   Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
#   Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
#   Flist=list(mean~lM,sd~1,t0~1,tau~1),
#   constants=c(t0=log(0)),
#   model=dasNormal)
# 
# 
# # All identical subjects
# p_vector <- log(c(mean=.7,mean_lMd=2,sd=.5,tau=1)) 
# simdata <- make_data(p_vector,design=design_Normal,trials=150)
# plot_defective_density(simdata,layout=c(2,3),factors=c("E","S"))
# 
# 
# dadm <- design_model(simdata,design_Normal)
# # log_likelihood_DAS(p_vector,dadm=dadm,pnames=c("mean","sd"))
# par(mfrow=c(2,3))
# profile_pmwg(pname="mean",p=p_vector,p_min=-.75,p_max=0,dadm=dadm,n_point=96,n_cores=32)
# profile_pmwg(pname="mean_lMd",p=p_vector,p_min=0,p_max=1,dadm=dadm,n_point=96,n_cores=32)
# profile_pmwg(pname="sd",p=p_vector,p_min=-1,p_max=0,dadm=dadm,n_point=96,n_cores=32)
# profile_pmwg(pname="tau",p=p_vector,p_min=-.5,p_max=.5,dadm=dadm,n_point=96,n_cores=32)
