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

####  Fit of winning (E on mu, sigma and t0) model ----

# post predict did not use first 1500, so inference based on 5000*3 samples
plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pplnr_musElMt0,layout=c(2,3),lpos="right")

# Very good account of RT but underestimates accuracy for speed left and 
# over-estimates for speed right.
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
tab <- plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
         stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
         stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),
         stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
tab <- plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
         stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")


#### Posterior parameter inference ----

### Population mean (mu) tests

# Priors all well dominated except some rate parameters
cirdm_Bvt0 <- plot_density(lnr_mu_sElM_t0,layout=c(3,6),selection="mu",subfilter=1500)
# On the natural scale it is evident this is because the mismatch (FALSE) rates
# are least well updated, due to the fairly low error rates (errors give the
# most information about FALSE rates).
cilnr_mu_sElM_t0_mapped <- plot_density(lnr_mu_sElM_t0,layout=c(3,6),
                                        selection="mu",mapped=TRUE)
# Looking at parameters both with and without mapping
round(cilnr_mu_sElM_t0,2)
round(cilnr_mu_sElM_t0_mapped,2)

# For simpler estimates:
# 1) t0 is longer than the LBA
# 2) Start point noise slightly larger relative to B than for the LBA.

### mu effects

# Recall the map used with 
get_map(lnr_mu_sElM_t0)$B

# B_lRd tests threshold right - threshold left, as for the LBA, although not 
# quite credible it indicates slightly higher right thresholds (i.e., a bias to 
# respond left)
p_test(x=lnr_mu_sElM_t0,x_name="B_lRd",subfilter=1500,digits=3)

# B_Ea-n and B_Ea-s measure differences in response caution (i.e., thresholds
# averaged over left and right accumulators), accuracy-neutral and accuracy-speed 
# respectively. Caution for accuracy is clearly higher than speed, but not 
# credibly greater than neutral.
p_test(x=lnr_mu_sElM_t0,x_name="B_Ea-n",subfilter=1500)
p_test(x=lnr_mu_sElM_t0,x_name="B_Ea-s",subfilter=1500)

# Here we construct a test on the natural scale showing caution is greater for 
# neutral than speed
p_test(x=lnr_mu_sElM_t0,mapped=TRUE,x_name="average B: neutral-speed",
  x_fun=function(x){mean(x[c("B_left_neutral","B_right_neutral")]) - 
                    mean(x[c("B_left_speed","B_right_speed")])})

# The remaining terms test interactions with bias (i.e., lR), with evidence of
# a small but credibly stronger bias to respond left (i.e., a lower threshold
# for the left accumulator) for speed than accuracy.
p_test(x=lnr_mu_sElM_t0,x_name="B_lRd:Ea-s",subfilter=1500,digits=2)

### sigma effects

# Again recall the map used with 
get_map(lnr_mu_sElM_t0)$v

# v_Ea-n v_Ea-s indicate that processing rate (the average of matching and 
# mismatching rates) is less in the accuracy condition than in neutral or speed.
p_test(x=lnr_mu_sElM_t0,x_name="v_Ea-n",subfilter=1500)
p_test(x=lnr_mu_sElM_t0,x_name="v_Ea-s",subfilter=1500)

# However, neutral and speed do not credibly differ
p_test(x=lnr_mu_sElM_t0,mapped=TRUE,x_name="average v: speed-neutral",
  x_fun=function(x){mean(x[c("v_TRUE_speed","v_FALSE_speed")]) - 
                    mean(x[c("v_TRUE_neutral","v_FALSE_neutral")])})

# v_lMd tests the quality of selective attention, rate match - rate mismatch
p_test(x=lnr_mu_sElM_t0,x_name="v_lMd",subfilter=1500)

# v_Ea-n indicates quality does not differ credibly between accuracy and neutral.
p_test(x=lnr_mu_sElM_t0,x_name="v_lMd:Ea-n",subfilter=1500)

# In contrast there is a strong difference for accuracy - speed
p_test(x=lnr_mu_sElM_t0,x_name="v_lMd:Ea-s",subfilter=1500)

# The neutral - speed difference is also highly credible, here is is tested
# in the mapped form
p_test(x=lnr_mu_sElM_t0,mapped=TRUE,x_name="quality: accuracy-neutral",
  x_fun=function(x){diff(x[c("v_FALSE_neutral","v_TRUE_neutral")]) - 
                    diff(x[c("v_FALSE_speed","v_TRUE_speed")])})

### t0 effects

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


