rm(list=ls())
source("emc/emc.R")
source("models/RACE/RDM/rdmB.R")

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

# Again we construct a series of plausible models to test. We do so in three sets,
# the first fixing moment-to-moment variability (s=1) and with no start point 
# variability (A), the second estimating A, and the third estimating s for 
# the matching accumulator and fixing it for the mismatching accumulator (as
# was done for sv in the LBA). Combining both A and s lead to sampling problems
# for more complex models and was not pursued. In each case we fit 4 models
# allowing only B to vary with emphasis, B & t0, B & v, and B, t0 & v.

# Each set of 4 models fit by a single file, using the run_emc
# function. This function is given the name of the file containing the samplers
# object and any arguments to its constituent functions (auto_burn, auto_adapt
# and auto_sample). 

# For example the following fits the first model in the first set with 8 cores 
# per chain and so 24 cores in total given the default of 3 chains. 
#
# run_emc("rdmPNAS_B.RData",cores_per_chain=8)
#
# If any stage fails check_run stops, saving the progress so far. The only 
# argument unique to this function is nsample, which determines the number of 
# iterations performed if the sample stage is reached (by default 1000).
#
# After the convergence checking reported below extra samples were obtained for
# some models that converged slowly, so all had 1000 converged samples. This was
# done by adding another run_emc call with nsamples set appropriately (run_emc
# recognizes that burn and adapt stages are complete and so just adds on the
# extra samples at the end). 
#
# Finally for the selected model extra samples were added to obtain 5000 
# iterations, so posterior inference was reliable, and posterior predictives 
# were obtained. (NB: if the file is run several times some lines must be 
# commented out or unwanted sample stage iteration will be added)


# A=0, s=1, fit by run_rdm.R

# Rates differ by latent match, B differs by emphasis and latent response
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1,s~1),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
rdm_B <- make_samplers(dat,design_B,type="standard",rt_resolution=.02)
save(rdm_B,file="rdmPNAS_B.RData")

design_Bt0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~E,s~1),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
rdm_Bt0 <- make_samplers(dat,design_Bt0,type="standard",rt_resolution=.02)
save(rdm_Bt0,file="rdmPNAS_Bt0.RData")

design_Bv <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~1,s~1),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
rdm_Bv <- make_samplers(dat,design_Bv,type="standard",rt_resolution=.02)
save(rdm_Bv,file="rdmPNAS_Bv.RData")

design_Bvt0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~E,s~1),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
rdm_Bvt0<- make_samplers(dat,design_Bvt0,type="standard",rt_resolution=.02)
save(rdm_Bvt0,file="rdmPNAS_Bvt0.RData")


# Add variation in A, , fit by run_rdm_A.R
design_B_A <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1,s~1),
  constants=c(s=log(1)),
  model=rdmB)
# rdm_B_A <- make_samplers(dat,design_B_A,type="standard",rt_resolution=.02)
# save(rdm_B_A,file="rdmPNAS_B_A.RData")

design_Bt0_A <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~E),
  model=rdmB)
# rdm_Bt0_A <- make_samplers(dat,design_Bt0_A,type="standard",rt_resolution=.02)
# save(rdm_Bt0_A,file="rdmPNAS_Bt0_A.RData")

# B and v only
design_Bv_A <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~1),
  model=rdmB)
# rdm_Bv_A <- make_samplers(dat,design_Bv_A,type="standard",rt_resolution=.02)
# save(rdm_Bv_A,file="rdmPNAS_Bv_A.RData")

design_Bvt0_A <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~E),
  model=rdmB)
# rdm_Bvt0_A <- make_samplers(dat,design_Bvt0_A,type="standard",rt_resolution=.02)
# save(rdm_Bvt0_A,file="rdmPNAS_Bvt0_A.RData")


# Add variation in s, , fit by run_rdm_s.R

design_B_s <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~1,s~lM),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
# rdm_B_s <- make_samplers(dat,design_B_s,type="standard",rt_resolution=.02)
# save(rdm_B_s,file="rdmPNAS_B_s.RData")

design_Bt0_s <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~E,s~lM),
  constants=c(s=log(1)),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
# rdm_Bt0_s <- make_samplers(dat,design_Bt0_s,type="standard",rt_resolution=.02)
# save(rdm_Bt0_s,file="rdmPNAS_Bt0_s.RData")

design_Bv_s <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~1,s~lM),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
# rdm_Bv_s <- make_samplers(dat,design_Bv_s,type="standard",rt_resolution=.02)
# save(rdm_Bv_s,file="rdmPNAS_Bv_s.RData")

design_Bvt0_s <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM*E,B~lR*E,A~1,t0~E,s~lM),
  constants=c(s=log(1),A=log(0)),
  model=rdmB)
# rdm_Bvt0_s <- make_samplers(dat,design_Bvt0_s,type="standard",rt_resolution=.02)
# save(rdm_Bvt0_s,file="rdmPNAS_Bvt0_s.RData")

#### Check convergence ----

# Load samples
print(load("models/RACE/RDM/examples/samples/rdmPNAS_B.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bt0.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bv.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0.RData"))

print(load("models/RACE/RDM/examples/samples/rdmPNAS_B_A.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bt0_A.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bv_A.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0_A.RData"))

print(load("models/RACE/RDM/examples/samples/rdmPNAS_B_s.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bt0_s.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bv_s.RData"))
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0_s.RData"))



check_run(rdm_B)
check_run(rdm_Bt0,subfilter=500)
check_run(rdm_Bv,subfilter=500)
check_run(rdm_Bvt0,subfilter=1500)

check_run(rdm_B_A)
check_run(rdm_Bt0_A,subfilter=500)
check_run(rdm_Bv_A,subfilter=1500)
check_run(rdm_Bvt0_A,subfilter=1500)

check_run(rdm_B_s,subfilter=500)
check_run(rdm_Bt0_s,subfilter=1000)
check_run(rdm_Bv_s,subfilter=4000)
check_run(rdm_Bvt0_s,subfilter=4000)


#### Model selection ----

# Bvt0s wins with Bvs second 
compare_IC(list(
  B=rdm_B,Bt0=rdm_Bt0,Bv=rdm_Bv,Bvt0=rdm_Bvt0,
  BA=rdm_B_A,Bt0A=rdm_Bt0_A,BvA=rdm_Bv_A,Bvt0A=rdm_Bvt0_A,
  Bs=rdm_B_s,Bt0s=rdm_Bt0_s,Bvs=rdm_Bv_s,Bvt0s=rdm_Bvt0_s),
           subfilter=list(0,500,500,1500,0,500,1500,1500,500,1000,4000,4001:5000))
# At the subject level, Bvs is a close second or equal
ICs <- compare_ICs(list(
  B=rdm_B,Bt0=rdm_Bt0,Bv=rdm_Bv,Bvt0=rdm_Bvt0,
  BA=rdm_B_A,Bt0A=rdm_Bt0_A,BvA=rdm_Bv_A,Bvt0A=rdm_Bvt0_A,
  Bs=rdm_B_s,Bt0s=rdm_Bt0_s,Bvs=rdm_Bv_s,Bvt0s=rdm_Bvt0_s),
           subfilter=list(0,500,500,1500,0,500,1500,1500,500,1000,4000,4001:5000))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$DIC)]})))
table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$BPIC)]})))

# Comparing with the best DDM (16 parameter) model to the best (15 parameter)
# LBA model and best (16 parameter) RDM the latter comes third. 
source("models/DDM/DDM/ddmTZD.R")
print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full.RData")) 
source("models/RACE/LBA/lbaB.R")
print(load("models/RACE/LBA/examples/samples/sPNAS_Bv_sv.RData"))
compare_IC(list(DDM_avt0=sPNAS_avt0_full,LBA_Bvsv=sPNAS_Bv_sv,RDM_Bvt0s=rdm_Bvt0_s),
           subfilter=list(501:1500,2001:3000,4001:5000))

# For the remaining analysis add 4000 iterations to Bvt0 model

####  Fit of winning model ----

# post predict did not use first 4000, so inference based on 5000*3 samples
plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pprdm_Bvt0_s,layout=c(2,3),lpos="right")

# Good fit
pc <- function(d) 100*mean(d$S==d$R)
plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),
         stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
tab <- plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),
                stat=function(d){mean(d$rt)},stat_name="Mean RT (s)",xlim=c(0.375,.625))
tab <- plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
                stat=function(d){quantile(d$rt,.1)},stat_name="10th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
                stat=function(d){quantile(d$rt,.9)},stat_name="90th Percentile (s)")
tab <- plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),
                stat=function(d){sd(d$rt)},stat_name="SD (s)",xlim=c(0.1,.355))
tab <- plot_fit(dat,pprdm_Bvt0_s,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
                stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")

#### Posterior parameter inference ----

### Population mean (mu) tests

# Priors all well dominated except some rate parameters
cirdm_Bvt0_s <- plot_density(rdm_Bvt0_s,layout=c(2,8),selection="mu",subfilter=4000)
# On the natural scale it is evident this is because the mismatch (FALSE) rates
# are least well updated, due to the fairly low error rates (errors give the
# most information about FALSE rates).
cirdm_Bvt0_s_mapped <- plot_density(rdm_Bvt0_s,layout=c(3,6),
                                  selection="mu",mapped=TRUE)
# Looking at parameters both with and without mapping
round(cirdm_Bvt0_s,2)
round(cirdm_Bvt0_s_mapped,2)

# For simpler estimates:
# 1) t0 is longer than the LBA
# 2) Start point noise slightly larger relative to B than for the LBA.

### B effects

# Recall the map used with 
get_map(rdm_Bvt0_s)$B

# B_lRd tests threshold right - threshold left, as for the LBA, although not 
# quite credible it indicates slightly higher right thresholds (i.e., a bias to 
# respond left)
p_test(x=rdm_Bvt0_s,x_name="B_lRd",subfilter=4000,digits=3)

# B_Ea-n and B_Ea-s measure differences in response caution (i.e., thresholds
# averaged over left and right accumulators), accuracy-neutral and accuracy-speed 
# respectively. Caution for accuracy is clearly higher than speed, but not 
# credibly greater than neutral.
p_test(x=rdm_Bvt0_s,x_name="B_Ea-n",subfilter=4000)
p_test(x=rdm_Bvt0_s,x_name="B_Ea-s",subfilter=4000)

# Here we construct a test on the natural scale showing caution is greater for 
# neutral than speed
p_test(x=rdm_Bvt0_s,mapped=TRUE,x_name="average B: neutral-speed",
       x_fun=function(x){mean(x[c("B_left_neutral","B_right_neutral")]) - 
           mean(x[c("B_left_speed","B_right_speed")])})

# The remaining terms test interactions with bias (i.e., lR), with evidence of
# a small but credibly stronger bias to respond left (i.e., a lower threshold
# for the left accumulator) for speed than accuracy.
p_test(x=rdm_Bvt0_s,x_name="B_lRd:Ea-s",subfilter=4000,digits=2)

### v effects

# Again recall the map used with 
get_map(rdm_Bvt0_s)$v

# v_Ea-n v_Ea-s indicate that processing rate (the average of matching and 
# mismatching rates) is less in the accuracy condition than in neutral or speed.
p_test(x=rdm_Bvt0_s,x_name="v_Ea-n",subfilter=4000)
p_test(x=rdm_Bvt0_s,x_name="v_Ea-s",subfilter=4000)

# However, neutral and speed do not credibly differ
p_test(x=rdm_Bvt0_s,mapped=TRUE,x_name="average v: speed-neutral",
       x_fun=function(x){mean(x[c("v_TRUE_speed","v_FALSE_speed")]) - 
           mean(x[c("v_TRUE_neutral","v_FALSE_neutral")])})

# v_lMd tests the quality of selective attention, rate match - rate mismatch
p_test(x=rdm_Bvt0_s,x_name="v_lMd",subfilter=4000)

# v_Ea-n indicates quality does not differ credibly between accuracy and neutral.
p_test(x=rdm_Bvt0_s,x_name="v_lMd:Ea-n",subfilter=4000)

# In contrast there is a strong difference for accuracy - speed
p_test(x=rdm_Bvt0_s,x_name="v_lMd:Ea-s",subfilter=4000)

# The neutral-speed difference is also highly credible, here is is tested
# in the mapped form
p_test(x=rdm_Bvt0_s,mapped=TRUE,x_name="quality: accuracy-neutral",
       x_fun=function(x){diff(x[c("v_FALSE_neutral","v_TRUE_neutral")]) - 
           diff(x[c("v_FALSE_speed","v_TRUE_speed")])})

# Finally we see that just as in the LBA (for sv) rate variability is lower for
# match than mismatch 
p_test(x=rdm_Bvt0_s,mapped=TRUE,x_name="rate variability: mismatch - match",
       x_fun=function(x){diff(x[c("s_TRUE","s_FALSE")])})


# -------------------------------------------------------------------------

# Here we visualize these threshold and rate effects
library(tidyverse)

# Get parameter data frame
ps <- parameters_data_frame(rdm_Bvt0_s, mapped = TRUE)
head(ps)

# Get means
pmeans <- data.frame(lapply(ps, mean))
pmeans

# Here we look at thresholds by emphasis and accumulator
B <- pmeans %>% transmute(Accuracy_Left = B_left_accuracy,
                          Accuracy_Right = B_right_accuracy,
                          Neutral_Left = B_left_neutral,
                          Neutral_Right = B_right_neutral,
                          Speed_Left = B_left_speed,
                          Speed_Right = B_right_speed) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis", "Accumulator"),
               names_pattern = "(.*)_(.*)",
               values_to = "B")
B

# Plot thresholds
B_plot <- ggplot(data = B, aes(Emphasis, B, col = Accumulator)) +
  geom_line(aes(group = Accumulator), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  ylim(0.5, 1) + 
  theme_minimal() + theme(text = element_text(size = 15))
B_plot


# Here instead we will calculate the difference between thresholds for left and 
# right responses to better illustrate the change in bias across emphasis conditions
Bmeans <- data.frame(
  a_bias = mean(ps$B_right_accuracy - ps$B_left_accuracy),
  n_bias = mean(ps$B_right_neutral - ps$B_left_neutral),
  s_bias = mean(ps$B_right_speed - ps$B_left_speed)
)
Bmeans

# Long format
B_bias <- Bmeans %>% transmute(Accuracy = a_bias,
                               Neutral = n_bias,
                               Speed = s_bias) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis"),
               names_pattern = "(.*)",
               values_to = "Bias")
B_bias

# Plot threshold bias
B_plot <- ggplot(data = B_bias, aes(Emphasis, Bias)) +
  geom_line(aes(group = 1), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  ylim(0, 0.12) + 
  theme_minimal() + theme(text = element_text(size = 15))
B_plot


# Now we will calculate the quality (match-mismatch difference) and quantity
# (match-mismatch average) of accumulation rates
vmeans <- data.frame(
  a_quality = mean(ps$v_TRUE_accuracy - ps$v_FALSE_accuracy),
  n_quality = mean(ps$v_TRUE_neutral - ps$v_FALSE_neutral),
  s_quality = mean(ps$v_TRUE_speed - ps$v_FALSE_speed),
  a_quantity = mean(c(ps$v_TRUE_accuracy, ps$v_FALSE_accuracy)),
  n_quantity = mean(c(ps$v_TRUE_neutral, ps$v_FALSE_neutral)),
  s_quantity = mean(c(ps$v_TRUE_speed, ps$v_FALSE_speed))
)
vmeans

# Long format
QQ <- vmeans %>% transmute(Accuracy_Quality = a_quality,
                           Neutral_Quality = n_quality,
                           Speed_Quality = s_quality,
                           Accuracy_Quantity = a_quantity,
                           Neutral_Quantity = n_quantity,
                           Speed_Quantity = s_quantity) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis", "Measure"),
               names_pattern = "(.*)_(.*)",
               values_to = "v")
QQ

# Plot quality and quantity
# Similar to what we saw with the LBA, quality degrades under speed emphasis
# but participants deploy more overall resources (quantity) to meet the
# extra demands of faster responding
QQ_plot <- ggplot(data = QQ, aes(Emphasis, v, col = Measure)) +
  geom_line(aes(group = Measure), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  ylim(1.8, 3.6) +
  theme_minimal() + theme(text = element_text(size = 15))
QQ_plot

# -------------------------------------------------------------------------

