rm(list=ls())
source("emc/emc.R")
source("models/RACE/LNR/lnrMS.R")

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

#### Models -----

# Mu varies by stimulus, emphasis, and latent match
design_mu <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~1,t0~1),
  model=lnrMS)
# lnr_mu <- make_samplers(dat,design_mu,type="standard",rt_resolution=.02)
# save(lnr_mu,file="lnrPNAS_mu.RData")

# Let sigma also vary by latent match
design_mu_slM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~lM,t0~1),
  model=lnrMS)
# lnr_mu_slM <- make_samplers(dat,design_mu_slM,type="standard",rt_resolution=.02)
# save(lnr_mu_slM,file="lnrPNAS_mu_slM.RData")

# Let t0 vary by speed-accuracy emphasis
design_mu_slM_t0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~lM,t0~E),
  model=lnrMS)
# lnr_mu_slM_t0 <- make_samplers(dat,design_mu_slM_t0,type="standard",rt_resolution=.02)
# save(lnr_mu_slM_t0,file="lnrPNAS_mu_slM_t0.RData")

# Let sigma vary by emphasis and latent match
design_mu_sElM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~E*lM,t0~1),
  model=lnrMS)
# lnr_mu_sElM <- make_samplers(dat,design_mu_sElM,type="standard",rt_resolution=.02)
# save(lnr_mu_sElM,file="lnrPNAS_mu_sElM.RData")

# Again let t0 vary by emphasis
design_mu_sElM_t0 <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~E*lM,t0~E),
  model=lnrMS)
# lnr_mu_sElM_t0 <- make_samplers(dat,design_mu_sElM_t0,type="standard",rt_resolution=.02)
# save(lnr_mu_sElM_t0,file="lnrPNAS_mu_sElM_t0.RData")


# Load saved samples
print(load("samples/models/RACE/LNR/examples/lnrPNAS_mu.RData"))
print(load("samples/models/RACE/LNR/examples/lnrPNAS_mu_slM.RData"))
print(load("samples/models/RACE/LNR/examples/lnrPNAS_mu_slM_t0.RData"))
print(load("samples/models/RACE/LNR/examples/lnrPNAS_mu_sElM.RData"))
print(load("samples/models/RACE/LNR/examples/lnrPNAS_mu_sElM_t0.RData"))

# Sampling checks
check_run(lnr_mu,layout=c(3,5),subfilter=500)
check_run(lnr_mu_slM,layout=c(3,5),subfilter=500)
check_run(lnr_mu_slM_t0,layout=c(3,5),subfilter=1000)
check_run(lnr_mu_sElM,layout=c(3,5),subfilter=1500)
check_run(lnr_mu_sElM_t0,layout=c(3,5),subfilter=2000:3000)


############ Model selection

# muElMt0 wins by DIC but mulMt0 by BPIC
compare_IC(list(mu=lnr_mu,mulM=lnr_mu_slM,mulMt0=lnr_mu_slM_t0,
                muElM=lnr_mu_sElM,muElMt0=lnr_mu_sElM_t0), subfilter=list(500,500,1000,1500,2000))
# On an individual basis muElMt0 also wins
ICs <- compare_ICs(list(mu=lnr_mu,mulM=lnr_mu_slM,mulMt0=lnr_mu_slM_t0,
                        muElM=lnr_mu_sElM,muElMt0=lnr_mu_sElM_t0), subfilter=list(0,0,0,500,1000))
sort(table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$DIC)]}))))
sort(table(unlist(lapply(ICs,function(x){row.names(x)[which.min(x$BPIC)]}))))

# Comparing with the best DDM (16 parameter) model to the best (15 parameter)
# LBA model and best (16 parameter) RDM and best (21 parameter) LNR, muElMt0 
# comes in second. 
source("models/DDM/DDM/ddmTZD.R")
print(load("samples/models/DDM/DDM/examples/sPNAS_avt0_full.RData")) 
source("models/RACE/LBA/lbaB.R")
print(load("samples/models/RACE/LBA/examples/sPNAS_Bv_sv.RData"))
source("models/RACE/RDM/rdmB.R")
print(load("samples/models/RACE/RDM/examples/rdmPNAS_Bvt0_s.RData"))

compare_IC(list(DDM=sPNAS_avt0_full,LBA=sPNAS_Bv_sv,RDM=rdm_Bvt0_s,
                LNRmulMt0=lnr_mu_sElM_t0),subfilter=list(500,2000,1500,4000))

# Here we add in the version of the winning DDM model that did not use cell
# coding. We see that there is very little effect on the ICs
print(load("samples/models/DDM/DDM/examples/sPNAS_avt0_full_nocell.RData")) 

compare_IC(list(DDM=sPNAS_avt0_full_nocell,LBA=sPNAS_Bv_sv,RDM=rdm_Bvt0_s,
                LNR=lnr_mu_sElM_t0),subfilter=list(500,2000,1500,4000))

# For the remaining analysis add 4000 iterations to Bvt0 model

####  Fit of winning (E on mu, sigma and t0) model ----

# post predict did not use first 2000, so inference based on 5000*3 samples
plot_fit(dat,pplnr_musElMt0,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))
plot_fits(dat,pplnr_musElMt0,layout=c(2,3),lpos="right")

# Very accurate fit
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

# Priors all well dominated 
cilnr_mu_sElM_t0 <- plot_density(lnr_mu_sElM_t0,layout=c(2,3),selection="mu",subfilter=2000)
# On the natural scale it is evident this is because the mismatch (FALSE) rates
# are least well updated, due to the fairly low error rates (errors give the
# most information about FALSE rates).
cilnr_mu_sElM_t0_mapped <- plot_density(lnr_mu_sElM_t0,layout=c(2,3),subfilter=2000,
                                        selection="mu",mapped=TRUE)
# Looking at parameters both with and without mapping
round(cilnr_mu_sElM_t0,2)
round(cilnr_mu_sElM_t0_mapped,2)

### m(ean) effects

get_map(lnr_mu_sElM_t0)$m

# lMd tests match minus mismatch, which is strongly negative, as the match
# racer is much faster than the mismatch racer
p_test(x=lnr_mu_sElM_t0,x_name="m_lMd",subfilter=2000)

# m_Sd tests main effect of stimulus (the mean for (right - left), and indicates 
# a smaller mean (faster processing) for right (-.92) than left (-0.88).
p_test(x=lnr_mu_sElM_t0,x_name="m_Sd",subfilter=2000)

# m_Sd:lMd tests right(match-mismatch) - left(match-mismatch) which is positive,
# indicating better discrimination for right than left. Rearranging this is also
# (right_match+left_mismatch) - (left_match+right_mismatch) =
# right accumulator - left accumulator, so it is also a measure of response bias
# the positive value indicates a bias to respond left (i.e., a lesser time for 
# the left racer to finish). These two interpretations cannot be disambiguated 
# from LNR parameter estimates.
p_test(x=lnr_mu_sElM_t0,x_name="m_Sd:lMd",subfilter=2000)

# m_Ea-n and m_Ea-s measure differences the time overall race time between
# accuracy and neutral/speed. Both are positive indicating the race time is
# slowest for accuracy
p_test(x=lnr_mu_sElM_t0,x_name="m_Ea-n",subfilter=2000)
p_test(x=lnr_mu_sElM_t0,x_name="m_Ea-s",subfilter=2000)
# We can test neutral - speed as follows, showing speed is quicker than neutral
p_test(x=lnr_mu_sElM_t0,x_name="m_Ea-s",y=lnr_mu_sElM_t0,y_name="m_Ea-n",subfilter=2000)

# We will not look at the remaining interactions here

### s(igma) effects

# Again recall the map used with 
get_map(lnr_mu_sElM_t0)$s

# This gives the overall mean on the natural scale
p_test(x=lnr_mu_sElM_t0,x_name="sigma",subfilter=2000,x_fun=function(x){exp(x["s"])})

# s_Ea-n s_Ea-s show sigma is credibly greater for accuracy than neutral, but
# not for speed
p_test(x=lnr_mu_sElM_t0,x_name="s_Ea-n",subfilter=2000)
p_test(x=lnr_mu_sElM_t0,x_name="s_Ea-s",subfilter=2000)

# Sigma is substantially less for match than mismatch
p_test(x=lnr_mu_sElM_t0,x_name="s_lMd",subfilter=2000)

# s_Ea-n:lMd s_Ea-s:lMd show no credible differences in the match vs. mismatch
# difference between accuracy and neutral or speed, suggesting an additive
# model would be better (E+lM)
p_test(x=lnr_mu_sElM_t0,x_name="s_Ea-n:lMd",subfilter=2000)
p_test(x=lnr_mu_sElM_t0,x_name="s_Ea-s:lMd",subfilter=2000)


### t0 effects

# Recall the map used with 
get_map(lnr_mu_sElM_t0)$t0

# From the tables we see that t0 for accuracy (0.26s) is credibly greater
# than for neutral and speed. The following test shows that neutral is 
# credibly greater than speed.
p_test(mapped=TRUE,subfilter=2000,x=lnr_mu_sElM_t0,x_name="t0_neutral",
       y=lnr_mu_sElM_t0,y_name="t0_speed")
p_test(mapped=TRUE,subfilter=2000,x=lnr_mu_sElM_t0,x_name="t0_accuracy",
       y=lnr_mu_sElM_t0,y_name="t0_neutral")


# -------------------------------------------------------------------------

# Make parameter plots
library(tidyverse)

# Get parameter data frame
ps <- parameters_data_frame(lnr_mu_sElM_t0, mapped = TRUE)
head(ps)

# Calculate means
pmeans <- data.frame(lapply(ps, mean))
pmeans

# Plot means
m <- pmeans %>% transmute(Accuracy_left_TRUE = m_left_accuracy_TRUE,
                          Accuracy_left_FALSE = m_left_accuracy_FALSE,
                          Accuracy_right_TRUE = m_right_accuracy_TRUE,
                          Accuracy_right_FALSE = m_right_accuracy_FALSE,
                          Neutral_left_TRUE = m_left_neutral_TRUE,
                          Neutral_left_FALSE = m_left_neutral_FALSE,
                          Neutral_right_TRUE = m_right_neutral_TRUE,
                          Neutral_right_FALSE = m_right_neutral_FALSE,
                          Speed_left_TRUE = m_left_speed_TRUE,
                          Speed_left_FALSE = m_left_speed_FALSE,
                          Speed_right_TRUE = m_right_speed_TRUE,
                          Speed_right_FALSE = m_right_speed_FALSE) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis", "Stimulus", "Match"),
               names_pattern = "(.*)_(.*)_(.*)",
               values_to = "m")
m

m_plot <- ggplot(data = m, aes(Emphasis, m, col = Match)) +
  geom_line(aes(group = Match), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  facet_grid(Stimulus ~ .) +
  ylim(-2, 1) +
  theme_minimal() + theme(text = element_text(size = 15))
m_plot


# Plot sigmas
s <- pmeans %>% transmute(Accuracy_TRUE = s_accuracy_TRUE,
                          Accuracy_FALSE = s_accuracy_FALSE,
                          Neutral_TRUE = s_neutral_TRUE,
                          Neutral_FALSE = s_neutral_FALSE,
                          Speed_TRUE = s_speed_TRUE,
                          Speed_FALSE = s_speed_FALSE) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis", "Match"),
               names_pattern = "(.*)_(.*)",
               values_to = "s")
s

s_plot <- ggplot(data = s, aes(Emphasis, s, col = Match)) +
  geom_line(aes(group = Match), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  ylim(0.3, 1) + 
  theme_minimal() + theme(text = element_text(size = 15))
s_plot


# Plot non-decision time
t0 <- pmeans %>% transmute(Accuracy = t0_accuracy,
                           Neutral = t0_neutral,
                           Speed = t0_speed) %>%
  pivot_longer(cols = everything(),
               names_to = c("Emphasis"),
               names_pattern = "(.*)",
               values_to = "t0")
t0

t0_plot <- ggplot(data = t0, aes(Emphasis, t0)) +
  geom_line(aes(group = 1), lty = "dashed", lwd = 1) +
  geom_point(size = 4, shape = 19) +
  ylim(0.18, 0.3) + 
  theme_minimal() + theme(text = element_text(size = 15))
t0_plot

# -------------------------------------------------------------------------

