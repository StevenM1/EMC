rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/exGaussian/exGaussian.R")


# Two choice normal
designEXG <- make_design(Flist=list(mu ~ lM, sigma ~ 1, tau ~ 1),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
  model=exGaussian)

p_vector <- sampled_p_vector(designEXG)
p_vector[1:2] <- log(c(.5,.5))
p_vector[3:4] <- log(c(.2,.4)) # sd

dataEXG <- make_data(p_vector,design=designEXG,trials=1000)
plot_defective_density(dataEXG,layout=c(1,2))

dadmEXG <- design_model(dataEXG,designEXG)
par(mfrow=c(1,4))
profile_pmwg(pname="mu",p=p_vector,p_min=log(.25),p_max=log(.75),dadm=dadmEXG)
profile_pmwg(pname="mu_lM1",p=p_vector,p_min=log(.25),p_max=log(.75),dadm=dadmEXG)
profile_pmwg(pname="sigma",p=p_vector,p_min=log(.1),p_max=log(.3),dadm=dadmEXG)
profile_pmwg(pname="tau",p=p_vector,p_min=log(.3),p_max=log(.5),dadm=dadmEXG)


# This example shows EMC works in a model without choice (Rlevels=1, no matchfun)
# Make a reasonably complex design. 
designEXG <- make_design(Flist=list(mu ~ A*B, sigma ~ B, tau ~ 1),
  Ffactors=list(subjects=50,A=1:3,B=1:2),Rlevels=1,model=exGaussian)

p_vector <- sampled_p_vector(designEXG)
p_vector[1:6] <- log(c(.5,1.2,1.4,1.5,1,1)) # Increasing in A and B, no interaction
p_vector[7:9] <- log(c(.1,2,.3)) # more variable for B2 than B1 


dataEXG <- make_data(p_vector,design=designEXG,trials=1000)
plot_defective_density(dataEXG,layout=c(2,3),rt="topleft")

dadmEXG <- design_model(dataEXG,designEXG)
par(mfrow=c(2,5))
for (i in names(p_vector))
  profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmEXG)

samplers <- make_samplers(dataEXG,designEXG,type="standard",rt_resolution=.001)

