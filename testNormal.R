rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/RACE/Normal/normal.R")


# Two choice normal
designNORMAL <- make_design(Flist=list(mean ~ lM, sd ~ 1),
                            Ffactors=list(subjects=1:10,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
                            Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
                            constants = c(sd = 0),
                            model=normal)

p_vector <- sampled_p_vector(designNORMAL)
p_vector[1:2] <- c(0,1)

Sigma <- matrix(c(.08, .02,
                  .02, .1), nrow = 2)

p_mat <- rmvnorm(length(designNORMAL$Ffactors$subjects),mean=p_vector,sigma=Sigma)
dimnames(p_mat)[[1]] <- designNORMAL$Ffactors$subjects
dataNORMAL <- make_data(p_vector= p_mat,design=designNORMAL,trials=200)
plot_defective_density(dataNORMAL,layout=c(1,2))

dadmNORMAL <- design_model(dataNORMAL,designNORMAL)
par(mfrow=c(1,3))
profile_pmwg(pname="mean",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmNORMAL)
profile_pmwg(pname="mean_lM1",p=p_vector,p_min=0,p_max=2,dadm=dadmNORMAL)
# profile_pmwg(pname="sd",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmNORMAL)

sampler <- make_samplers(dadmNORMAL, designNORMAL, type = "single")
#samples <- run_chains(sampler, iter = c(50, 0, 0), cores_per_chain = 6, cores_for_chains = 1, verbose_run_stage = T)
source("emc/emc.R")

# debug(auto_burn)
burned <- auto_burn(sampler,nstart = 100, ndiscard = 20, cores_per_chain = 10, cores_for_chains = 1, verbose_run_stage = T, step_size = 50)
adapted <- auto_adapt(burned, cores_for_chains = 1, cores_per_chain = 10, verbose_run_stage = T)
sampled <- auto_sample(adapted, iter = 100, cores_for_chains = 1, cores_per_chain = 10, verbose_run_stage = T)
