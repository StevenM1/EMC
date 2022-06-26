rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/RACE/Normal/normal.R")


# Two choice normal
design1 <- make_design(Flist=list(mean ~ lM, sd ~ 1),
                            Ffactors=list(subjects=1:10,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
                            Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
                            constants = c(sd = 0),
                            model=normal)

p_vector <- sampled_p_vector(design1)
p_vector[1:2] <- c(0,1)

Sigma <- matrix(c(.08, .02,
                  .02, .1), nrow = 2)

p_mat1 <- rmvnorm(length(design1$Ffactors$subjects),mean=p_vector,sigma=Sigma)
dimnames(p_mat1)[[1]] <- design1$Ffactors$subjects
data1 <- make_data(p_vector= p_mat1,design=design1,trials=200)

dadm1 <- design_model(data1,design1)

# Two choice normal
design2 <- make_design(Flist=list(mean ~ lM*S, sd ~ 1),
                            Ffactors=list(subjects=1:10,S=1:2),Rlevels=1:2,matchfun=function(d)d$S==d$lR,
                            Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
                            constants = c(sd = 0),
                            model=normal)

p_vector <- sampled_p_vector(design2)
p_vector[1:4] <- c(0,1,2,3)

Sigma <- matrix(c(.12, .02, .03, .05,
                  .02, .1, .04, .05,
                  .03, .04, .2, .01,
                  .05, .05, .01, .25), nrow = 4)

p_mat2 <- rmvnorm(length(design2$Ffactors$subjects),mean=p_vector,sigma=Sigma)
dimnames(p_mat2)[[1]] <- design2$Ffactors$subjects
data2 <- make_data(p_vector= p_mat2,design=design2,trials=200)


dadm2 <- design_model(data2,design2)
sampler <- make_samplers(list(dadm1, dadm2), list(design1, design2), type = "standard")

samples <- run_chains(sampler, iter = c(25, 0, 0), cores_per_chain = 6, cores_for_chains = 1, verbose_run_stage = T)
plot_chains(samples, filter = "burn", selection = "mu", plot_acf = T)

gd_pmwg(samples, filter = "burn")

ppNormals <- post_predict(samples,n_cores=6, filter = "burn")


ciNormals <- plot_density(samples,layout=c(2,4),selection="mu", filter = "burn", do_plot=FALSE)
debug(mapped_par)
mapped_par(ciNormals[2,1:2], design1)

plot_fit(data1,ppNormals[[1]],layout=c(2,3))

source("emc/emc.R")

# debug(auto_burn)
burned <- auto_burn(sampler,nstart = 100, ndiscard = 20, cores_per_chain = 10, cores_for_chains = 1, verbose_run_stage = T, step_size = 50)
adapted <- auto_adapt(burned, cores_for_chains = 1, cores_per_chain = 10, verbose_run_stage = T)
sampled <- auto_sample(adapted, iter = 100, cores_for_chains = 1, cores_per_chain = 10, verbose_run_stage = T)
