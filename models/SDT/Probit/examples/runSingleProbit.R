rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")

load("probitIndividual.RData")
samples <- run_chains(samplers,iter=c(1000,0,0))
save(samples,file="probitIndividual.RData")
