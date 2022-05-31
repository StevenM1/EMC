rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")

load("probitFW1.RData")
samples <- run_chains(samplers,iter=c(2000,0,0))
save(samples,file="probitFW1.RData")
