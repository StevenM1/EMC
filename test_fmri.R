rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/fMRI/Normal/normal_whitened.R")

# Load in the data
load("~/Documents/UVA/2022/EMC/Data/VM2011/fmri_dm_VM2011.RData")
load("~/Documents/UVA/2022/EMC/Data/VM2011/fmri_data_VM2011.RData")

# subset data to only have one roi
data_lSTR <- fmri_data[, c("lSTR", "subjects")]
design1 <- make_design_fmri(design_matrix = fmri_dm, data = data_lSTR, model = normal, whiten = T)

sampler <- make_samplers(list(data_lSTR), list(design1))

samples <- run_emc(sampler, nsample = 500, cores_per_chain = 8, cores_for_chains = 1, verbose_run_stage = T)
