require(coda)
require(mvtnorm)
require(MCMCpack)
source("emc/objects.R")    # pmwg object manipulation + mcmc creation functions
source("emc/statistics.R") # chain statistics
source("emc/plotting.R")   # chain and data plotting
source("emc/data.R")       # data generation 
source("emc/map.R")        # parameter mapping
source("emc/likelihood.R") # likelihoods
source("emc/fitting.R")    # fitting automation
source("emc/design.R")     # make design matrix and view parameters & mapping
source("emc/priors.R")     # sampling from priors
source("emc/pmwg.R")       # interface with pmwg
