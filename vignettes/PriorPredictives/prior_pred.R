rm(list=ls())
source("emc/emc.R")


prior_predict <- function(n,samples,model=NULL,n_cores=1)
  
{
  mu <- prior_samples(samples[[1]],"mu",n)
  sigma <- prior_samples(samples[[1]],"sigma",n)
  design <- attr(samples,"design_list")[[1]]
  if (is.null(model)) model <- attr(samples,"model_list")[[1]]
  # Just do for one subject as all identical
  design$Ffactors$subjects <- design$Ffactors$subjects[1]  
  do.call(rbind,mclapply(1:n,function(i){
    alpha <- setNames(
      as.vector(mvtnorm::rmvnorm(1, mean = mu[i,], sigma = sigma[,,i])), 
      samples[[1]]$par_names)
    make_data(alpha,design,model,trials=1)[,-1]
  },mc.cores=n_cores))
}

plot_prior_predict <- function(prp,slow=100,layout=NULL) {
  bad <- is.na(prp$rt)
  print(c(isNA=100*mean(bad)))
  prp <- prp[!bad,]
  bad <- prp$rt>slow
  print(c(isSlow=100*mean(bad)))
  prp <- prp[!bad,]
  print(c(n=dim(prp)[1]))
  prp$rt <- as.numeric(prp$rt) # kluge
  plot_defective_density(prp,layout=layout)
  invisible(prp)
}


source("models/RACE/LBA/lbaB.R")
print(load("models/RACE/LBA/examples/samples/sPNAS_Bvt0_sv_NOa_n.RData"))
prp_lba <- prior_predict(n=1000,samples=sPNAS_Bvt0_sv_NOa_n,n_cores=10)
plot_prior_predict(prp_lba,slow=100,layout=c(2,3))


source("models/RACE/RDM/rdmB.R")
print(load("models/RACE/RDM/examples/samples/rdmPNAS_Bvt0.RData"))
prp_rdm <- prior_predict(n=1000,samples=rdm_Bvt0,n_cores=10)
plot_prior_predict(prp,layout=c(2,3))

source("models/RACE/LNR/lnrMS.R")
print(load("models/RACE/LNR/examples/samples/lnrPNAS_mu_sElM_t0.RData"))
prp_lnr <- prior_predict(n=1000,samples=lnr_mu_sElM_t0,n_cores=10)
plot_prior_predict(prp,layout=c(2,3))

source("models/DDM/DDM/ddmTZD.R")

print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full.RData"))
prp_avt0_full <- prior_predict(n=10000,samples=sPNAS_avt0_full,model=ddmTZD,n_cores=32)
plot_prior_predict(prp_avt0_full,layout=c(2,3))
gc()


print(load("models/DDM/DDM/examples/samples/sPNAS_avt0_full_nocell.RData"))
prp_avt0_full_nocell <- prior_predict(n=10000,samples=sPNAS_avt0_full_nocell,model=ddmTZD,n_cores=32)
plot_prior_predict(prp_avt0_full_nocell,layout=c(2,3))
gc()

save(prp_avt0_full,prp_avt0_full_nocell,prp_lba,prp_lnr,prp_rdm,file="prp.RData")
