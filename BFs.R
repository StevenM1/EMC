#### DDM ----

#### LBA ----

print(load("models/RACE/LBA/examples/samples/is2_lba_B.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bt0_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bv_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0_sv_NOa_n.RData"))

compare_MLL(list(B=is2_lba_B,Bvt0=is2_lba_Bvt0,Bt0sv=is2_lba_Bt0_sv,Bvsv=is2_sPNAS_Bv_sv,Bvt0sv=is2_lba_Bvt0_sv_NOa_n))
# Bvt0sv   Bvsv  Bt0sv   Bvt0      B 
#   0.60   0.32   0.08   0.00   0.00 

#### RDM ----

print(load("models/RACE/RDM/examples/samples/is2_rdm_B.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bv.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bvt0.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bt0.RData"))

compare_MLL(list(B=is2_rdm_B,Bv=is2_rdm_Bv,Bt0=is2_rdm_Bt0,Bvt0=is2_rdm_Bvt0))
#    B   Bv  Bt0 Bvt0 
# 0.82 0.08 0.08 0.01  

#### LNR ----

print(load("models/RACE/LNR/examples/samples/is2_lnr_mu.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM_t0.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_sElM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_mu_sElM_t0.RData"))

compare_MLL(list(mu=is2_lnr_mu,muslM=is2_lnr_mu_slM,muslMt0=is2_lnr_mu_slM_t0,
                 musElm=is2_lnr_mu_sElM,musElMt0=is2_lnr_mu_sElM_t0))
 # muslMt0    muslM   musElm musElMt0       mu 
 #    0.56     0.41     0.03     0.01     0.00
    
#### DIC/BPIC best ----

print(load("is2/is2_ddm_avt0_full.RData")) 
print(load("is2/is2_ddm_avt0_full_nocell.RData"))
print(load("is2/is2_lba_Bv_sv.RData"))
print(load("is2/is2_mu_sElM_t0.RData"))
print(load("is2/is2_rdm_Bvt0.RData"))


mll <- list(ddm=is2_ddm_avt0_full,ddm_cell=is2_sPNAS_avt0_full_nocell,
            lba=is2_sPNAS_Bv_sv,rdm=is2_rdm_Bvt0,lnr=is2_lnr_mu_sElM_t0)

compare_MLL(mll)
    #  ddm      lba ddm_cell      rdm      lnr 
    # 0.59     0.32     0.09     0.01     0.00
    
#### BF best

compare_MLL(list(lba=is2_lba_Bvt0_sv_NOa_n,
                 rdm=is2_rdm_B,lnr=is2_lnr_mu_slM_t0))
#  lba  rdm  lnr 
# 0.77 0.16 0.07 

### Simulation study ----

model_p <- function(mll,nboot=100000) 
  # mll is a list of vectors of marginal log-likelihoods for a set of models
  # picks a vector of mlls for each model in the list randomly with replacement
  # nboot times, calculates model probabilities and averages, the default
  # nboot seems good for stable results at 2 decimal places.
{
  pmp <- function(x) 
    # posterior model probability for a vector of marginal log-likelihoods  
  {
    x <- exp(x-max(x))
    x/sum(x)
  }
  
  sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    x[sample(length(x),nboot,replace=TRUE)]})),2,pmp),1,mean),decreasing=TRUE)
}

round(model_p(mll[-2]),2)
round(model_p(mll[-1]),2)




round(model_p(mll[1:2]),2)
round(model_p(mll[c(1,3)]),2)
round(model_p(mll[c(2,3)]),2)
round(model_p(mll[c(2,4)]),2)

sizeCheck <- function(mll1,mll2) {
  r40 <- numeric(24)
  for (i in 1:length(r40)) {
    indx <- (1+(i-1)*40):(i*40)
    r40[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r40),sd(r40)),2)) 
  r80 <- numeric(12)
  for (i in 1:length(r80)) {
    indx <- (1+(i-1)*80):(i*80)
    r80[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r80),sd(r80)),2)) 
  r120 <- numeric(8)
  for (i in 1:length(r120)) {
    indx <- (1+(i-1)*120):(i*120)
    r120[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r120),sd(r120)),2)) 
  r160 <- numeric(6)
  for (i in 1:length(r160)) {
    indx <- (1+(i-1)*160):(i*160)
    r160[i] <- model_p(list(mll1[[1]][indx],mll2[[1]][indx]))[1]
  }
  print(round(c(range(r160),sd(r160)),2)) 
}

sizeCheck(mll[1],mll[2])
sizeCheck(mll[1],mll[3])
sizeCheck(mll[2],mll[3])
sizeCheck(mll[2],mll[4])






