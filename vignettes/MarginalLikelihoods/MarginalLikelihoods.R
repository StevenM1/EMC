rm(list=ls())
source("emc/emc.R")

# It takes quite a bit of compute to get the IS2 estimates. Here is the head of
# the script used to run 15*64 = 960 chains, each producing a single estimate.
#
# rm(list=ls())
# source("emc/emc.R")
# source("models/DDM/DDM/ddmTZD.R")
# n_cores=64
# 
# print(load("DDM/sPNAS_a.RData"))
# is2_ddm_a <- run_IS2(sPNAS_a, 
#     n_cores = n_cores,IS_samples = 15*n_cores, max_particles = 5000)$finished
# save(is2_ddm_a,file="is2_ddm_a.RData")
#  ...


# First lets apply marginal likelihoods to selecting a parameterization
# within each model type.

#### DDM ----

print(load("models/DDM/DDM/examples/samples/is2_ddm_a.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_a_full.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_at0_full.RData"))
# following is a replicate to check how much they differ
print(load("models/DDM/DDM/examples/samples/is2_ddm_av_full_1.RData"))
is2_ddm_av_full_1 <- is2_ddm_av_full
print(load("models/DDM/DDM/examples/samples/is2_ddm_av_full.RData"))
print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full.RData"))
# Same model as before but no cell coding
print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full_nocell.RData"))

# Results are a vector of 960 marginal log-likelihood estimates
head(is2_ddm_a)
length(is2_ddm_a)

# Remember the IC results ...
#               DIC  wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# avt0nocell -17039 0.126 -16806 0.885        233 -17271 -17449 -17504
# avt0       -17043 0.874 -16802 0.115        240 -17283 -17468 -17523
# at0        -16948 0.000 -16745 0.000        203 -17151 -17314 -17354
# av         -16410 0.000 -16210 0.000        199 -16609 -16777 -16808
# a          -16343 0.000 -16176 0.000        167 -16510 728446 -16677

# Now lets look at model comparison with marginal likelihoods (i.e., using 
# Bayes Factors). This is done using posterior model probabilities. 
# The compare_MLL function takes argument mll, a list of vectors of marginal 
# log-likelihoods for a set of models and picks a vector of mlls for each model 
# in the list randomly with replacement nboot times, calculates model 
# probabilities and averages, the default nboot seems good for stable results 
# at 2 decimal places.

compareMLL <- function(mll,nboot=100000,digits=2,print_summary=TRUE) {
  
  pmp <- function(x) 
    # posterior model probability for a vector of marginal log-likelihoods  
  {
    x <- exp(x-max(x))
    x/sum(x)
  }
  
  out <- sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    x[sample(length(x),nboot,replace=TRUE)]})),2,pmp),1,mean),decreasing=TRUE)
  print(round(out,digits))
  invisible(out)
}


compare_MLL(list(a=is2_ddm_a,afull=is2_ddm_a_full,at0full=is2_ddm_at0_full,
                 av=is2_ddm_av_full,avt0full=is2_ddm_avt0_full,
                 avt0fullnocell=is2_sPNAS_avt0_full_nocell))
      # avt0full        at0full avt0fullnocell             av              a          afull 
      #     0.65           0.24           0.10           0.01           0.00           0.00 
          
# We see that the full avt0 model wins as was the case for DIC, and that 
# cell coding vs. no cell coding now makes quite a difference. This is 
# because the same prior was used in both cases, but means different things in
# the two parameterizations. This was for illustration, in practice priors
# should be matched to paramerizations

# A little sanity check to see that IS2 sampling noise is not to great.
avrep <- list(av=is2_ddm_av_full,av=is2_ddm_av_full_1)
compare_MLL(avrep)    
#  av  av 
# 0.5 0.5 
# Good agreement for central tendency
lapply(avrep,median)
lapply(avrep,mean)
# Bigger differences in variability
lapply(avrep,sd)
lapply(avrep,IQR)


#### LBA ----

print(load("models/RACE/LBA/examples/samples/is2_lba_B.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bt0_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bv_sv.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bvt0_sv_NOa_n.RData"))

# Remember the IC results ...
#           DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# B      -15243    0 -15100 0.000        143 -15386 -15485 -15529
# Bvt0   -15946    0 -15711 0.000        236 -16182 -16336 -16418
# Bt0sv  -16870    0 -16667 0.000        203 -17073 -17206 -17276
# Bvsv   -17200    1 -16986 0.999        214 -17414 -17569 -17627
# Bvt0sv -17182    0 -16973 0.001        209 -17391 -17542 -17601

compare_MLL(list(B=is2_lba_B,Bvt0=is2_lba_Bvt0,Bt0sv=is2_lba_Bt0_sv,Bvsv=is2_sPNAS_Bv_sv,Bvt0sv=is2_lba_Bvt0_sv_NOa_n))
# Bvt0sv   Bvsv  Bt0sv   Bvt0      B 
#   0.60   0.32   0.08   0.00   0.00 
# Similar, but liked the more complex model

#### RDM (only models with s variable) ----

print(load("models/RACE/RDM/examples/samples/is2_rdm_B_s.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bv_s.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bvt0_s.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bt0_s.RData"))

# Remember the IC results ...
#          DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# Bs    -16717    0 -16512 0.000        205 -16922 -17063 -17127
# Bt0s  -16738    0 -16551 0.000        186 -16924 -17064 -17110
# Bvs   -16856    0 -16629 0.005        227 -17083 -17235 -17309
# Bvt0s -16878    1 -16640 0.995        238 -17116 -17276 -17355

compare_MLL(list(Bs=is2_rdm_B_s,Bvs=is2_rdm_Bv_s,Bt0s=is2_rdm_Bt0_s,Bvt0s=is2_rdm_Bvt0_s))
 #   Bs   Bvs  Bt0s Bvt0s 
 # 0.81  0.12  0.05  0.01
 
# !!! Big discrepancy, simpler models clearly preferred. It seems likley there are
#     some issues with the prior here.

#### LNR ----

print(load("models/RACE/LNR/examples/samples/is2_lnr_mu.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_slM_t0.RData"))
print(load("models/RACE/LNR/examples/samples/is2_lnr_mu_sElM.RData"))
print(load("models/RACE/LNR/examples/samples/is2_mu_sElM_t0.RData"))

# Remember the IC results ...
#            DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# mu      -15080    0 -14869     0        211 -15290 -15451 -15501
# mulM    -16840    0 -16633     0        207 -17047 -17217 -17254
# mulMt0  -17102    0 -16879     1        223 -17325 -17509 -17548
# muElM   -17053    0 -16809     0        245 -17298 -17501 -17542
# muElMt0 -17129    1 -16846     0        283 -17411 -17609 -17694

compare_MLL(list(mu=is2_lnr_mu,mulM=is2_lnr_mu_slM,mulMt0=is2_lnr_mu_slM_t0,
                 muElm=is2_lnr_mu_sElM,muElMt0=is2_lnr_mu_sElM_t0))
 # mulMt0    mulM   muElm muElMt0      mu 
 #   0.56    0.41    0.03    0.00    0.00 

# In this case agrees with BPIC.
    
#### Compare DIC/BPIC best models of each type ----

print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full.RData")) 
print(load("models/DDM/DDM/examples/samples/is2_ddm_avt0_full_nocell.RData"))
print(load("models/RACE/LBA/examples/samples/is2_lba_Bv_sv.RData"))
print(load("models/RACE/LNR/examples/samples/is2_mu_sElM_t0.RData"))
print(load("models/RACE/RDM/examples/samples/is2_rdm_Bvt0.RData"))

# Remember the IC results ...
#              DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# DDM_avt0  -17043    0 -16802     0        240 -17283 -17468 -17523
# LBA_Bvsv  -17188    1 -16964     1        224 -17412 -17564 -17636
# RDM_Bvt0  -16648    0 -16403     0        245 -16893 -17054 -17137
# LNRmulMt0 -17129    0 -16846     0        283 -17411 -17609 -17694

mll <- list(ddm=is2_ddm_avt0_full,ddm_cell=is2_sPNAS_avt0_full_nocell,
            lba=is2_sPNAS_Bv_sv,rdm=is2_rdm_Bvt0,lnr=is2_lnr_mu_sElM_t0)
compare_MLL(mll)
    #  ddm      lba ddm_cell      rdm      lnr 
    # 0.59     0.32     0.09     0.01     0.00 
# DDM overtakes the LBA    

#### Comparing BF best models

compare_MLL(list(ddm=is2_ddm_avt0_full,lba=is2_lba_Bvt0_sv_NOa_n,
                 rdm=is2_rdm_B,lnr=is2_lnr_mu_slM_t0))
#  ddm  lba  rdm  lnr 
# 0.50 0.45 0.04 0.02 

# Now results are more equivocal

#### No cell DDM 

compare_MLL(list(ddm=is2_sPNAS_avt0_full_nocell,lba=is2_lba_Bvt0_sv_NOa_n,
                 rdm=is2_rdm_B,lnr=is2_lnr_mu_slM_t0))
#  lba  ddm  rdm  lnr 
# 0.71 0.22 0.05 0.02 

# But LBA beats the no_cell DDM, again showing prior influence

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






