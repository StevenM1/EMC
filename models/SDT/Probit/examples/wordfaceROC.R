rm(list=ls())
source("emc/emc.R")
source("models/SDT/Probit/SDTgaussian.R")

print(load("Data/wordfaceROC.RData"))
# remove RT for ROC fitting
wordfaceROC$rt <- NA

wordfaceROC1 <- wordfaceROC[wordfaceROC$subjects=="104",]
wordfaceROC1$subjects <- factor(as.character(wordfaceROC1$subjects))

# To allow different thresholds for words and faces we must nest the lR factor
# within the FW factor (FW/lR). This ensures that thresholds are always 
# increasing within each level of FW. FW*lR will not enforce the increasing 
# constraint (but FW + lR will). If you wanted to allow different but increasing
# thresholds across combinations of factors, say A and B, use (A*B)/lR etc.
designFW1 <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=104,S=c("new","old"),FW=c("faces","words")),Rlevels=1:6, 
  matchfun=function(d) as.numeric(d$S) == (1+as.numeric(d$lR)>3),
  constants=c(mean=0,sd=0),model=probit)

p_vector <- sampled_p_vector(designFW1)
# constant mean new face = 0 
p_vector[1:3] <- c(-.5,1,1) # new words, old faces, old words 
# constant new face sd = 1
p_vector[4:6] <- log(c(1,1.25,1)) # old sd = 1.25, no item effect on sd 
p_vector[7:8] <- c(-0.5,-0.5)  # first threshold for faces, shift down by 0.5 for words
p_vector[9:16] <- log(rep(c(.5,.75),4)) # .5 threshold spacing for faces, 0.75 spacing for words

# Check mapping
mapped_par(p_vector,designFW1) 

# Simulate data
simFW1 <- make_data(p_vector,design=designFW1,trials=10000)
par(mfrow=c(2,2))
plot_roc(simFW1[simFW1$FW=="faces",],main="Faces")
plot_roc(simFW1[simFW1$FW=="faces",],zROC=TRUE,qfun=qnorm,main="Faces",lim=c(-2.5,2))
plot_roc(simFW1[simFW1$FW=="words",],main="Words")
plot_roc(simFW1[simFW1$FW=="words",],zROC=TRUE,qfun=qnorm,main="Words",lim=c(-2.5,2))

# profiles
dadmFWsim <- design_model(data=simFW1,design=designFW1)
par(mfrow=c(2,8))
for (i in names(p_vector))
  print(profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.25,p_max=p_vector[i]+.25,dadm=dadmFWsim))

samplers <- make_samplers(simFW1,designFW1,type="single")
save(samplers,file="probitFWsim.RData")
# runSingleProbitFWsim.R to get 1000 samples

samplers <- make_samplers(wordfaceROC1,designFW1,type="single")
save(samplers,file="probitFW1.RData")
# runSingleProbitFW1.R to get 1000 samples

# Look at simulation, clearly not close at 1000, added another 2000
print(load("probitFWsim.RData"))
plotChains(samples,subfilter=0,layout=c(4,4)) 
gd_pmwg(samples,subfilter=400)
plotACFs(samples,subfilter=400,layout=c(4,4))
round(100*es_pmwg(samples,subfilter=400)/1800) 
iat_pmwg(samples,subfilter=400) 
tabs <- plotDensity(samples,subfilter=400,layout=c(2,4),pars=p_vector)

# Look at data, quick convergence 
print(load("probitFW1.RData")) 
plotChains(samples,subfilter=100,layout=c(4,4)) # words parameters very autocorrelated
gd_pmwg(samples,subfilter=100)
plotACFs(samples,subfilter=100,layout=c(2,4))
iat_pmwg(samples,subfilter=100) 
tabs <- plotDensity(samples,subfilter=100,layout=c(4,4)) # returning prior


designFW1 <- make_design(Flist=list(mean ~ FW*S, sd ~ FW*S,threshold ~ FW/lR),
  Ffactors=list(subjects=levels(wordfaceROC$subjects),S=levels(wordfaceROC$S),
                FW=levels(wordfaceROC$FW)),Rlevels=1:6, 
  matchfun=function(d) as.numeric(d$S) == (as.numeric(d$lR)>3),
  constants=c(mean=0,sd=0),model=probit)

