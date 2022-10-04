rm(list=ls())
source("emc/emc.R")

############### EXG
source("models/RACE/StopSignal/SSexGaussian.R")

# NB: Stop *must* be first level of Rlevels
designEXG <- make_design(Flist=list(mu ~ lM, sigma ~ 1, tau ~ 1, muS~1, sigmaS~1, tauS~1, gf~ 1, tf~1),
  Ffactors=list(subjects=1,S=c("left","right")),Rlevels=c("stop","left","right"),
  matchfun=function(d)as.character(d$S)==as.character(d$lR),
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert),
  Fcovariates="SSD",model=SSexGaussian)


p_vector <- sampled_p_vector(design=designEXG,model=SSexGaussian)
p_vector[1:2] <- log(c(.6,.75)) # Intercept and (multiplicative) quality    
p_vector[3:4] <- c(log(.05),log(0.1)) # sigma and tau
p_vector[5:7] <- c(log(.2),log(.05),log(0.1)) # muS, sigmaS and tauS
p_vector[8:9] <- c(qnorm(0),qnorm(0))
p_vector[8:9] <- c(qnorm(.05),qnorm(0))
p_vector[8:9] <- c(qnorm(0),qnorm(.1))
p_vector[8:9] <- c(qnorm(.05),qnorm(.1))

# Staircase run where SSD is NA, starting point, step, and min and max SSD 
# values controlled by defaults set in model 
dataEXG <- make_data(p_vector,design=designEXG,trials=10000)
# Mixture of three cases
dataEXG$SSD <- rep(c(rep(Inf,300),rep(c(.1,.2,.3),100),rep(NA,400)),times=2)
# All staircase
dataEXG$SSD <- rep(NA,dim(dataEXG)[1])
dataEXG <- make_data(p_vector,design=designEXG,data=dataEXG)

# Need to specify factors to avoid including SSD
plot_defective_density(dataEXG,layout=c(1,2),factors="S")

# 30 subjects
designEXG <- make_design(Flist=list(mu ~ lM, sigma ~ 1, tau ~ 1, muS~1, sigmaS~1, tauS~1, gf~ 1, tf~1),
  Ffactors=list(subjects=1:30,S=c("left","right")),Rlevels=c("stop","left","right"),
  matchfun=function(d)as.character(d$S)==as.character(d$lR),
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert),
  Fcovariates="SSD",model=SSexGaussian)


ss <- length(designEXG$Ffactors$subjects)
p_mat <- t(matrix(
  rnorm(length(p_vector)*ss,mean=rep(p_vector,ss),sd=rep(abs(p_vector)/10,ss)),
  ncol=ss,dimnames=list(names(p_vector),1:ss)))
par(mfrow=c(2,5))
for (i in names(p_vector)) {
  if (i %in% c("tf","gf")) p <- pnorm(p_mat[,i]) else p <- exp(p_mat[,i])
  hist(p,breaks="fd",main=paste(i,"=",round(mean(p),3)))
}


dataEXG <- make_data(p_vector,design=designEXG,trials=200,
  Fcovariates=data.frame(SSD=rep(c(rep(Inf,150),rep(NA,50)),times=60))) # 25% staircae

SSD <- dataEXG$SSD[is.finite(dataEXG$SSD)]
mean(SSD)
mean(is.finite(dataEXG$SSD) & is.na(dataEXG$rt)) / mean(is.finite(dataEXG$SSD))
plot(SSD[1:400],type="l")

test_exgss <- make_samplers(dataEXG,designEXG)
save(test_exgss,file="test_exgss.RData")

dadmEXG <- design_model(data=dataEXG,design=designEXG)
# log_likelihood_race_ss(p_vector,dadm=dadmEXG)
par(mfrow=c(2,5))
profile_pmwg(pname="mu",p=p_vector,p_min=log(.25),p_max=log(.75),dadm=dadmEXG,cores=10)
profile_pmwg(pname="mu_lM1",p=p_vector,p_min=log(.5),p_max=log(1),dadm=dadmEXG,cores=10)
profile_pmwg(pname="sigma",p=p_vector,p_min=log(.025),p_max=log(.1),dadm=dadmEXG,cores=10)
profile_pmwg(pname="tau",p=p_vector,p_min=log(.05),p_max=log(.15),dadm=dadmEXG,cores=10)
profile_pmwg(pname="muS",p=p_vector,p_min=log(.1),p_max=log(.3),dadm=dadmEXG,cores=10)
profile_pmwg(pname="sigmaS",p=p_vector,p_min=log(.025),p_max=log(.075),dadm=dadmEXG,cores=10)
profile_pmwg(pname="tauS",p=p_vector,p_min=log(.05),p_max=log(.15),dadm=dadmEXG,cores=10)
profile_pmwg(pname="tf",p=p_vector,p_min=qnorm(.05),p_max=qnorm(.15),dadm=dadmEXG,cores=10)
profile_pmwg(pname="gf",p=p_vector,p_min=qnorm(.025),p_max=qnorm(.1),dadm=dadmEXG,cores=10)

print(load("test_exgss.RData"))
plot_chains(test_exgss,filter="burn",subfilter=80)

################ RDEX
#NB: t0 should *not* differ between choice accumulators unless this is due to
#    encoding time (as t0 included in race)

source("models/RACE/StopSignal/SSrdexB.R")

# NB: Stop *must* be first level of Rlevels
designRDEX <- make_design(Flist=list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s~1, 
  muS~1,sigmaS~1,tauS~1,gf~1,tf~1),
  Ffactors=list(subjects=1,S=c("left","right")),Rlevels=c("stop","left","right"),
  matchfun=function(d)as.character(d$S)==as.character(d$lR),
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert),
  Fcovariates="SSD",constants = c(s=log(1)),model=SSrdexB)

p_vector <- sampled_p_vector(designRDEX,SSrdexB)
p_vector[1:2] <- log(c(4,1.25)) # Intercept and (multiplicative) quality    
p_vector[3:4] <- c(log(1),log(0.1)) # B and A
p_vector[5] <- log(c(0.3))   # log scale t0
p_vector[6:8] <- log(c(.2,.05,.1)) # log scale mu,sigma,tau
p_vector[9:10] <- c(qnorm(0),qnorm(0))
p_vector[9:10] <- c(qnorm(.05),qnorm(0))
p_vector[9:10] <- c(qnorm(0),qnorm(.1))
p_vector[9:10] <- c(qnorm(.05),qnorm(.1))

dataRDEX <- make_data(p_vector,design=designRDEX,trials=1000)
# Mix of all
dataRDEX$SSD <- rep(c(rep(Inf,300),rep(c(.1,.2,.3),100),rep(NA,400)),times=2)
# All staircase
dataEXG$SSD <- rep(NA,dim(dataEXG)[1])
dataRDEX <- make_data(p_vector,design=designRDEX,data=dataRDEX)

plot_defective_density(dataRDEX,layout=c(1,2),factors="S")

# 30 subjects
designRDEX <- make_design(Flist=list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s~1, 
  muS~1,sigmaS~1,tauS~1,gf~1,tf~1),
  Ffactors=list(subjects=1:30,S=c("left","right")),Rlevels=c("stop","left","right"),
  matchfun=function(d)as.character(d$S)==as.character(d$lR),
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert),
  Fcovariates="SSD",constants = c(s=log(1)),model=SSrdexB)


ss <- length(designRDEX$Ffactors$subjects)
sds <- abs(p_vector)/10; sds["B"] <- .1
p_mat <- t(matrix(
  rnorm(length(p_vector)*ss,mean=rep(p_vector,ss),sd=rep(sds,ss)),
  ncol=ss,dimnames=list(names(p_vector),1:ss)))
par(mfrow=c(2,5))
for (i in names(p_vector)) {
  if (i %in% c("tf","gf")) p <- pnorm(p_mat[,i]) else p <- exp(p_mat[,i])
  hist(p,breaks="fd",main=paste(i,"=",round(mean(p),3)))
}


dataRDEX <- make_data(p_vector,design=designRDEX,trials=200,
  Fcovariates=data.frame(SSD=rep(c(rep(Inf,150),rep(NA,50)),times=60))) # 25% staircae

SSD <- dataRDEX$SSD[is.finite(dataRDEX$SSD)]
mean(SSD)
mean(is.finite(dataRDEX$SSD) & is.na(dataRDEX$rt)) / mean(is.finite(dataRDEX$SSD))
plot(SSD[1:400],type="l")

test_rdex <- make_samplers(dataRDEX,designRDEX)
save(test_rdex,file="test_rdex.RData")



dadmRDEX <- design_model(data=dataRDEX,design=designRDEX)
# head(log_likelihood_race_ss(p_vector,dadm=dadmRDEX))
par(mfrow=c(2,5))
profile_pmwg(pname="v",p=p_vector,p_min=log(3),p_max=log(5),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=log(1),p_max=log(1.5),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="B",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="A",p=p_vector,p_min=log(.05),p_max=log(.15),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="t0",p=p_vector,p_min=log(.2),p_max=log(.4),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="muS",p=p_vector,p_min=log(.1),p_max=log(.3),dadm=dadmRDEX)
profile_pmwg(pname="sigmaS",p=p_vector,p_min=log(.025),p_max=log(.1),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="tauS",p=p_vector,p_min=log(.05),p_max=log(.15),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="tf",p=p_vector,p_min=qnorm(.05),p_max=qnorm(.15),dadm=dadmRDEX,cores=10)
profile_pmwg(pname="gf",p=p_vector,p_min=qnorm(.025),p_max=qnorm(.15),dadm=dadmRDEX,cores=10)

