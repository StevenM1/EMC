rm(list=ls())
source("emc/emc.R")
source("models/RACE/RDM/rdmB.R")

print( load("Data/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

add_trials <- function(dat) {
  n <- table(dat$subjects)
  dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}

dat <- add_trials(dat)
  

ses <- factor(paste(dat$subjects,dat$E,dat$S))
cellMean <- tapply(dat$rt,ses,mean)
levels(ses) <- cellMean[levels(ses)]
cellMean <- as.numeric(as.character(ses))
residual <- dat$rt-cellMean
# Half of the residuals
dat$MT <- residual/2
# at 20ms resolution
dat$MT <- round(5*dat$MT,1)/5

# ss <- levels(dat$subjects)[1]
# isin <- dat$subjects==ss
# plot(dat$trial[isin],residual[isin],type="l",xlab="Trials",ylab="Residuals",main=ss)
# lines(dat$trial[isin],dat$MT[isin],col="red")
# abline(h=0)
# sd(residual[isin])
# sd(dat$MT[isin])
# sd(residual[isin]-dat$MT[isin])
  
# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))  

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))

#### Models ----

design_B_MT <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Fcovariates = c("MT"),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,B~lR*E,A~1,t0~MT),
  model=rdmB)

rdm_B <- make_samplers(dat,design_B_MT,type="standard",rt_resolution=.02)
save(rdm_B,file="rdm_B_MT.RData")

print(load("models/RACE/RDM/examples/samples/rdmPNAS_B.RData"))

# Simulate 1000 times using posterior mean
pp100 <- post_predict(rdm_B,n_post=1000,use_par="mean",n_cores=19)
# get per-trial predictions
predtrials <- tapply(pp100$rt,pp100[,c("subjects","E","S","trials")],mean)
pt <- as.data.frame.table(predtrials)
pt <- pt[!is.na(pt$Freq),]
names(pt)[5] <- "rt"

#
dat <- dat[order(dat$subjects,dat$E,dat$S),]
pt <- pt[order(pt$subjects,pt$E,pt$S),]
