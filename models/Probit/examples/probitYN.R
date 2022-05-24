rm(list=ls())
library(mvtnorm)
source("emc/emc.R")
source("models/Probit/probitYN.R")



# Two choice probit

# Not that for models of type="SDT" (such as probit) lR is automatically given
# a special "contr.increasing" contrast, and any terms involving the threshold
# threshold for the last level of lR are set to constants.

# assumes binary noise/signal factor and even number of confidence ratings
matchfun <- function(d) as.numeric(d$S) == (1+as.numeric(d$lR)>2)


# Simple EVSD 
designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ lR),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:2,matchfun=matchfun,
  constants=c(mean=0,sd=0), # noise mean = 0 and log(sd) = 0 (sd=1)
  model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 2         # signal mean
p_vector[2] <- log(1)    # signal sd 
p_vector[3] <- 1         # threshold, natural scale

mapped_par(p_vector,designPROBIT) 

dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=1000)
tapply(as.numeric(dataPROBIT$R)-1,dataPROBIT$S,mean)

# 3 level confidence
designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ lR),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,sd=0), # noise mean = 0 and log(sd) = 0 (sd=1)
  model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 1         # signal mean
p_vector[2] <- log(1.25)    # signal sd = 1.25 (treatment coding)
p_vector[3] <- -.5       # first threshold untransformed
# other thresholds exponentiated so > 0 then added to previous
# this is done in p_vector transform function
p_vector[4:7] <- log(rep(.5,4))  

mapped_par(p_vector,designPROBIT) 

dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)


plot_roc <- function(data,zROC=FALSE,main="",lim=NULL) {
  tab <- table(data$R,data$S)
  ctab <- 1-apply(tab/apply(tab,2,sum),2,cumsum)[-dim(tab)[1],] # p(Signal)
  if (!zROC) {
    if (!is.null(lim)) {xlim <- lim; ylim <- lim} else
                       {xlim <- c(0,1); ylim <- c(0,1)} 
    plot(ctab[,1],ctab[,2],xlab="p(FA)",ylab="p(H)",xlim=xlim,ylim=ylim,
         main=paste("ROC",main))
    lines(ctab[,1],ctab[,2])
    abline(a=0,b=1,lty=3)
  } else {
    ctab <- qnorm(ctab)
    if (is.null(lim)) lim <- c(min(ctab),max(ctab))
    plot(ctab[,1],ctab[,2],main=paste("zROC",main),
         xlab="z(FA)",ylab="z(H)",xlim=lim,ylim=lim)
     lines(ctab[,1],ctab[,2])
     abline(a=0,b=1,lty=3)
  }
}



par(mfrow=c(1,2))
plot_roc(dataPROBIT)
plot_roc(dataPROBIT,zROC=TRUE)

# 3 level confidence, factor A shifts threshold up

# Flist=list(mean ~ S, sd ~ S,threshold ~ A*lR)
# Ffactors=list(subjects=1,S=1:2,A=1:2);Rlevels=1:6
# constants=c(mean=0,sd=0);model=probit; Clist=NULL

designPROBIT <- make_design(Flist=list(mean ~ S, sd ~ S,threshold ~ A+lR),
  Ffactors=list(subjects=1,S=1:2,A=1:2),Rlevels=1:6, matchfun=matchfun,
  constants=c(mean=0,sd=0), # noise mean = 0 and log(sd) = 0 (sd=1)
  model=probit)

p_vector <- sampled_p_vector(designPROBIT)
p_vector[1] <- 1         # signal mean
p_vector[2] <- log(1)    # signal sd = 1.25 (treatment coding)
p_vector[3] <- -1        # first threshold untransformed
p_vector[4] <- log(.5)   # A2 shift up
# other thresholds exponentiated so > 0 then added to previous
# this is done in p_vector transform function
p_vector[5:8] <- log(rep(.5,4))  

mapped_par(p_vector,designPROBIT) 

dataPROBIT <- make_data(p_vector,design=designPROBIT,trials=10000)

par(mfrow=c(2,2))
plot_roc(dataPROBIT[dataPROBIT$A==1,],main="A=1")
plot_roc(dataPROBIT[dataPROBIT$A==1,],zROC=TRUE,main="A=1",lim=c(-1.5,2))
plot_roc(dataPROBIT[dataPROBIT$A==2,],main="A=2")
plot_roc(dataPROBIT[dataPROBIT$A==2,],zROC=TRUE,main="A=2",lim=c(-1.5,2))


tab <- table(dataPROBIT[,c("R","S","A")])
# p(Hit)
ctab <- 1-apply(tab/apply(tab,2:3,sum),2,cumsum)

par(mfrow=c(1,2))
plot(ctab[,1],ctab[,2],xlab="p(FA)",ylab="p(H)",main="ROC",xlim=c(0,1),ylim=c(0,1))
lines(ctab[,1],ctab[,2])
abline(a=0,b=1,lty=3)
lim <- c(min(qnorm(ctab)),max(qnorm(ctab)))
plot(qnorm(ctab[,1]),qnorm(ctab[,2]),main="zROC",
     xlab="z(FA)",ylab="z(H)",xlim=lim,ylim=lim)
lines(qnorm(ctab[,1]),qnorm(ctab[,2]))
abline(a=0,b=1,lty=3)


plot_defective_density(dataNORMAL,layout=c(1,2))

dadmNORMAL <- design_model(dataNORMAL,designNORMAL)
par(mfrow=c(1,3))
profile_pmwg(pname="mean",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmNORMAL)
profile_pmwg(pname="mean_lM1",p=p_vector,p_min=0,p_max=2,dadm=dadmNORMAL)
profile_pmwg(pname="sd",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmNORMAL)


# This example shows EMC works in a model without choice (Rlevels=1, no matchfun)
# Make a reasonably complex design. 
designNORMAL <- make_design(Flist=list(mean ~ A*B, sd ~ B),
  Ffactors=list(subjects=50,A=1:3,B=1:2),Rlevels=1,model=normal)

p_vector <- sampled_p_vector(designNORMAL)
p_vector[1:6] <- c(1,1,2,2,0,0) # Increasing in A and B, no interaction
p_vector[7:8] <- c(log(1,2)) # more variable for B2 than B1 


dataNORMAL <- make_data(p_vector,design=designNORMAL,trials=100)
plot_defective_density(dataNORMAL,layout=c(2,3),rt="topleft")

dadmNORMAL <- design_model(dataNORMAL,designNORMAL)
par(mfrow=c(2,4))
for (i in names(p_vector))
  profile_pmwg(pname=i,p=p_vector,p_min=p_vector[i]-.5,p_max=p_vector[i]+.5,dadm=dadmNORMAL)

samplers <- make_samplers(dataNORMAL,designNORMAL,type="standard",rt_resolution=.001)

