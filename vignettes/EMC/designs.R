### Examples of different EMC design and model combinations
rm(list=ls())
source("emc/emc.R")

# LBA Models

source("models/RACE/LBA/lbaB.R")
# Define LBA "B" parameterization to use making available a list called lbaB
# with elements 
# 0) type" "RACE"
# 1) p_types: names of different parameter types)
#             here ("v","sv","B","A","t0")
# 2) transform: transform applied to sampled parameter vector before mapping
#               here just the identity
# 3) Ntransform: takes parameter vector (p_vector) or matrix back to natural 
#                scale according to type defined by p_vector (letters in p_vector
#                or matrix column names before a "_", e.g. v_fred.
#                Here exp for all except v which is already on natural scale.
# 4) Ttransforms applied after mapping by design matrix and transformation to
#                natural scale, here b = B + A 
# 5) rfun: random function
# 6) dfun: Density function (PDF)
# 7) pfun: Probability function (CDF)
# 8) log_likelihood: log_likelihood_race function
# Also sources lba.R in same directory that defines rfun/dfun/pfun (in this
# case rLBA, dLBA and pLBA). 

# Define experiment design

# Match scoring, needed to create lM factor
matchfun <-  function(d)d$S==d$lR


ns <- 1 # Number of subjects, 1 here for initial examples

# "Formula factors": subjects factor and all factors used in design formulas that
# are defined as factors in the data frame to be fit.
# NB: must use "subjects", for others don't use "trial", "R", "rt", "lR" or "lM" 
#     or parameter names or "_" character in a column names.
Ffactors = list(subjects=1:ns,S=1:2) 

# "Formula list": Linear models for each p_type 
Flist=list(v ~ S*lM, sv ~ 1, B ~ lR, A ~ 1, t0 ~ 1) 

# Example of a custom contrast, used for lM
cmat <- matrix(c(-1/2,1/2),ncol=1)  

# Contrast list: default same for each p_type, but to make them differ specify 
# as a p_type named list. If not specified for a factor default contr.treatment 
Clist <- list(lM=cmat,S=contr.helmert,lR=contr.helmert) # can also use contr functions

# Response (R) factor levels
Rlevels = 1:2

# LBA scaling constant, can be anything named by sampled_p_vector below
constants <- c(sv=log(1))


# Clist is optional (default contr.treatment), matchfun is optional if Flist
# does not use latent match factor, lM
design <- make_design(Flist=Flist,Ffactors=Ffactors,Rlevels=Rlevels,
                      matchfun=matchfun,Clist=Clist,constants=constants,model=lbaB)

# Make parameters

# Make an appropriately named parameter vector to simulate data
# Get an empty vector (names corresponding to linear models)
p_vector <- sampled_p_vector(design,lbaB)
# natural scale v intercept, stimulus effect, match effect and their interaction 
p_vector[1:2] <- c(3,0)                # Match = 3, mismatch = 1              
p_vector[3:4] <- c(1,0)                # Natural scale match effect and stimulus effect
p_vector[5:6] <- c(log(3),log(1.05))   # log scale B intercept and bias (B/1.05, B*1.05)
p_vector[7:8] <- log(c(0.5,0.3))       # log scale A and t0 

# Data where every subject has identical parameters when p_vector is a vector,
# internally replicated to a matrix with one row per subject (with corresponding
# row names) and columns for each parameter, here
#   v v_S1 v_lM1 v_S1:lM1        B      B_lR1          A        t0
# 1 3    0     1        0 1.098612 0.04879016 -0.6931472 -1.203973

data <- make_data(p_vector=p_vector,design=design,model=lbaB,trials=2)
data

# Sample different parameters per subject

# Individual differences standard deviations 
sd_vector <- p_vector
sd_vector[1:length(sd_vector)] <- c(0.3,0.1,0.1,0.1, 0.15,0.1, 0.2, 0.2) 

# Independent variance-covariance matrix 
sigma <- diag(sd_vector^2)

# p_mat is subjects x parameters
p_mat <- rmvnorm(length(Ffactors$subjects),mean=p_vector,sigma=sigma)
dimnames(p_mat)[[1]] <- Ffactors$subjects
round(p_mat,2)

# Tiny data makes it easy to look at internal structure 
data <- make_data(p_vector=p_mat,design=design,model=lbaB,trials=2)
data

# How it works

# 1) Make "augmented" data (rows replicated for each accumulator).
da <- add_accumulators(data,matchfun)
da

# Make a compressed design matrix (with attribute "expand" 
# that can be used to make it into full form) 
# Contrasts can be either Clist functions or matrices (levels will be 
# automatically added as row names). 
dm <- make_dm(form=v ~ S*lM,da=da,Clist=Clist)
dm
dm[attr(dm,"expand"),] # twice as large
# Another contrast
dm <- make_dm(B ~ lR,da,Clist)
dm
dm[attr(dm,"expand"),] # 4x as large
# Intercept only
dm <- make_dm(A ~ 1,da) # Clist not needed for intercept only model
dm
dm[attr(dm,"expand"),,drop=FALSE] # 8x as large

mapped_par(p_vector,design,model=lbaB)
                       
# Add design and model to augmented data
# Note that default rounds RT to nearest 0.02 (argument rt_resolution=.02)
dadm <- design_model(data,design)
dadm

# This can be split up into a list with an entry for each subject
dadm_list <- dm_list(dadm)
dadm_list

# calculate likelihoods, same p_vector or all subjects
attr(dadm,"model")$log_likelihood(p_vector,dadm)

######### Profile plots to check it works for 1 subject with 20k trials

# Same design and p_vector as above but just 1 subject
design1 <- make_design(Flist=Flist,Ffactors=list(subjects=1,S=1:2),Rlevels=Rlevels,
                      matchfun=matchfun,Clist=Clist,constants=c(sv=log(1)),model=lbaB)

p_vector <- sampled_p_vector(design,lbaB)
p_vector[1:8] <- c(3,0,1,0,log(c(3,1.05,0.5,0.3)))

data1 <- make_data(p_vector,design1,lbaB,trials=2000)
plot_defective_density(data1,layout=c(1,2))

# Some examples of different compression levels
dadm1 <- design_model(data1,design1,rt_resolution=NULL)
dadm1 <- design_model(data1,design1,rt_resolution=.01)
dadm1 <- design_model(data1,design1) # Default .001 (ms) resolution

par(mfrow=c(2,5))
profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadm1)
profile_pmwg(pname="v_S1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm1)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm1)
profile_pmwg(pname="v_S1:lM1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm1)
profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm1)
profile_pmwg(pname="B_lR1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm1)
profile_pmwg(pname="A",p=p_vector,p_min=-1.5,p_max=0,dadm=dadm1)
profile_pmwg(pname="t0",p=p_vector,p_min=-1.5,p_max=-1,dadm=dadm1)



################ Three choice
design2 <- make_design(Flist=list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                      Ffactors=list(subjects=1,S=1:3),Rlevels=1:3,model=lbaB,
                      matchfun=matchfun,Clist=Clist,constants=c(sv=log(1)))

p_vector <- sampled_p_vector(design2,lbaB)
p_vector[1:2] <- c(3,1)   
p_vector[3:4] <- c(log(3),log(0.5)) # B and A
p_vector[5] <- log(c(0.3))   # log scale t0  

data2 <- make_data(p_vector,design2,lbaB,trials=1000)
plot_defective_density(data2,layout=c(1,3))

dadm2 <- design_model(data2,design2)
par(mfrow=c(2,3))
profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadm2)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm2)
profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm2)
profile_pmwg(pname="A",p=p_vector,p_min=-1.5,p_max=0,dadm=dadm2)
profile_pmwg(pname="t0",p=p_vector,p_min=-1.5,p_max=-1,dadm=dadm2)

################ Three choice with stimulus effect on v
Clist <- list(lM=cmat,S=contr.helmert,lR=contr.helmert) 

design3 <- make_design(Flist=list(v ~ S*lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                      Ffactors=list(subjects=1,S=1:3),Rlevels=1:3,model=lbaB,
                      matchfun=matchfun,Clist=Clist,constants=c(sv=log(1)))

p_vector <- sampled_p_vector(design3,lbaB)
p_vector[1:6] <- c(3,-1/2,1/2,1,0,0)
p_vector[7:8] <- c(log(3),log(0.5)) # B and A
p_vector[9] <- log(c(0.3))   # log scale t0

# Getting complicated to figure out parameters, this shows the mapping
mapped_par(p_vector,design3,lbaB)

data3 <- make_data(p_vector,design3,lbaB,trials=1000)
plot_defective_density(data3,layout=c(1,3))

dadm3 <- design_model(data3,design3)
par(mfrow=c(2,5))
profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadm3)
profile_pmwg(pname="v_S1",p=p_vector,p_min=-1,p_max=0,dadm=dadm3)
profile_pmwg(pname="v_S2",p=p_vector,p_min=0,p_max=1,dadm=dadm3)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm3)
profile_pmwg(pname="v_S1:lM1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm3)
profile_pmwg(pname="v_S2:lM1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm3)
profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm3)
profile_pmwg(pname="A",p=p_vector,p_min=-1.5,p_max=0,dadm=dadm3)
profile_pmwg(pname="t0",p=p_vector,p_min=-1.5,p_max=-1,dadm=dadm3)


################ Extra factor called FA

design4 <- make_design(Flist=list(v ~ S*lM, sv ~ 1, B ~ FA, A ~ 1, t0 ~ 1),
  Ffactors=list(subjects=1,FA=1:2,S=1:2),Rlevels=Rlevels,constants=c(sv=log(1)),
  matchfun=matchfun,Clist=list(lM=cmat,S=contr.helmert,lR=contr.helmert),model=lbaB )

p_vector <- sampled_p_vector(design4,lbaB)
p_vector[1:4] <- c(3,0,1,0)   
p_vector[5:7] <- c(log(2.5),log(1.4),log(0.5)) # B1 B1*1.4 and and A
p_vector[8] <- log(c(0.3))   # log scale t0  

# Note no contrasts specified for new factor A so uses default treatment coding
# B Intercept = first level of A = 2.5, 2nd level is intercept*1.4 = 3.5 
mapped_par(p_vector,design4,lbaB)

data4 <- make_data(p_vector,design4,lbaB,trials=1000)
plot_defective_density(data4,layout=c(2,2),xlim=c(0,3),rt="topleft")

dadm4 <- design_model(data4,design4)
par(mfrow=c(2,5))
profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadm4)
profile_pmwg(pname="v_S1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm4)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm4)
profile_pmwg(pname="v_S1:lM1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadm4)
profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm4)
profile_pmwg(pname="B_FA2",p=p_vector,p_min=0,p_max=0.5,dadm=dadm4)
profile_pmwg(pname="A",p=p_vector,p_min=-2,p_max=0,dadm=dadm4)
profile_pmwg(pname="t0",p=p_vector,p_min=-1.5,p_max=-1,dadm=dadm4)


################ Different t0 for each accumulator
design5 <- make_design(Flist=list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ lR),
                      Ffactors=list(subjects=1,S=1:2),Rlevels=Rlevels,model=lbaB,
                      matchfun=matchfun,Clist=Clist,constants=c(sv=log(1)))

p_vector <- sampled_p_vector(design5,lbaB)
p_vector[1:2] <- c(3,1)   
p_vector[3:4] <- c(log(3),log(0.5)) # B and A
p_vector[5:6] <- log(c(0.3,1.2))   # log scale t0 and t0_lR 

data5 <- make_data(p_vector,design5,lbaB,trials=1000)
plot_defective_density(data5,layout=c(1,2),xlim=c(0,3))

dadm5 <- design_model(data5,design5)
par(mfrow=c(2,4))
profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadm5)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm5)
profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadm5)
profile_pmwg(pname="A",p=p_vector,p_min=-1.5,p_max=0,dadm=dadm5)
profile_pmwg(pname="t0",p=p_vector,p_min=-1.5,p_max=-1,dadm=dadm5)
profile_pmwg(pname="t0_lR1",p=p_vector,p_min=0,p_max=0.4,dadm=dadm5)


#####  Multiple subjects

ns <- 3 
Ffactors = list(subjects=1:ns,S=1:2) 
Flist=list(v ~ S*lM, sv ~ 1, B ~ lR, A ~ 1, t0 ~ 1) 

design6 <- make_design(Flist=Flist,Ffactors=Ffactors,Rlevels=Rlevels,model=lbaB,
                      matchfun=matchfun,Clist=Clist,constants=c(sv=log(1)))

p_vector <- sampled_p_vector(design6,lbaB)
p_vector[1:8] <- c(3,0,1,0,log(c(1,3,1.05,.5,.3)))[-5] 

sd_vector <- p_vector
sd_vector[1:8] <- c(0.3,0.1,0.1,0.1,0.1, 0.15,0.1, 0.2, 0.2)[-5] 
sigma <- diag(sd_vector^2)
p_mat <- rmvnorm(length(Ffactors$subjects),mean=p_vector,sigma=sigma)
dimnames(p_mat)[[1]] <- Ffactors$subjects
head(round(p_mat,2))

# Range of sampled values
round(lbaB$Ntransform(apply(p_mat,2,range)),2)

# range of mapped parameters
mapped_p <- make_data(p_vector=p_mat,design=design6,model=lbaB,mapped_p=TRUE)
round(apply(mapped_p[names(mapped_p) %in% lbaB$p_types],2,range),2)

data6 <- make_data(p_mat,design6,lbaB,trials=10000)
par(mfcol=c(2,3))
sacc <- plot_defective_density(data6,correct_fun=function(data) data$S == data$R,xlim=c(0,3))
round(sort(sacc),2)

dadm6 <- design_model(data6,design6)
# individual subject likelihood
dadml <- dm_list(dadm6)
for (i in Ffactors$subjects)
 print(attr(dadm,"model")$log_likelihood(p_mat[i,],dadml[[i]]))


for (i in names(dadml)) {
  dadmi <- dadml[[i]]
  p_vector <- p_mat[i,]
  par(mfrow=c(2,5))
  profile_pmwg(pname="v",p=p_vector,p_min=2,p_max=4,dadm=dadmi,main=i)
  profile_pmwg(pname="v_S1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadmi,main=i)
  profile_pmwg(pname="v_lM1",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadmi,main=i)
  profile_pmwg(pname="v_S1:lM1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadmi,main=i)
  profile_pmwg(pname="B",p=p_vector,p_min=0.5,p_max=1.5,dadm=dadmi,main=i)
  profile_pmwg(pname="B_lR1",p=p_vector,p_min=-0.5,p_max=0.5,dadm=dadmi,main=i)
  profile_pmwg(pname="A",p=p_vector,p_min=-1.5,p_max=0,dadm=dadmi,main=i)
  profile_pmwg(pname="t0",p=p_vector,p_min=-1.75,p_max=-0.5,dadm=dadmi,main=i)
}


################ Extra factor with cell coding, single subject 
cmat2 <- diag(2)[,2:1]
Clist7 <- list(lM = cmat2, lR = cmat2, S = cmat2, FA = cmat2)
  
# Different constant to go with cell code mismatch sv = 1 (i.e., sv_lmFALSE)   
design7 <- make_design(Flist=list(v ~ 0 + lM, sv ~ 0 + lM, B ~ 0 + lR, A ~ 0 + S, t0 ~ 0 + FA),
  Ffactors=list(subjects=1,FA=1:2,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=matchfun,Clist=Clist7,constants=c(sv_lMFALSE=log(1)),model=lbaB)

p_vector <- sampled_p_vector(design7,lbaB)
p_vector[1:3] <- c(2,3,log(0.5))  # v mismatch, match, sv match   
p_vector[4:5] <- log(c(3,2.5))    # B left = 3, right = 2
p_vector[6:7] <- log(c(0.5,1))    # A for left and right 
p_vector[8:9] <- log(c(0.2,0.3))  # larger t0 for FA=2 than FA=1


# Simple cell parameters as expected
mapped_par(p_vector,design7,lbaB)

data7 <- make_data(p_vector,design7,lbaB,trials=1000)
plot_defective_density(data7,layout=c(2,2),xlim=c(0,3),rt="topleft")

dadm7 <- design_model(data7,design7)
par(mfrow=c(2,5))
profile_pmwg(pname="v_lMFALSE",p=p_vector,p_min=1,p_max=3,dadm=dadm7)
profile_pmwg(pname="v_lMTRUE",p=p_vector,p_min=2,p_max=4,dadm=dadm7)
profile_pmwg(pname="sv_lMTRUE",p=p_vector,p_min=log(0.25),p_max=log(0.75),dadm=dadm7)
profile_pmwg(pname="B_lRleft",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm7)
profile_pmwg(pname="B_lRright",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm7)
profile_pmwg(pname="A_Sleft",p=p_vector,p_min=log(0.25),p_max=log(1),dadm=dadm7)
profile_pmwg(pname="A_Sright",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadm7)
profile_pmwg(pname="t0_FA1",p=p_vector,p_min=log(.15),p_max=log(.25),dadm=dadm7)
profile_pmwg(pname="t0_FA2",p=p_vector,p_min=log(.25),p_max=log(.35),dadm=dadm7)

################ Redundant factor to data (FB, using data7 as a template)

Clist8 <- list(lM = cmat2, lR = cmat2, S = cmat2, FA = cmat2, FB = cmat2)

design8 <- make_design(Flist=list(v ~ 0 + lM, sv ~ 0 + lM, B ~ 0 + lR, A ~ 0 + S, t0 ~ 0 + FB),
  Ffactors=list(subjects=1,FA=1:2,FB=1:2,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=matchfun,Clist=Clist8,constants=c(sv_lMFALSE=log(1)),model=lbaB)

p_vector <- sampled_p_vector(design8,lbaB)
p_vector[1:3] <- c(2,3,log(0.5)) # v mismatch, match, sv match   
p_vector[4:5] <- log(c(3,2.5))   # B left = 3, right = 2
p_vector[6:7] <- log(c(0.5,1))   # A for left and right 
p_vector[8:9] <- log(c(0.2,1)) # MUCH Larger t0 for FB=2 (FA=1 and S=left or FA=2 and S=right) than FA=1


# Use data7 as a template
data8 <- data7
# Add interaction between FA and S and call it FB.
data8$FB <- factor(data8$FA==1 & data8$S=="left" | data8$FA==2 & data8$S=="right",
                   labels=1:2)
# make_data will fill in data8 with new R and RT based on formula
data8 <- make_data(p_vector,design8,lbaB,data=data8)
# Have to supply factors so doesnt include FB in the design
plot_defective_density(data8,layout=c(2,2),xlim=c(0,3),rt="topleft",factors=c("FA","S"))

# If we exclude FA can see the mixture
plot_defective_density(data8,layout=c(1,2),xlim=c(0,3),rt="topleft",factors=c("S"))


dadm8 <- design_model(data8,design8)
par(mfrow=c(2,5))
profile_pmwg(pname="v_lMFALSE",p=p_vector,p_min=1,p_max=3,dadm=dadm8)
profile_pmwg(pname="v_lMTRUE",p=p_vector,p_min=2,p_max=4,dadm=dadm8)
profile_pmwg(pname="sv_lMTRUE",p=p_vector,p_min=log(0.25),p_max=log(0.75),dadm=dadm8)
profile_pmwg(pname="B_lRleft",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm8)
profile_pmwg(pname="B_lRright",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm8)
profile_pmwg(pname="A_Sleft",p=p_vector,p_min=log(0.25),p_max=log(1),dadm=dadm8)
profile_pmwg(pname="A_Sright",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadm8)
profile_pmwg(pname="t0_FB1",p=p_vector,p_min=log(.15),p_max=log(.25),dadm=dadm8)
profile_pmwg(pname="t0_FB2",p=p_vector,p_min=log(.95),p_max=log(1.05),dadm=dadm8)


################ Continuous covariate (again using data7 as a template)

# Include an additive CV effect on rates (e.g., trial to trial fluctuations in mental speed)
# In Ffactors specify as NULL because no levels, no need to mention in Clist
design9 <- make_design(Flist=list(v ~ 0 + lM + CV, sv ~ 0 + lM, B ~ 0 + lR, A ~ 0 + S, t0 ~ 0 + FA),
  Ffactors=list(subjects=1,FA=1:2,S=c("left","right")),Fcovariates = c("CV"),
  Rlevels=c("left","right"),
  matchfun=matchfun,Clist=Clist7,constants=c(sv_lMFALSE=log(1)),model=lbaB)

p_vector <- sampled_p_vector(design9,lbaB)
p_vector[1:4] <- c(2,3,1,log(0.5))  # v mismatch, match, CV slope, sv match   
p_vector[5:6] <- log(c(3,2.5))    # B left = 3, right = 2
p_vector[7:8] <- log(c(0.5,1))    # A for left and right 
p_vector[9:10] <- log(c(0.2,0.3))  # larger t0 for FA=2 than FA=1


# Add covariate to existing data
CV <- round(runif(4000),2) # Increase between 0-1 x CV slope parameter
data9 <- cbind(data7,CV=CV)

# Simulate data
data9 <- make_data(p_vector,design9,lbaB,data=data9)
# Again have to ignore covariate
plot_defective_density(data9,layout=c(2,2),xlim=c(0,3),rt="topleft",factors=c("FA","S"))

dadm9 <- design_model(data9,design9)
par(mfrow=c(2,5))
profile_pmwg(pname="v_lMFALSE",p=p_vector,p_min=1,p_max=3,dadm=dadm9)
profile_pmwg(pname="v_lMTRUE",p=p_vector,p_min=2,p_max=4,dadm=dadm9)
profile_pmwg(pname="v_CV",p=p_vector,p_min=0.75,p_max=1.25,dadm=dadm9)
profile_pmwg(pname="sv_lMTRUE",p=p_vector,p_min=log(0.25),p_max=log(0.75),dadm=dadm9)
profile_pmwg(pname="B_lRleft",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm9)
profile_pmwg(pname="B_lRright",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadm9)
profile_pmwg(pname="A_Sleft",p=p_vector,p_min=log(0.25),p_max=log(1),dadm=dadm9)
profile_pmwg(pname="A_Sright",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadm9)
profile_pmwg(pname="t0_FA1",p=p_vector,p_min=log(.15),p_max=log(.25),dadm=dadm9)
profile_pmwg(pname="t0_FA2",p=p_vector,p_min=log(.25),p_max=log(.35),dadm=dadm9)

################ DDM
source("models/DDM/DDM/ddmTZD.R")

designDDM <- make_design(
  Flist=list(v~0+S, a~1,sv~1, t0~1, st0~1, s~1, Z~1, SZ~1, DP~1),
  Ffactors=list(subjects=1,S=1:2),Rlevels=1:2,Clist=list(S=diag(2)),
  constants=c(s=log(1)),model=ddmTZD)

p_vector <- sampled_p_vector(design=designDDM,model=ddmTZD)
p_vector[1:2] <- c(-1,1) # rates for lower (S=1) and upper (S=2) choices
p_vector[3:6] <- log(c(2,1,.3,.1)) # a, sv, t0, st0
p_vector[7:9] <- qnorm(c(0.5,0.2,0.5))   # Z, SZ, DP 

# Large number of trials but makes profiles accurate enough
dataDDM <- make_data(p_vector,design=designDDM,model=ddmTZD,trials=5000)
plot_defective_density(dataDDM,layout=c(1,2),xlim=c(0,3))

dadmDDM <- design_model(data=dataDDM,design=designDDM) 
# head(ddmTZD$log_likelihood(p_vector,dadmDDM))
par(mfrow=c(2,5))
profile_pmwg(pname="v_S1",p=p_vector,p_min=-1.5,p_max=-.5,dadm=dadmDDM)
profile_pmwg(pname="v_S2",p=p_vector,p_min=.5,p_max=1.5,dadm=dadmDDM)
profile_pmwg(pname="a",p=p_vector,p_min=log(1.5),p_max=log(2.5),dadm=dadmDDM)
profile_pmwg(pname="sv",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmDDM)
profile_pmwg(pname="t0",p=p_vector,p_min=log(.2),p_max=log(0.4),dadm=dadmDDM)
profile_pmwg(pname="st0",p=p_vector,p_min=log(.05),p_max=log(.15),dadm=dadmDDM)
profile_pmwg(pname="Z",p=p_vector,p_min=qnorm(.25),p_max=qnorm(.75),dadm=dadmDDM)
profile_pmwg(pname="SZ",p=p_vector,p_min=qnorm(.1),p_max=qnorm(.3),dadm=dadmDDM)
profile_pmwg(pname="DP",p=p_vector,p_min=qnorm(.4),p_max=qnorm(.6),dadm=dadmDDM)

################ Three choice RDM
source("models/RACE/RDM/rdmB.R")

# Diffusive variability, s, sets scale
designRDM <- make_design(Flist=list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s~1),
  Ffactors=list(subjects=1,S=1:3),Rlevels=1:3,matchfun=function(d)d$S==d$lR,
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert),
  model=rdmB,constants = c(s=log(1)))

p_vector <- sampled_p_vector(designRDM,rdmB)
p_vector[1:2] <- log(c(2,1.25)) # Intercept and (multiplicative) quality    
p_vector[3:4] <- c(log(3),log(0.5)) # B and A
p_vector[5] <- log(c(0.3))   # log scale t0  


dataRDM <- make_data(p_vector,design=designRDM,trials=1000)
plot_defective_density(dataRDM,layout=c(1,3))

dadmRDM <- design_model(dataRDM,designRDM)
par(mfrow=c(2,3))
profile_pmwg(pname="v",p=p_vector,p_min=log(1),p_max=log(3),dadm=dadmRDM)
profile_pmwg(pname="v_lM1",p=p_vector,p_min=log(1),p_max=log(1.5),dadm=dadmRDM)
profile_pmwg(pname="B",p=p_vector,p_min=log(2),p_max=log(4),dadm=dadmRDM)
profile_pmwg(pname="A",p=p_vector,p_min=log(.25),p_max=log(.75),dadm=dadmRDM)
profile_pmwg(pname="t0",p=p_vector,p_min=log(.2),p_max=log(.4),dadm=dadmRDM)

################ Three choice LNR
source("models/RACE/LNR/lnrMS.R")

#LNR intrinsically does not require scaling constant
designLNR <- make_design(Flist=list(m ~ lM, s ~ 1, t0 ~ 1),model=lnrMS,
  Ffactors=list(subjects=1,S=1:3),Rlevels=1:3,matchfun=function(d)d$S==d$lR,
  Clist=list(lM=matrix(c(-1/2,1/2),ncol=1),S=contr.helmert,lR=contr.helmert))

p_vector <- sampled_p_vector(designLNR,lnrMS)
p_vector[1:2] <- c(0,-1) # Intercept and quality of meanlog effectiely inverse
                         # of rate so match is LESS than mismatch
p_vector[3] <- c(log(1)) # sdlog
p_vector[4] <- log(c(0.3)) # log scale t0  


dataLNR <- make_data(p_vector,design=designLNR,model=lnrMS,trials=1000)
plot_defective_density(dataLNR,layout=c(1,3))

dadmLNR <- design_model(dataLNR,designLNR)
par(mfrow=c(2,2))
profile_pmwg(pname="m",p=p_vector,p_min=-.5,p_max=.5,dadm=dadmLNR)
profile_pmwg(pname="m_lM1",p=p_vector,p_min=-2,p_max=0,dadm=dadmLNR)
profile_pmwg(pname="s",p=p_vector,p_min=log(.5),p_max=log(1.5),dadm=dadmLNR)
profile_pmwg(pname="t0",p=p_vector,p_min=log(.15),p_max=log(.35),dadm=dadmLNR)

