design_Bv_sv <- make_design(
  Ffactors=list(subjects=unique(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~E*lM,sv~lM,B~lR*E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)