# Normal race

#### distribution functions


dNORMAL <- function(rt,pars)
  # density for single accumulator
{
  dnorm(rt,mean=pars[,"mean"],sd=pars[,"sd"])
}

