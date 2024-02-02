
library(ggplot2)     
library(bayesplot)   # Package for MCMC traceplot
library(coda)


##------------------------------------------------ Load data sets----------------------------------------- ##
load("combinedchains.DEP.RData")


theta.MCMC.DEP <-log(c( 
  P_skin_plasma_DEP		=7.29,
  P_restbody_plasma_DEP	=1.66,
  kurine_MEP		=0.25,  #0.25                             # 1/h
  kskin_absorption=500,
  
  sig2  = 0.5,                                           # Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  sig_P_skin_plasma_DEP		=0.3,
  sig_P_restbody_plasma_DEP	=0.3,
  sig_kurine_MEP		=0.3,                            
  sig_kskin_absorption=0.3
))


##------------------------------------------------ Figure 5A ------------------------------------------------ ##
color_scheme_set("mix-blue-red")
mcmc_trace (
  combinedchains.DEP,
  pars =names(theta.MCMC.DEP),
  size = 0.5) 


##------------------------------------------------ Figure 5B ------------------------------------------------ ##
## gelman plot
gelman.plot(combinedchains.DEP)          


