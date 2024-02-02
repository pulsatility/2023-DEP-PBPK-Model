library(deSolve)
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(readxl)

# Bayesian MCMC simulation results for human exposure to 250 μg/m3 DEP in room air for 6 h without wearing hood. 
# Weschler 2015: 250 ug/m^3 exposed for 6 hours top body naked

#Exposure Parameters
exposure_duration = 6                      # h
oral_dose = 0                              # ug

#Physiolocal Parameters

Fgut = 0.0171
Fliver = 0.0257
Fplasma = 0.0424
Fgonads = 0.0005
Ffat	= 0.2142
Fskin	= 0.0371 
Frestbody = 0.91- Fliver- Fplasma- Fgut- Fgonads- Ffat- Fskin

QP=5*60			                              # ventilation rate (L/h)

FQgut = 0.178                             # based on <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844798/> (Table 2)
FQliver = 0.053                           # fractional hepatic liver blood flow
FQgonads = 0.0005                    		  # Fractional gonads blood flow for male
FQskin = 0.058                            # Fractional skin blood flow
FQfat = 0.052                         		# Fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin

ave_creatinine_prod_rate = 20             # mg/kg BW/day

FXskin =0.91                              # 0.91(no hood Weschler), 

#Fraction of skin exposed
CLGint =991  
CLLint =1642 

#Chemical-specific Parameters
P_gut_plasma_DEP		=5.92
P_liver_plasma_DEP		=5.95
P_gonads_plasma_DEP		=3.26
P_fat_plasma_DEP		=64.42
P_plasma_air_DEP		= 5e7
P_gut_plasma_MEP		=1.88	                # This is a parameter we added to account for MEP gut blood flow.
P_liver_plasma_MEP		=1.89
P_gonads_plasma_MEP	=1.16
P_skin_plasma_MEP		=2.25
P_fat_plasma_MEP		=17.42
P_restbody_plasma_MEP	=0.75

fu_DEP			= 0.30		                    # fu_DEP is the fraction unbound constant for DEP in plasma
fu_MEP			=0.598 	                      # fu_MEP is the fraction unbound constant for MEP in plasma

k_glucuronidation		= 1000                # k_glucuronidation is rate constant of MEP glucuronidation in the liver

#scaled Physiolocal Parameters for MCMC with individual data
BW.DEP <- function(BW,Skin_area){
  return(c(Vgut= Fgut*BW,
           Vliver= Fliver*BW,
           Vplasma= Fplasma*BW,
           Vgonads= Fgonads*BW,
           Vfat= Ffat*BW,
           Vskin= Fskin*BW,
           Vrestbody= Frestbody*BW,
           
           QC = 15*BW^0.74,
           Qgut= FQgut* 15*BW^0.74, #R does not allow to use QC directly
           Qliver= FQliver* 15*BW^0.74,
           Qgonads= FQgonads* 15*BW^0.74,
           Qskin= FQskin* 15*BW^0.74,
           Qfat= FQfat* 15*BW^0.74,
           Qrestbody= FQrestbody* 15*BW^0.74,
           
           creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000),
           
           Skin_area_expo = Skin_area*FXskin
  ))
}


############################################# Model Calibration with MCMC ###################################################
## One data sets was used in model evaluation                                                                               #
## A1: Measured pooled Plasma from 1977 - 2006 (Haugh et al., 2009)                                                         #
#############################################################################################################################

## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)


theta.MCMC.DEP <-log(c( 
           P_skin_plasma_DEP		=7.29,
           P_restbody_plasma_DEP	=1.66,
           kurine_MEP		=0.25,
           kskin_absorption=500,
           
           sig2  = 0.5,                                           # Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
           
           sig_P_skin_plasma_DEP		=0.3,
           sig_P_restbody_plasma_DEP	=0.3,
           sig_kurine_MEP		=0.3,                            
           sig_kskin_absorption=0.3
           ))


which_sig.DEP <- grep("sig", names(theta.MCMC.DEP))

state <- c(
  Agut  = 0,					
  Aliver  = 0,
  Afat  = 0,
  Agonads  = 0,
  Askin  = 0,
  Aabsorbed = 0,
  Arestbody  = 0,
  Aplasma  = 0,
  Aexhaled = 0,
  Ainhaled = 0,
  Agut_MEP  = 0,
  Aliver_MEP  = 0,
  Afat_MEP  = 0,
  Agonads_MEP  = 0,
  Askin_MEP  = 0,
  Arestbody_MEP  = 0,
  Aplasma_MEP  = 0,
  Aurine_MEP  = 0,
  Aglucuronidated_MEP =0
)


DEPMEP <- function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         #Exposure Parameters
         CI = CI_Air;

         FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin;
         Qrestbody= FQrestbody* QC;

         Cplasma = Aplasma/Vplasma;    
         Cplasma_MEP = Aplasma_MEP/Vplasma;
         Cgut = Agut/Vgut;       
         
         RAMG = CLGint * (Cgut*fu_DEP/P_gut_plasma_DEP);        
         dAgut = Qgut*Cplasma - Qgut*Cgut/P_gut_plasma_DEP - RAMG;
      
         Cliver = Aliver/Vliver;
  
         RAML=CLLint * (Cliver*fu_DEP/P_liver_plasma_DEP)
         
         dAliver = Qliver*Cplasma + Qgut*Cgut/P_gut_plasma_DEP - (Qliver+Qgut)*Cliver/P_liver_plasma_DEP - RAML;

         Cfat = Afat/Vfat;
         dAfat = Qfat*Cplasma - Qfat*Cfat/P_fat_plasma_DEP;
         
         Cgonads = Agonads/Vgonads;
         dAgonads = Qgonads*Cplasma - Qgonads*Cgonads/P_gonads_plasma_DEP;
         
         Cskin = Askin/Vskin;
         dAskin = Qskin*Cplasma - Qskin*Cskin/P_skin_plasma_DEP + CI*Skin_area_expo*kskin_absorption;
         dAabsorbed = CI*Skin_area_expo*kskin_absorption
         
         Crestbody = Arestbody/Vrestbody;
         dArestbody = Qrestbody*Cplasma - Qrestbody*Crestbody/P_restbody_plasma_DEP;        
         
         dAplasma = Qfat*Cfat/P_fat_plasma_DEP + (Qliver+ Qgut)*Cliver/P_liver_plasma_DEP + Qgonads*Cgonads/P_gonads_plasma_DEP + Qskin *Cskin/P_skin_plasma_DEP + Qrestbody*Crestbody/P_restbody_plasma_DEP - QC*Cplasma + QP*(CI - Cplasma/P_plasma_air_DEP);
         
         dAexhaled = QP*Cplasma/P_plasma_air_DEP;        
         dAinhaled = QP*CI;
         
         Cgut_MEP = Agut_MEP/Vgut;
         dAgut_MEP = Qgut*Cplasma_MEP - Qgut*Cgut_MEP/P_gut_plasma_MEP + RAMG;        
         
         Cliver_MEP = Aliver_MEP/Vliver;
         RGlucuronidation = k_glucuronidation* Cliver_MEP * (fu_MEP/P_liver_plasma_MEP);
         
         dAliver_MEP = Qliver*Cplasma_MEP + Qgut*Cgut_MEP/P_gut_plasma_MEP - (Qliver+Qgut)*Cliver_MEP/P_liver_plasma_MEP + RAML  - RGlucuronidation;
         
         Cfat_MEP = Afat_MEP/Vfat;
         dAfat_MEP = Qfat*Cplasma_MEP - Qfat*Cfat_MEP/P_fat_plasma_MEP;
         
         Cgonads_MEP = Agonads_MEP/Vgonads;
         dAgonads_MEP = Qgonads*Cplasma_MEP - Qgonads*Cgonads_MEP/P_gonads_plasma_MEP;
         
         Cskin_MEP = Askin_MEP/Vskin;
         dAskin_MEP = Qskin*Cplasma_MEP - Qskin*Cskin_MEP/P_skin_plasma_MEP;
         
         Crestbody_MEP = Arestbody_MEP/Vrestbody;
         dArestbody_MEP = Qrestbody*Cplasma_MEP - Qrestbody*Crestbody_MEP/P_restbody_plasma_MEP;
                 
         dAplasma_MEP = Qfat*Cfat_MEP/P_fat_plasma_MEP + (Qliver+Qgut)*Cliver_MEP/P_liver_plasma_MEP + Qgonads*Cgonads_MEP/P_gonads_plasma_MEP + Qskin*Cskin_MEP/P_skin_plasma_MEP + Qrestbody*Crestbody_MEP/P_restbody_plasma_MEP - QC*Cplasma_MEP;
         
         dAurine_MEP = kurine_MEP*Aglucuronidated_MEP;
         
         dAglucuronidated_MEP = RGlucuronidation - kurine_MEP*Aglucuronidated_MEP
                 
         list(c(	  dAgut,
                   dAliver,
                   dAfat,
                   dAgonads,
                   dAskin,
                   dAabsorbed,
                   dArestbody,
                   dAplasma,
                   dAexhaled,
                   dAinhaled,
                   dAgut_MEP,
                   dAliver_MEP,
                   dAfat_MEP,
                   dAgonads_MEP,
                   dAskin_MEP,
                   dArestbody_MEP,
                   dAplasma_MEP,
                   dAurine_MEP,
                   dAglucuronidated_MEP),dAurine_MEP=dAurine_MEP,RAML=RAML, RAMG=RAMG,Cplasma=Cplasma);
       })
}

#BW
BW1 <- 99
BW2 <- 90
BW3 <- 63
BW4 <- 74
BW5 <- 80
BW6 <- 77

#Body surface area
Skin_area1 <-2.21 
Skin_area2 <-2.16 
Skin_area3 <-1.73 
Skin_area4 <-1.93 
Skin_area5 <-2.03 
Skin_area6 <-1.96 

mcmc.fun.DEP <- function (pars.DEP){

  pars.data <- function(BW,Skin_area,dose){
    return(c(CI_Air = dose,BW.DEP(BW,Skin_area),lapply(pars.DEP[-which_sig.DEP],exp)))
  }

  call_1_ODE_solver<- function(BW,Skin_area,dose){
    y <- lsoda(y = state, times = times, func = DEPMEP, parms = pars.data(BW,Skin_area,dose), jacfunc = NULL)
    return(y)
  }
  
  times <- seq(0, exposure_duration, by = 0.1)
  CI_Air <- 0.25
  
  out.1.1<- call_1_ODE_solver(BW1,Skin_area1,CI_Air)
  out.1.1.dataframe<- as.data.frame(out.1.1)
  
  out.1.2<- call_1_ODE_solver(BW2,Skin_area2,CI_Air)
  out.1.2.dataframe<- as.data.frame(out.1.2)
  
  out.1.3<- call_1_ODE_solver(BW3,Skin_area3,CI_Air)
  out.1.3.dataframe<- as.data.frame(out.1.3)
  
  out.1.4<- call_1_ODE_solver(BW4,Skin_area4,CI_Air)
  out.1.4.dataframe<- as.data.frame(out.1.4)
  
  out.1.5<- call_1_ODE_solver(BW5,Skin_area5,CI_Air)
  out.1.5.dataframe<- as.data.frame(out.1.5)
  
  out.1.6<- call_1_ODE_solver(BW6,Skin_area6,CI_Air)
  out.1.6.dataframe<- as.data.frame(out.1.6)
  
  state.time6hr <- function(X){
    return(as.numeric(tail(X,n=1))[-1])
  }
  
  state.time6hr.vector <- function(X){
                   return(c(  Agut  =state.time6hr(X)[1],
                              Aliver  = state.time6hr(X)[2],
                              Afat  = state.time6hr(X)[3],
                              Agonads  = state.time6hr(X)[4],
                              Askin  = state.time6hr(X)[5],
                              Aabsorbed = state.time6hr(X)[6],
                              Arestbody  = state.time6hr(X)[7],
                              Aplasma  = state.time6hr(X)[8],
                              Aexhaled = state.time6hr(X)[9],
                              Ainhaled = state.time6hr(X)[10],
                              Agut_MEP  = state.time6hr(X)[11],
                              Aliver_MEP  =state.time6hr(X)[12],
                              Afat_MEP  = state.time6hr(X)[13],
                              Agonads_MEP  = state.time6hr(X)[14],
                              Askin_MEP  = state.time6hr(X)[15],
                              Arestbody_MEP  = state.time6hr(X)[16],
                              Aplasma_MEP  = state.time6hr(X)[17],
                              Aurine_MEP  = state.time6hr(X)[18],
                              Aglucuronidated_MEP =state.time6hr(X)[ 19]
                               ))
    
  }
  
  call_2_ODE_solver <- function(X,BW,Skin_area,dose){
    return(lsoda(y = state.time6hr.vector(X) , times = times.2, func = DEPMEP, parms = pars.data(BW,Skin_area,dose), jacfunc=NULL ))
  }

  times.2 <- seq(exposure_duration,60 , by = 0.1)
  #terminating exposure
  CI_Air <- 0
  
  out.2.1<- call_2_ODE_solver(out.1.1.dataframe,BW1,Skin_area1,CI_Air)
  out.2.1.dataframe<- as.data.frame(out.2.1)
  outdf.1<-as.data.frame(rbind(out.1.1,out.2.1[-1,]))# -1: remove the first time point of out.ave.2 which is the same as the last time point as out.ave.1
  
  out.2.2<- call_2_ODE_solver(out.1.2.dataframe,BW2,Skin_area2,CI_Air)
  out.2.2.dataframe<- as.data.frame(out.2.2)
  outdf.2<-as.data.frame(rbind(out.1.2,out.2.2[-1,]))
  
  out.2.3<- call_2_ODE_solver(out.1.3.dataframe,BW3,Skin_area3,CI_Air)
  out.2.3.dataframe<- as.data.frame(out.2.3)
  outdf.3<-as.data.frame(rbind(out.1.3,out.2.3[-1,]))
  
  out.2.4<- call_2_ODE_solver(out.1.4.dataframe,BW4,Skin_area4,CI_Air)
  out.2.4.dataframe<- as.data.frame(out.2.4)
  outdf.4<-as.data.frame(rbind(out.1.4,out.2.4[-1,]))
  
  out.2.5<- call_2_ODE_solver(out.1.5.dataframe,BW5,Skin_area5,CI_Air)
  out.2.5.dataframe<- as.data.frame(out.2.5)
  outdf.5<-as.data.frame(rbind(out.1.5,out.2.5[-1,]))
  
  out.2.6<- call_2_ODE_solver(out.1.6.dataframe,BW6,Skin_area6,CI_Air)
  out.2.6.dataframe<- as.data.frame(out.2.6)
  outdf.6<-as.data.frame(rbind(out.1.6,out.2.6[-1,]))
  
  outdf.1$cre_adjusted_Auriune_MEP <- outdf.1[,21]/BW.DEP(BW1,Skin_area1)["creatinine_prod_rate"] # outdf.1[,21]:dAurine_MEP
  outdf.2$cre_adjusted_Auriune_MEP <- outdf.2[,21]/BW.DEP(BW2,Skin_area2)["creatinine_prod_rate"]
  outdf.3$cre_adjusted_Auriune_MEP <- outdf.3[,21]/BW.DEP(BW3,Skin_area3)["creatinine_prod_rate"]
  outdf.4$cre_adjusted_Auriune_MEP <- outdf.4[,21]/BW.DEP(BW4,Skin_area4)["creatinine_prod_rate"]
  outdf.5$cre_adjusted_Auriune_MEP <- outdf.5[,21]/BW.DEP(BW5,Skin_area5)["creatinine_prod_rate"]
  outdf.6$cre_adjusted_Auriune_MEP <- outdf.6[,21]/BW.DEP(BW6,Skin_area6)["creatinine_prod_rate"]
 
  Wes.data.nohood.individual <- readxl::read_excel("Weschler data.xlsx",4)
  
  #match the simulation sampling time with data sampling time
  outdf.1 = outdf.1 [which(round(outdf.1$time,1) %in% round(subset(Wes.data.nohood.individual,subject==1)$'time(h)',1)),]
  outdf.2 = outdf.2 [which(round(outdf.2$time,1) %in% round(subset(Wes.data.nohood.individual,subject==2)$'time(h)',1)),]
  outdf.3 = outdf.3 [which(round(outdf.3$time,1) %in% round(subset(Wes.data.nohood.individual,subject==3)$'time(h)',1)),]
  outdf.4 = outdf.4 [which(round(outdf.4$time,1) %in% round(subset(Wes.data.nohood.individual,subject==4)$'time(h)',1)),]
  outdf.5 = outdf.5 [which(round(outdf.5$time,1) %in% round(subset(Wes.data.nohood.individual,subject==5)$'time(h)',1)),]
  outdf.6 = outdf.6 [which(round(outdf.6$time,1) %in% round(subset(Wes.data.nohood.individual,subject==6)$'time(h)',1)),]
  
  #log model output of cre_adjusted_Auriune_MEP
  log.yhat      <-c(log(outdf.1$cre_adjusted_Auriune_MEP),
                    log(outdf.2$cre_adjusted_Auriune_MEP),
                    log(outdf.3$cre_adjusted_Auriune_MEP),
                    log(outdf.4$cre_adjusted_Auriune_MEP),
                    log(outdf.5$cre_adjusted_Auriune_MEP),
                    log(outdf.6$cre_adjusted_Auriune_MEP))
  
  ## log.Observed data
  log.y         <-c(log(subset(Wes.data.nohood.individual,subject==1)$`UrineMEP(ug/g)`),
                    log(subset(Wes.data.nohood.individual,subject==2)$`UrineMEP(ug/g)`),
                    log(subset(Wes.data.nohood.individual,subject==3)$`UrineMEP(ug/g)`),
                    log(subset(Wes.data.nohood.individual,subject==4)$`UrineMEP(ug/g)`),
                    log(subset(Wes.data.nohood.individual,subject==5)$`UrineMEP(ug/g)`),
                    log(subset(Wes.data.nohood.individual,subject==6)$`UrineMEP(ug/g)`)
  )

  
  # The method of likelihood
  sig2            <- as.numeric((exp(pars.DEP[which_sig.DEP][1])))
  log_likelihood  <- -2*sum ((dnorm (log.y,
                                     mean =log.yhat,
                                     sd = sqrt(sig2),
                                     log = TRUE))) #dnorm outputs log-likelihood

  
  return(log_likelihood)
  
}

mcmc.fun.DEP(theta.MCMC.DEP)


## Define the Prior distributions: either normal or lognormal distribution
## nomral distribution

Prior.DEP <- function(pars.DEP) {
  
  ## Population level
  
  # The likelihood for population mean (parameters)
  pars.data = exp(pars.DEP[-which_sig.DEP])
  sig  <- as.numeric(exp(pars.DEP[which_sig.DEP])[2:length(which_sig.DEP)])# Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric(exp(pars.DEP[which_sig.DEP][1]))                  # error variances from model residual
  
  mean           = exp(theta.MCMC.DEP[-which_sig.DEP])
  CV             = 0.5                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  sd             = mean*CV

  logμ=log(mean^2/sqrt(sd^2+mean^2))
  logsigma=sqrt(log(sd^2/mean^2+1))
  
  print(logsigma)
  
  # Calculate likelihoods of each parameters; P(u|M,S)
  prior_pars     = dlnorm(pars.data, 
                              meanlog  = logμ, sdlog  = logsigma)
  
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  CV.sig         = exp(theta.MCMC.DEP[which_sig.DEP])[2:length(which_sig.DEP)] # Singmal0
  alpha          = (2+1)/(CU^2)                                   # Shape parametrer of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # Prior distribution for population variances; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # Error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  

  
  ## individual level; P(theta|u,sigmal^2)
  mean_i        = prior_pars 
  sd_i_2           = prior_sig

  logμ_i=log(mean_i^2/sqrt(sd_i_2+mean_i^2))
  logsigma_i=sqrt(log(sd_i_2/mean_i^2+1))
  
  prior_pars_i     = dlnorm(prior_pars,  
                                 meanlog  = logμ_i, sdlog  = logsigma_i) 
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig    = log (prior_sig)
  log.pri.pars.i = log (prior_pars_i)
  log.pri.sig2   = log (prior_sig2)
  
  # Maximau likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}

Prior.DEP(theta.MCMC.DEP)

###### random select initial values for MCMC chains from prior distribution

mean           = exp(theta.MCMC.DEP[-which_sig.DEP])
CV             = 0.5                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
sd             = mean*CV


prior_pars_dataframe <-data.frame(rbind(mean-sd,mean,mean+sd))

colnames(prior_pars_dataframe) <- names(theta.MCMC.DEP[-which_sig.DEP])
rownames(prior_pars_dataframe) <- c("Chain1","Chain2","Chain3")

prior_pars_dataframe$sig2<- exp(theta.MCMC.DEP[which_sig.DEP]["sig2"])
prior_pars_dataframe$sig_P_skin_plasma_DEP<- exp(theta.MCMC.DEP[which_sig.DEP]["sig_P_skin_plasma_DEP"])
prior_pars_dataframe$sig_P_restbody_plasma_DEP<- exp(theta.MCMC.DEP[which_sig.DEP]["sig_P_restbody_plasma_DEP"])
prior_pars_dataframe$sig_kurine_MEP<- exp(theta.MCMC.DEP[which_sig.DEP]["sig_kurine_MEP"])
prior_pars_dataframe$sig_kskin_absorption<- exp(theta.MCMC.DEP[which_sig.DEP]["sig_kskin_absorption"])
prior_pars_dataframe

prior_pars_dataframe_log <- log(prior_pars_dataframe)
prior_pars_dataframe_log


###################### MCMC simulation
cl <- makeCluster(2)
registerDoParallel(cl)
on.exit(stopCluster(cl))
# start time
strt<-Sys.time()

system.time(
MCMC.DEP <- apply(prior_pars_dataframe_log, MARGIN = 1, FUN = function(prior_pars_dataframe_log) {
  
  modMCMC(     f             = mcmc.fun.DEP, 
               p             = prior_pars_dataframe_log,  ## Qiang:this the place where different chains may start from different initila pramaeter values
               niter         = 1200,                      ## iteration number  #default 500000
               jump          = 0.01,                      ## jump function generation new parameters distribution using covariate matrix
               prior         = Prior.DEP,                 ## prior function
               updatecov     = 50,                        ## Adaptative Metropolis
               var0          = NULL,                      ## initial model variance;
               #wvar0         = 0.01,                     ## "Weight" for the initial model variance
               ntrydr        = 2,                         ## Delayed Rejection
               burninlength  = 0,                         #250000 ## number of initial iterations to be removed from output.
               outputlength  = 1100)                      ## delayed rejection (RD)
})
)

#end time
print(Sys.time()-strt)

stopCluster(cl) 


## Performance four chains to check the convergences

MC.H.1.DEP = as.mcmc (MCMC.DEP[[1]]$pars)     # first  chain
MC.H.2.DEP = as.mcmc (MCMC.DEP[[2]]$pars)     # second chain
MC.H.3.DEP = as.mcmc (MCMC.DEP[[3]]$pars)     # third  chain

combinedchains.DEP = mcmc.list(MC.H.1.DEP,MC.H.2.DEP,MC.H.3.DEP) ## combine all chains

gelman.diag (combinedchains.DEP)          # Gelman convergence diagnosis
heidel.diag (combinedchains.DEP)          # covergence diagnosis/Heidelberger and Welch's convergence diagnostic
gelman.plot (combinedchains.DEP[,1:4])    # gelman plot

# Save the posterior parameters (95% CI)
quan.human = exp(summary(MC.H.1.DEP)$quantiles) 

saveRDS(MCMC.DEP[[1]],file ='Human.DEP.MCMC.rds')
save(combinedchains.DEP,file = "combinedchains.DEP.RData")
save(MCMC.DEP,file = "MCMC.DEP.RData")
write.csv(quan.human,file="Human.summary_pos.csv")
write.csv(MC.H.1.DEP,file="Human.pos.csv")



