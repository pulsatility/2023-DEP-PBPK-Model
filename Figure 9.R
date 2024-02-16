# Comparison of PBPK model simulations with experimental data of # Krais 2018: 300 ug/m^3 exposed for 3 hours wearing clothes.

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


#Exposure Parameters
exposure_duration = 3                    # h

oral_dose = 0                            # ug

#Physiological Parameters
BW =75                                   # kg

Fgut = 0.0171
Fliver = 0.0257
Fplasma = 0.0424
Fgonads = 0.0005
Ffat	= 0.2142
Fskin	= 0.0371 
Frestbody = 0.91- Fliver- Fplasma- Fgut- Fgonads- Ffat- Fskin

QP=5*60			                              # ventilation rate (L/h) (O for hood)
QC = 15*BW^0.74                           # (L/h) cardiac blood output

FQgut = 0.178                             # based on <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844798/> (Table 2)
FQliver = 0.053                           # fractional hepatic liver blood flow
FQgonads = 0.0005                    		  # Fractional gonads blood flow for male
FQskin = 0.058                            # Fractional skin blood flow
FQfat = 0.052                         		# Fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin

Vgut= Fgut*BW
Vliver= Fliver*BW
Vplasma= Fplasma*BW
Vgonads= Fgonads*BW
Vfat= Ffat*BW
Vskin= Fskin*BW
Vrestbody= Frestbody*BW

ave_creatinine_prod_rate = 20             # mg/kg BW/day
creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000)

Skin_area =1.91  #1.91 total skin area
FXskin =0.04     #0.82 (hood Weschler), 0.91(no hood Weschler), 0.04 (hood Krais), 0.13 (no hood Krais)

#Fraction of skin exposed
kgut_abs <-0  # assuming half of DEP is absorbed in 20 mins
CLGint =991 #991  
CLLint =1642 #1642

Skin_area_expo = Skin_area*FXskin

#Chemical-specific Parameters
P_gut_plasma_DEP		=5.92
P_liver_plasma_DEP		=5.95
P_gonads_plasma_DEP		=3.26
P_fat_plasma_DEP		=64.42
P_plasma_air_DEP		= 5e7

P_gut_plasma_MEP		=1.88	                 # This is a parameter we added to account for MEP gut blood flow.
P_liver_plasma_MEP		=1.89
P_gonads_plasma_MEP	=1.16
P_skin_plasma_MEP		=2.25
P_fat_plasma_MEP		=17.42
P_restbody_plasma_MEP	=0.75

fu_DEP			= 0.30	#0.30		               # fu_DEP is the fraction unbound constant for DEP in plasma
fu_MEP			=0.598 #0.598		               # fu_MEP is the fraction unbound constant for MEP in plasma

k_glucuronidation		= 1000                 # k_glucuronidation is rate constant of MEP glucuronidation in the liver

#scaled Physiolocal Parameters for MCMC with individual data
parameters <- c(
  #Exposure Parameters
  
  oral_dose = oral_dose, 
  
  #Physiolocal Parameters
  Vgut= Vgut,
  Vliver= Vliver,
  Vplasma= Vplasma,
  Vgonads= Vgonads,
  Vfat= Vfat,
  Vskin= Vskin,
  Vrestbody= Vrestbody,
  
  Qgut= FQgut* QC,
  Qliver= FQliver* QC,
  Qgonads= FQgonads* QC,
  Qskin= FQskin* QC,
  Qfat= FQfat* QC,
  Qrestbody= FQrestbody* QC,
  
  ave_creatinine_prod_rate = 20,                                 # mg/kg BW/day
  creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000), # g/h
  
  Skin_area_expo = Skin_area*FXskin,
  
  
  #Chemical-specific Parameters
  P_gut_plasma_DEP		=5.92,
  P_liver_plasma_DEP		=5.95,
  P_gonads_plasma_DEP		=3.26,
  P_fat_plasma_DEP		=64.42,
  P_plasma_air_DEP		= 5e7,
  
  P_gut_plasma_MEP		=1.88,	                                   # This is a parameter we added to account for MEP gut blood flow.
  P_liver_plasma_MEP		=1.89,
  P_gonads_plasma_MEP	=1.16,
  P_skin_plasma_MEP		=2.25,
  P_fat_plasma_MEP		=17.42,
  P_restbody_plasma_MEP	=0.75,
  
  fu_DEP			= 0.30,	#0.30		                                   # fu_DEP is the fraction unbound constant for DEP in plasma
  fu_MEP			=0.598, #0.598		                                 # fu_MEP is the fraction unbound constant for MEP in plasma
  
  CLLint =CLLint,
  CLGint=CLGint,
  
  k_glucuronidation		= 1000                                     # k_glucuronidation is rate constant of MEP glucuronidation in the liver
)



## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)


#use kurine_MEP
theta.MCMC.DEP <-log(c(

  P_skin_plasma_DEP		=7.29,
  P_restbody_plasma_DEP	=1.66,
  kurine_MEP		=0.25,                             # 1/h
  kskin_absorption=500,

  sig2  = 0.5,                                     # Model error (residuals); mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
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
         CI = CI_Air; # DEP concentration in inhaled air, ug/L 
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
         
         
         list(c(	 dAgut,
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


pred.DEP <- function (pars.DEP){
  
  pars.data <- function(dose){
  return(c(CI_Air = dose,parameters,lapply(pars.DEP[-which_sig.DEP],exp)))
  }
  CI_Air <- 0.3 # DEP concentration in inhaled air, ug/L 
  times <- seq(0, exposure_duration, by = 0.1)
  
  out<- lsoda(y = state, times = times, func = DEPMEP, parms = pars.data(CI_Air), jacfunc = NULL)
  
  out.dataframe <- as.data.frame(out)
  state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]
  
  state.time.2.vector <- c( Agut  =state.time.2[1],
                            Aliver  = state.time.2[2],
                            Afat  = state.time.2[3],
                            Agonads  = state.time.2[4],
                            Askin  = state.time.2[5],
                            Aabsorbed = state.time.2[6],
                            Arestbody  = state.time.2[7],
                            Aplasma  = state.time.2[8],
                            Aexhaled = state.time.2[9],
                            Ainhaled = state.time.2[10],
                            Agut_MEP  = state.time.2[11],
                            Aliver_MEP  = state.time.2[12],
                            Afat_MEP  = state.time.2[13],
                            Agonads_MEP  = state.time.2[14],
                            Askin_MEP  = state.time.2[15],
                            Arestbody_MEP  = state.time.2[16],
                            Aplasma_MEP  = state.time.2[17],
                            Aurine_MEP  = state.time.2[18],
                            Aglucuronidated_MEP =state.time.2[19]
                            
  )
  
  CI_Air <- 0
  times.2 <- seq(exposure_duration,60 , by = 0.01)
  out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = pars.data(CI_Air), jacfunc=NULL )
  out.final <-as.data.frame(rbind(out,out.2[-1,]))
  out.final$cre_adjusted_Auriune_MEP <- out.final$dAurine_MEP/creatinine_prod_rate #ug/g
  
  return (out.final[,c("time","cre_adjusted_Auriune_MEP")])
}

pred.DEP(theta.MCMC.DEP)


##------------------------------------------------ create dataset ------------------------------------------------ ##

Newtime.h       = pred.DEP(theta.MCMC.DEP)$time
nrow.h          = length (Newtime.h)

## Create the matrix 
MC.human    = matrix(nrow = nrow.h, ncol = 500) 

Human.DEP.MCMC<- readRDS(file = "Human.DEP.MCMC.rds")

## Read in Krais data
#krais no hood
krais.data.nohood.se<-read.table("Krais - No Hood average-ug-g_se.txt",sep = "\t",col.names = c("time","UrineMEP","se_concentration","se_time"),
                                 fill = TRUE,strip.white = TRUE,na.strings = "#NA")
#krais hood
krais.data.hood.se<-read.table("Krais - Hood average-ug-g_se.txt",sep = "\t",col.names = c("time","UrineMEP","se_concentration","se_time"),
                               fill = TRUE,strip.white = TRUE,na.strings = "#NA")

## Input paramters
for(i in 1:500){
  
  j = i *100  
  pars.human              = Human.DEP.MCMC$pars[j,]
  
  MCdata                  = pred.DEP(pars.human)
  MC.human[,i]  = MCdata$cre_adjusted_Auriune_MEP
  
  cat("iteration = ", i , "\n")
}

# no hood
write.csv(MC.human, file = "MC.human.csv")
# hood
write.csv(MC.human, file = "MC.human.hood.csv")


##------------------------------------------------ load the datasets ------------------------------------------------ ##
# no hood
MC.krais.nohood <- read.csv("MC.human.csv")[,-1]
# hood
MC.krais.hood <- read.csv("MC.human.hood.csv")[,-1]

MC.krais.nohood.data <- cbind(
  Time = Newtime.h, 
  as.data.frame(t(apply(log10(MC.krais.nohood), 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.krais.hood.data <- cbind(
  Time = Newtime.h, 
  as.data.frame(t(apply(log10(MC.krais.hood), 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)


##------------------------------------------------ plot ------------------------------------------------ ##
##------- Figure9 A --------##
# hood
p.MC.krais.hood <- 
  ggplot() + 
  geom_ribbon(data = MC.krais.hood.data, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue", alpha=0.3) +
  geom_line(data= MC.krais.hood.data, aes(x = Time, y = median_est), 
            color="blue",size=1) +
  geom_point(data = krais.data.hood.se, aes(x=time, y= log10(UrineMEP)),
             shape=1,size=2)+
  geom_errorbar(data = krais.data.hood.se, aes(x=time,
                                                 ymin= log10(UrineMEP-se_concentration), 
                                                 ymax = log10(UrineMEP+se_concentration)), size = 0.5,width=0.5,
                position=position_dodge(0.05),colour = "black")+
  labs (x ="Time (h)", y="Log10 Urine MEP/Creatinine (ug/g)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  ylim(c(-0.5,2))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 12),
                       axis.text.x = element_text(size = 12,colour="black"),
                       axis.text.y = element_text(size = 12,colour="black"))
p.MC.krais.hood

##------- Figure9 B --------##
# no hood
p.MC.krais.nohood <- 
  ggplot() + 
  geom_ribbon(data = MC.krais.nohood.data , aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue", alpha=0.3) +
  geom_line(data= MC.krais.nohood.data, aes(x = Time, y = median_est), 
            color="blue",size=1) +
  geom_point(data = krais.data.nohood.se, aes(x=time, y= log10(UrineMEP)),
             shape=1,size=2)+
  geom_errorbar(data = krais.data.nohood.se, aes(x=time,
                                                 ymin= log10(UrineMEP-se_concentration), 
                                                 ymax = log10(UrineMEP+se_concentration)), size = 0.5,width=0.5,
                position=position_dodge(0.05),colour = "black")+
  labs (x ="Time (h)", y="Log10 Urine MEP/Creatinine (ug/g)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  ylim(c(0,3.5))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12,colour="black"),
        axis.text.y = element_text(size = 12,colour="black"))

p.MC.krais.nohood
