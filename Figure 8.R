# Comparisons of model predictions with Urine MEP data under no-hood DEP exposure from Weschler et al., (2015): 250 ug/m^3 exposed for 6 hours body top exposed to air.

library(deSolve)
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library (DescTools)
library(ggExtra)
library(gridExtra)
library(ggpubr)


#Exposure Parameters
exposure_duration = 6 #h

oral_dose = 0 #ug

#Physiolocal Parameters
Fgut = 0.0171
Fliver = 0.0257
Fplasma = 0.0424
Fgonads = 0.0005
Ffat	= 0.2142
Fskin	= 0.0371 
Frestbody = 0.91- Fliver- Fplasma- Fgut- Fgonads- Ffat- Fskin

QP=5*60			                          # ventilation rate (L/h)

FQgut = 0.178                             # based on <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844798/> (Table 2)
FQliver = 0.053                           # fractional hepatic liver blood flow
FQgonads = 0.0005                    		  # Fractional gonads blood flow for male
FQskin = 0.058                            # Fractional skin blood flow
FQfat = 0.052                         		# Fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin

ave_creatinine_prod_rate = 20 #mg/kg BW/day

FXskin =0.91     #0.82 (hood Weschler), 0.91(no hood Weschler), 0.04 (hood Krais), 0.13 (no hood Krais)

CLGint =991  
CLLint =1642 

#Chemical-specific Parameters
P_gut_plasma_DEP		=5.92
P_liver_plasma_DEP		=5.95
P_gonads_plasma_DEP		=3.26

P_fat_plasma_DEP		=64.42

P_plasma_air_DEP		= 5e7

P_gut_plasma_MEP		=1.88	#This is a parameter we added to account for MEP gut blood flow.
P_liver_plasma_MEP		=1.89
P_gonads_plasma_MEP	=1.16
P_skin_plasma_MEP		=2.25
P_fat_plasma_MEP		=17.42
P_restbody_plasma_MEP	=0.75

fu_DEP			= 0.30	#0.30		  # fu_DEP is the fraction unbound constant for DEP in plasma
fu_MEP			=0.598 #0.598		# fu_MEP is the fraction unbound constant for MEP in plasma


k_glucuronidation		= 1000 #1000 #k_glucuronidation is rate constant of MEP glucuronidation in the liver

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
         
         # Aplasma is the amount of chemical in plasma
         
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
               Aliver  = state.time6hr(X)[ 2],
               Afat  = state.time6hr(X)[ 3],
               Agonads  = state.time6hr(X)[4 ],
               Askin  = state.time6hr(X)[ 5],
               Aabsorbed = state.time6hr(X)[6 ],
               Arestbody  = state.time6hr(X)[7 ],
               Aplasma  = state.time6hr(X)[ 8],
               Aexhaled = state.time6hr(X)[9 ],
               Ainhaled = state.time6hr(X)[ 10],
               Agut_MEP  = state.time6hr(X)[11 ],
               Aliver_MEP  =state.time6hr(X)[12 ],
               Afat_MEP  = state.time6hr(X)[13 ],
               Agonads_MEP  = state.time6hr(X)[ 14],
               Askin_MEP  = state.time6hr(X)[15 ],
               Arestbody_MEP  = state.time6hr(X)[16 ],
               Aplasma_MEP  = state.time6hr(X)[ 17],
               Aurine_MEP  = state.time6hr(X)[18 ],
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
  
  return (list("outdf.1"  = outdf.1, 
               "outdf.2"  = outdf.2, 
               "outdf.3"  = outdf.3, 
               "outdf.4"  = outdf.4, 
               "outdf.5"  = outdf.5, 
               "outdf.6"  = outdf.6))
  
  
}

pred.DEP(theta.MCMC.DEP)

##----------------------------- create dataset ---------------------------------------- ##

Human.DEP.MCMC<- readRDS(file = "Human.DEP.MCMC.rds")
Wes.data.nohood.individual <- readxl::read_excel("Weschler data.xlsx",4)

predf.DEP      = pred.DEP(Human.DEP.MCMC$bestpar[-which_sig.DEP])

predf.DEP.1   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.1$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==1)$`UrineMEP(ug/g)`,
  participant="P1")

predf.DEP.2   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.2$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==2)$`UrineMEP(ug/g)`,
  participant="P2")

predf.DEP.3   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.3$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==3)$`UrineMEP(ug/g)`,
  participant="P3")

predf.DEP.4   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.4$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==4)$`UrineMEP(ug/g)`,
  participant="P4")

predf.DEP.5   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.5$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==5)$`UrineMEP(ug/g)`,
  participant="P5")

predf.DEP.6   = cbind.data.frame(
  pre.value   = predf.DEP$outdf.6$cre_adjusted_Auriune_MEP,
  obs.value   = subset(Wes.data.nohood.individual,subject==6)$`UrineMEP(ug/g)`,
  participant="P6")

predf.DEP = rbind(predf.DEP.1,predf.DEP.2,predf.DEP.3,
                  predf.DEP.4,predf.DEP.5,predf.DEP.6)


### model fit 

predf.DEP$log.obs = log10(predf.DEP$obs.value)
predf.DEP$log.pre = log10(predf.DEP$pre.value)

fit <- lm(log.obs ~ log.pre, data=predf.DEP)                   # Fit the model
predf.DEP$residuals = residuals(fit)                           # Save the residual values
predf.DEP$predicted = predict(fit)                             # Save the predicted values
predf.DEP$OPratio = predf.DEP$pre.value/predf.DEP$obs.value    # Estimated the predicted-to-observed ratio

##---------------------------------- plot ------------------------------------- ##
# par(ask=F)
# quartz()

#### ----- Figure 8a: predicted vs.observed plot ----- ####
p1 <- 
  ggplot(predf.DEP, aes(log.obs, log.pre)) + 
  geom_segment(aes(xend    = log.obs,                   # connect the actual data points with their corresponding predicted value    
                    yend    = predicted),
                    alpha   = .2) +
  geom_point(aes(shape   = as.factor(participant)),size = 4)  +
  geom_abline(intercept = 0, 
               slope     = 1,
               color     ="black",size = 1)+
  annotation_logticks() +
  scale_y_continuous(limits = c(1,4), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(1,4),labels = scales::math_format(10^.x))+ 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="#f7f7f7"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  xlab("Observed value")+
  ylab("Predicted value")

p1

#### -----  Figure 8b: predicted-to-observed vs. prediction plot ----- ####
p2 <-
  ggplot(predf.DEP, aes(log.pre, log10(OPratio))) +
  geom_hline(yintercept = log10(2),linetype = 3,color   = "red", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3,color   = "red", size =1) +
  #geom_ref_line(colour = "white", h = 0, size =0.5) +
  geom_point(color   = "lightblue", 
             aes(shape= as.factor(participant)),size = 3) +
  geom_smooth(se = FALSE) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(1,4),labels = scales::math_format(10^.x)) +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="#f7f7f7"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none')+
  xlab("Predicted value")+
  ylab("Predicted/Observed")

p3 <-ggMarginal(p2, type = "histogram", margins = "y",  
                yparams = list(binwidth = 0.1, fill = "lightblue"))

ggarrange(p1, p3, nrow = 1)

