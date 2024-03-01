
# Global parameter sensitivity analyses on Plasma DEP and Gonad DEP, Plasma MEP, Gonad MEP, and Urine MEP for different kinetic metrics: Cmax, AUC, and half-lives of the first and second clearance phases for simulations of human exposure to 250 μg/m3 DEP in room air for 6 h wearing hood with body top exposed to air as in Weschler 2015.

library(deSolve)
library(sensitivity)
library(ggplot2)

#Exposure Parameters
exposure_duration = 6                    # h
CI_Air <- 0.25	                         # DEP concentration in inhaled air, ug/L

oral_dose = 0                            # ug

#Physiolocal Parameters
BW =75                                   # kg
Density_liver=1.051                      # kg/L https://www.ncbi.nlm.nih.gov/pubmed/3579513 
Density_gut=1.045                        # kg/L

Fgut = 0.0171
Fliver = 0.0257
Fplasma = 0.0424
Fgonads = 0.0005
Ffat	= 0.2142
Fskin	= 0.0371  
Frestbody = 0.91- Fliver- Fplasma- Fgut- Fgonads- Ffat- Fskin

QP= 0*5*60		                           	# ventilation rate (L/h)
QC = 15*BW^0.74                           # (L/h/kg) Cardiac output

FQgut = 0.178                             # based on <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844798/> (Table 2)
FQliver = 0.053                           # fractional hepatic liver blood flow
FQgonads = 0.0005                    		  # Fractional gonads blood flow for male
FQskin = 0.058                            # Fractional skin blood flow
FQfat = 0.052                         		# Fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin

Qgut= FQgut* QC
Qliver= FQliver* QC
Qgonads= FQgonads* QC
Qskin= FQskin* QC
Qfat= FQfat* QC
Qrestbody= FQrestbody* QC

Vgut= Fgut*BW
Vliver= Fliver*BW
Vplasma= Fplasma*BW
Vgonads= Fgonads*BW
Vfat= Ffat*BW
Vskin= Fskin*BW
Vrestbody= Frestbody*BW

Skin_area =1.9  
FXskin = 0.82                             # 0.82 (hood Weschler)
#area of skin exposed
Skin_area_expo = Skin_area*FXskin

kgut_abs <- 0                             # assuming half of DEP is absorbed in 20 mins

CLGint =991 #991  
CLLint =1642 #1642

ave_creatinine_prod_rate = 20             # mg/kg BW/day
creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000)

parameters <- c(
  #Exposure Parameters
  CI = CI_Air,
  
  oral_dose = oral_dose, 
  
  #Physiolocal Parameters
  Vgut= Vgut,
  Vliver= Vliver,
  Vplasma= Vplasma,
  Vgonads= Vgonads,
  Vfat= Vfat,
  Vskin= Vskin,
  Vrestbody= Vrestbody,
  
  QC = 15*BW^0.74,                                    # Cardiac output
  Qgut= FQgut* QC,
  Qliver= FQliver* QC,
  Qgonads= FQgonads* QC,
  Qskin= FQskin* QC,
  Qfat= FQfat* QC,
  Qrestbody= FQrestbody* QC,
  
  Skin_area_expo = Skin_area_expo,
  
  #parameters below used for sensitivity analysis
  kskin_absorption = 500, 	                          # Skin absorption rate constant, length/h

  #Chemical-specific Parameters
  P_gut_plasma_DEP		=5.92,
  P_liver_plasma_DEP		=5.95,
  P_gonads_plasma_DEP		=3.26,
  P_skin_plasma_DEP		=7.29, 
  P_fat_plasma_DEP		=64.42, 
  P_restbody_plasma_DEP	=1.66,	
  P_plasma_air_DEP		= 5e7,
  P_gut_plasma_MEP		=1.88,	                        # This is a parameter we added to account for MEP gut blood flow.
  P_liver_plasma_MEP		=1.89,
  P_gonads_plasma_MEP	=1.16,
  P_skin_plasma_MEP		=2.25,
  P_fat_plasma_MEP		=17.42, 
  P_restbody_plasma_MEP	=0.75,
  
  fu_DEP			= 0.30,			                            # fu_DEP is the fraction unbound constant for DEP in plasma
  fu_MEP			=0.598, 		                            # fu_MEP is the fraction unbound constant for MEP in plasma
  
  kurine_MEP		=0.25,                                # 1/h
  
  CLLint =CLLint,
  CLGint=CLGint,
  
  k_glucuronidation		= 1000                          # k_glucuronidation is rate constant of MEP glucuronidation in the liver 
)

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
         Cplasma = Aplasma/Vplasma;
         Cplasma_MEP = Aplasma_MEP/Vplasma;
         Cgut = Agut/Vgut;
         RAMG = CLGint * (Cgut*fu_DEP/P_gut_plasma_DEP);
         dAgut = Qgut*Cplasma - Qgut*Cgut/P_gut_plasma_DEP - RAMG 

         Cliver = Aliver/Vliver;
         RAML=CLLint * (Cliver*fu_DEP/P_liver_plasma_DEP);  
         dAliver = Qliver*Cplasma + Qgut*Cgut/P_gut_plasma_DEP - (Qliver+Qgut)*Cliver/P_liver_plasma_DEP - RAML;
         Cfat = Afat/Vfat;
         dAfat = Qfat*Cplasma - Qfat*Cfat/P_fat_plasma_DEP;
         
         Cgonads = Agonads/Vgonads;
         dAgonads = Qgonads*Cplasma - Qgonads*Cgonads/P_gonads_plasma_DEP;
         
         Cskin = Askin/Vskin;
         dAskin = Qskin*Cplasma - Qskin*Cskin/P_skin_plasma_DEP + CI*Skin_area_expo*kskin_absorption;
         dAabsorbed = CI*Skin_area_expo*kskin_absorption;
         
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
 
         list(c(	  
           dAgut,
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
           dAglucuronidated_MEP),dAurine_MEP=dAurine_MEP,RAML=RAML, RAMG=RAMG);
       })
}

#half life function
half_life_function <- function(x){
  return(out.2$time[min(which(x<max(x)/2))])
}

#half life 2 function
half_life_2_function <- function(x){
  
  half_life_2 <- log(2) / ((log(out.final[which(out.final[, 1]==40),][,x])-
                              log(out.final[which(out.final[, 1]==50),][,x]))/(50-40))
  
  return(half_life_2)
}



factors <- names(parameters)[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)] #remove oral dose, CI,BW, QC, creatinine prod rate, ave creatinine prod rate

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit

# Define the lower and upper limits that will be the input to the Morris function
binf <- parameters[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]*LL
bsup <- parameters[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]*UL

sample <- seq(from = 1, to = 1, by = 1)

for (i in 1:length(sample)) {
  set.seed(12345)
  
  x <- morris(model = NULL, factors = factors, r = sample[i],
              design = list(type = "oat", levels = 5, grid.jump = 1), 
              binf = binf, bsup = bsup, scale = TRUE)
  
  for (iteration in 1:nrow(x$X)) { 
    
    parameters=x$X[iteration, ]
    
    times <- seq(0,6,0.01)
    CI<- 0.25

    tmp <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters, jacfunc = NULL)

    parameters.2<-parameters
    CI <- 0
   
    out.dataframe <- as.data.frame(tmp)
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
    
    times.2 <- seq(6, 200, by = 0.01)
    out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
    out.final <-rbind(tmp,out.2[-1,])
    
    out.final <- as.data.frame(out.final)

    AUC <- apply(out.final[,-1], 2, sum)*0.1/2
    Cmax <- apply(out.final[,-1], 2, max)
    out.2 <- as.data.frame(out.2)
    half_life <- apply(out.2[,-1], 2, half_life_function) 
    half_life_2 <- sapply(names(out.final[,-1]),half_life_2_function)
    #half_life_2 <- apply(out.final[,-1],2,half_life_2_function)
    
    if (iteration == 1) { # initialize
      results = AUC
      results.2 = Cmax
      results.3=half_life
      results.4=half_life_2
      sampled.parms=parameters
      
    } else { # accumulate
      results = rbind(results, AUC) # this is the place we can change to AUC as opposed to one time point response 
      results.2 = rbind(results.2,Cmax)
      results.3=rbind(results.3,half_life)
      results.4=rbind(results.4,half_life_2)
      sampled.parms=rbind(sampled.parms,parameters)
      
    }
  }
  
  ###Aurine_MEP
  #AUC
  tell(x, results[ ,c("dAurine_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X1 <- apply(x$ee, 2, mean)
    mu.star1 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma1<- apply(x$ee, 2, sd)
  } else {
    X1 <- rbind(X1, apply(x$ee, 2, mean))
    mu.star1 <- rbind(mu.star1,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma1 <-rbind(sigma1,apply(x$ee, 2, sd)) 
  }
  
  #Cmax
  tell(x, results.2[ ,c("dAurine_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X1.2 <- apply(x$ee, 2, mean)
    mu.star1.2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma1.2<- apply(x$ee, 2, sd)
  } else {
    X1.2 <- rbind(X1.2, apply(x$ee, 2, mean))
    mu.star1.2 <- rbind(mu.star1.2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma1.2 <-rbind(sigma1.2,apply(x$ee, 2, sd)) 
  }
  
  #half life
  tell(x, results.3[ ,c("dAurine_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X1.3 <- apply(x$ee, 2, mean)
    mu.star1.3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma1.3<- apply(x$ee, 2, sd)
  } else {
    X1.3 <- rbind(X1.3, apply(x$ee, 2, mean))
    mu.star1.3 <- rbind(mu.star1.3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma1.3 <-rbind(sigma1.3,apply(x$ee, 2, sd)) 
  }
  
  #half life_2
  tell(x, results.4[ ,c("dAurine_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X1.4 <- apply(x$ee, 2, mean)
    mu.star1.4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma1.4<- apply(x$ee, 2, sd)
  } else {
    X1.4 <- rbind(X1.4, apply(x$ee, 2, mean))
    mu.star1.4<- rbind(mu.star1.4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma1.4 <-rbind(sigma1.4,apply(x$ee, 2, sd)) 
  }
  
  ###Aplasma
  #AUC
  tell(x, results[ ,c("Aplasma")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X2 <- apply(x$ee, 2, mean)
    mu.star2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma2 <- apply(x$ee, 2, sd)
  } else {
    X2 <- rbind(X2, apply(x$ee, 2, mean))
    mu.star2 <- rbind(mu.star2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma2<-rbind(sigma2,apply(x$ee, 2, sd)) 
  }
  
  # Cmax
  tell(x, results.2[ ,c("Aplasma")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X2.2 <- apply(x$ee, 2, mean)
    mu.star2.2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma2.2 <- apply(x$ee, 2, sd)
  } else {
    X2.2 <- rbind(X2.2, apply(x$ee, 2, mean))
    mu.star2.2 <- rbind(mu.star2.2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma2.2<-rbind(sigma2.2,apply(x$ee, 2, sd)) 
  }
  
  tell(x, results.3[ ,c("Aplasma")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X2.3 <- apply(x$ee, 2, mean)
    mu.star2.3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma2.3 <- apply(x$ee, 2, sd)
  } else {
    X2.3 <- rbind(X2.3, apply(x$ee, 2, mean))
    mu.star2.3 <- rbind(mu.star2.3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma2.3<-rbind(sigma2.3,apply(x$ee, 2, sd)) 
  }
  
  
  ##half life_2
  tell(x, results.4[ ,c("Aplasma")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X2.4 <- apply(x$ee, 2, mean)
    mu.star2.4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma2.4 <- apply(x$ee, 2, sd)
  } else {
    X2.4 <- rbind(X2.4, apply(x$ee, 2, mean))
    mu.star2.4 <- rbind(mu.star2.4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma2.4<-rbind(sigma2.4,apply(x$ee, 2, sd)) 
  }
  
  
  ### Aplasma_MEP
  #AUC
  tell(x, results[ ,c("Aplasma_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X3 <- apply(x$ee, 2, mean)
    mu.star3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma3 <- apply(x$ee, 2, sd)
  } else {
    X3 <- rbind(X3, apply(x$ee, 2, mean))
    mu.star3 <- rbind(mu.star3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma3<-rbind(sigma3,apply(x$ee, 2, sd)) 
  }
  
  #Cmax
  tell(x, results.2[ ,c("Aplasma_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X3.2 <- apply(x$ee, 2, mean)
    mu.star3.2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma3.2 <- apply(x$ee, 2, sd)
  } else {
    X3.2 <- rbind(X3.2, apply(x$ee, 2, mean))
    mu.star3.2 <- rbind(mu.star3.2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma3.2<-rbind(sigma3.2,apply(x$ee, 2, sd)) 
  }
  
  #half life
  tell(x, results.3[ ,c("Aplasma_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X3.3 <- apply(x$ee, 2, mean)
    mu.star3.3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma3.3 <- apply(x$ee, 2, sd)
  } else {
    X3.3 <- rbind(X3.3, apply(x$ee, 2, mean))
    mu.star3.3 <- rbind(mu.star3.3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma3.3<-rbind(sigma3.3,apply(x$ee, 2, sd)) 
  }
  
  
  #half life_2
  tell(x, results.4[ ,c("Aplasma_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X3.4 <- apply(x$ee, 2, mean)
    mu.star3.4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma3.4<- apply(x$ee, 2, sd)
  } else {
    X3.4 <- rbind(X3.4, apply(x$ee, 2, mean))
    mu.star3.4<- rbind(mu.star3.4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma3.4 <-rbind(sigma3.4,apply(x$ee, 2, sd)) 
  }
  
  
  ### Agonads
  #AUC
  tell(x, results[ ,c("Agonads")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X4 <- apply(x$ee, 2, mean)
    mu.star4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma4 <- apply(x$ee, 2, sd)
  } else {
    X4 <- rbind(X4, apply(x$ee, 2, mean))
    mu.star4 <- rbind(mu.star4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma4<-rbind(sigma4,apply(x$ee, 2, sd)) 
  }
  
  #Cmax
  tell(x, results.2[ ,c("Agonads")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X4.2 <- apply(x$ee, 2, mean)
    mu.star4.2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma4.2 <- apply(x$ee, 2, sd)
  } else {
    X4.2<- rbind(X4.2, apply(x$ee, 2, mean))
    mu.star4.2 <- rbind(mu.star4.2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma4.2<-rbind(sigma4.2,apply(x$ee, 2, sd)) 
  }
  
  #half life
  tell(x, results.3[ ,c("Agonads")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X4.3 <- apply(x$ee, 2, mean)
    mu.star4.3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma4.3 <- apply(x$ee, 2, sd)
  } else {
    X4.3<- rbind(X4.3, apply(x$ee, 2, mean))
    mu.star4.3 <- rbind(mu.star4.3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma4.3<-rbind(sigma4.3,apply(x$ee, 2, sd)) 
  }
  
  
  #half life
  tell(x, results.4[ ,c("Agonads")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X4.4 <- apply(x$ee, 2, mean)
    mu.star4.4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma4.4 <- apply(x$ee, 2, sd)
  } else {
    X4.4<- rbind(X4.4, apply(x$ee, 2, mean))
    mu.star4.4 <- rbind(mu.star4.4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma4.4<-rbind(sigma4.4,apply(x$ee, 2, sd)) 
  }
  
  
  ### Agonads_MEP
  #AUC
  tell(x, results[ ,c("Agonads_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X5<- apply(x$ee, 2, mean)
    mu.star5 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma5 <- apply(x$ee, 2, sd)
  } else {
    X5 <- rbind(X5, apply(x$ee, 2, mean))
    mu.star5 <- rbind(mu.star5,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma5<-rbind(sigma5,apply(x$ee, 2, sd)) 
  }
  
  #Cmax
  tell(x, results.2[ ,c("Agonads_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X5.2<- apply(x$ee, 2, mean)
    mu.star5.2 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma5.2 <- apply(x$ee, 2, sd)
  } else {
    X5.2<- rbind(X5.2, apply(x$ee, 2, mean))
    mu.star5.2 <- rbind(mu.star5.2,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma5.2<-rbind(sigma5.2,apply(x$ee, 2, sd)) 
  }
  
  #half life
  tell(x, results.3[ ,c("Agonads_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X5.3<- apply(x$ee, 2, mean)
    mu.star5.3 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma5.3 <- apply(x$ee, 2, sd)
  } else {
    X5.3<- rbind(X5.3, apply(x$ee, 2, mean))
    mu.star5.3 <- rbind(mu.star5.3,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma5.3<-rbind(sigma5.3,apply(x$ee, 2, sd)) 
  }
  
  #half life
  tell(x, results.4[ ,c("Agonads_MEP")]) # We focus on parameter effect on blood concentration
  
  if (i == 1){
    X5.4<- apply(x$ee, 2, mean)
    mu.star5.4 <- apply(x$ee, 2, function(x) mean(abs(x)))
    sigma5.4 <- apply(x$ee, 2, sd)
  } else {
    X5.4<- rbind(X5.4, apply(x$ee, 2, mean))
    mu.star5.4 <- rbind(mu.star5.4,apply(x$ee, 2, function(x) mean(abs(x))))
    sigma5.4<-rbind(sigma5.4,apply(x$ee, 2, sd)) 
  }
  
}

# set the order
order_sensitivity <-c("P_skin_plasma_MEP",
                                   "P_liver_plasma_DEP",
                                   "P_liver_plasma_MEP",
                                   "P_gut_plasma_DEP",
                                   "P_gut_plasma_MEP",
                                   "P_plasma_air_DEP",
                                   "P_gonads_plasma_MEP",
                                   "P_fat_plasma_MEP",
                                   "P_restbody_plasma_MEP",
                                   "fu_DEP",
                                   "kurine_MEP",
                                   "CLLint",            
                                   "CLGint",
                                   "fu_MEP", 
                                   "k_glucuronidation",
                                   "P_gonads_plasma_DEP",
                                   "P_fat_plasma_DEP",  
                                   "P_restbody_plasma_DEP",
                                   "P_skin_plasma_DEP",
                                   "kskin_absorption")

##------------------------------------------------ plot ------------------------------------------------ ##

###Aurine_MEP
#AUC (no rank)
X1 <- as.data.frame(X1)
X1$var <- rownames(X1)
X1$var <- factor(X1$var,levels = order_sensitivity)

p1.1<-ggplot(data=X1, aes(x=var,y=X1)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aurine_MEP AUC")
p1.1

#Cmax
X1.2<- as.data.frame(X1.2)
X1.2$var <- rownames(X1.2)
X1.2$var <- factor(X1.2$var,levels = order_sensitivity)

p1.2<-ggplot(data=X1.2, aes(x=var,y=X1.2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aurine_MEP Cmax")
p1.2

#half_life
X1.3<- as.data.frame(X1.3)
X1.3$var <- rownames(X1.3)
X1.3$var <- factor(X1.3$var,levels = order_sensitivity)

p1.3<-ggplot(data=X1.3, aes(x=var,y=X1.3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aurine_MEP Half life")
p1.3

#half_life.2
X1.4<- as.data.frame(X1.4)
X1.4$var <- rownames(X1.4)
X1.4$var <- factor(X1.4$var,levels = order_sensitivity)

p1.4<-ggplot(data=X1.4, aes(x=var,y=X1.4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aurine_MEP Half life.2")
p1.4


### Aplasma
#AUC
X2 <- as.data.frame(X2)
X2$var <- rownames(X2)
X2$var <- factor(X2$var,levels = order_sensitivity)

p2.1<-ggplot(data=X2, aes(x=var,y=X2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma AUC")
p2.1

#Cmax
X2.2<- as.data.frame(X2.2)
X2.2$var <- rownames(X2.2)
X2.2$var <- factor(X2.2$var,levels = order_sensitivity)

p2.2<-ggplot(data=X2.2, aes(x=var,y=X2.2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma Cmax")
p2.2

#half_life
X2.3<- as.data.frame(X2.3)
X2.3$var <- rownames(X2.3)
X2.3$var <- factor(X2.3$var,levels = order_sensitivity)

p2.3<-ggplot(data=X2.3, aes(x=var,y=X2.3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma Half life")
p2.3

#half_life.2
X2.4<- as.data.frame(X2.4)
X2.4$var <- rownames(X2.4)
X2.4$var <- factor(X2.4$var,levels = order_sensitivity)

p2.4<-ggplot(data=X2.4, aes(x=var,y=X2.4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma Half life.2")
p2.4


### Aplasma_MEP
#AUC 
X3 <- as.data.frame(X3)
X3$var <- rownames(X3)
X3$var <- factor(X3$var,levels = order_sensitivity)

p3.1<-ggplot(data=X3, aes(x=var,y=X3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma_MEP AUC")
p3.1

#Cmax
X3.2<- as.data.frame(X3.2)
X3.2$var <- rownames(X3.2)
X3.2$var <- factor(X3.2$var,levels = order_sensitivity)

p3.2<-ggplot(data=X3.2, aes(x=var,y=X3.2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma_MEP Cmax")
p3.2

#half_life
X3.3<- as.data.frame(X3.3)
X3.3$var <- rownames(X3.3)
X3.3$var <- factor(X3.3$var,levels = order_sensitivity)

p3.3<-ggplot(data=X3.3, aes(x=var,y=X3.3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma_MEP Half life")
p3.3

#half_life.2
X3.4<- as.data.frame(X3.4)
X3.4$var <- rownames(X3.4)
X3.4$var <- factor(X3.4$var,levels = order_sensitivity)

p3.4<-ggplot(data=X3.4, aes(x=var,y=X3.4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Aplasma_MEP Half life.2")
p3.4


### Agonads
#AUC 
X4 <- as.data.frame(X4)
X4$var <- rownames(X4)
X4$var <- factor(X4$var,levels = order_sensitivity)

p4.1<-ggplot(data=X4, aes(x=var,y=X4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads AUC")
p4.1

#Cmax
X4.2<- as.data.frame(X4.2)
X4.2$var <- rownames(X4.2)
X4.2$var <- factor(X4.2$var,levels = order_sensitivity)

p4.2<-ggplot(data=X4.2, aes(x=var,y=X4.2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads Cmax")
p4.2

#half_life
X4.3<- as.data.frame(X4.3)
X4.3$var <- rownames(X4.3)
X4.3$var <- factor(X4.3$var,levels = order_sensitivity)

p4.3<-ggplot(data=X4.3, aes(x=var,y=X4.3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads Half life")
p4.3

#half_life.2
X4.4<- as.data.frame(X4.4)
X4.4$var <- rownames(X4.4)
X4.4$var <- factor(X4.4$var,levels = order_sensitivity)

p4.4<-ggplot(data=X4.4, aes(x=var,y=X4.4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads Half life.2")
p4.4

### Agonads MEP

#AUC 
X5<- as.data.frame(X5)
X5$var <- rownames(X5)
X5$var <- factor(X5$var,levels = order_sensitivity)

p5.1<-ggplot(data=X5, aes(x=var,y=X5)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads MEP AUC")
p5.1

#Cmax
X5.2<- as.data.frame(X5.2)
X5.2$var <- rownames(X5.2)
X5.2$var <- factor(X5.2$var,levels = order_sensitivity)

p5.2<-ggplot(data=X5.2, aes(x=var,y=X5.2)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads MEP Cmax")
p5.2

#half_life
X5.3<- as.data.frame(X5.3)
X5.3$var <- rownames(X5.3)
X5.3$var <- factor(X5.3$var,levels = order_sensitivity)

p5.3<-ggplot(data=X5.3, aes(x=var,y=X5.3)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads MEP Half life")
p5.3

#half_life.2
X5.4<- as.data.frame(X5.4)
X5.4$var <- rownames(X5.4)
X5.4$var <- factor(X5.4$var,levels = order_sensitivity)

p5.4<-ggplot(data=X5.4, aes(x=var,y=X5.4)) +
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Agonads MEP Half life.2")
p5.4

