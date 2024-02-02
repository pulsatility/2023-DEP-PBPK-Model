#Exposure Parameters
exposure_duration = 6                     # h 
CI_Air <- 0.25	                          # DEP concentration in inhaled air, ug/L 
oral_dose = 0                             # ug

#Physiolocal Parameters
BW =75                                    # kg

Fgut = 0.0171
Fliver = 0.0257
Fplasma = 0.0424
Fgonads = 0.0005
Ffat	= 0.2142
Fskin	= 0.0371 
Frestbody = 0.91- Fliver- Fplasma- Fgut- Fgonads- Ffat- Fskin

QP=0*5*60			                            # ventilation rate (L/h) 
QC = 15*BW^0.74                           # (L/h) cardiac blood output

FQgut = 0.178                             # based on <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4844798/> (Table 2)
FQliver = 0.053                           # fractional hepatic liver blood flow
FQgonads = 0.0005                    	    # Fractional gonads blood flow for male
FQskin = 0.058                            # Fractional skin blood flow
FQfat = 0.052                             # Fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin

Vgut= Fgut*BW
Vliver= Fliver*BW
Vplasma= Fplasma*BW
Vgonads= Fgonads*BW
Vfat= Ffat*BW
Vskin= Fskin*BW
Vrestbody= Frestbody*BW

Skin_area =1.91                           # total skin area
FXskin =0.82                              # Fraction of skin exposed, 0.82 (hood Weschler)

kgut_abs <-0  # assuming half of DEP is absorbed in 20 mins

CLGint =991   
CLLint =1642 

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
  
  Qgut= FQgut* QC,
  Qliver= FQliver* QC,
  Qgonads= FQgonads* QC,
  Qskin= FQskin* QC,
  Qfat= FQfat* QC,
  Qrestbody= FQrestbody* QC,
  
  # ave_creatinine_prod_rate = 20,                                 # mg/kg BW/day
  # creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000), # g/h
  
  Skin_area_expo = Skin_area*FXskin,
  kskin_absorption = 500,                                          # Skin absorption rate constant, length/h
  
  #Chemical-specific Parameters
  P_gut_plasma_DEP		=5.92,
  P_liver_plasma_DEP		=5.95,
  P_gonads_plasma_DEP	=3.26,
  P_skin_plasma_DEP		= 7.29,
  P_fat_plasma_DEP		=64.42,
  P_restbody_plasma_DEP	=1.66, 
  P_plasma_air_DEP		= 5e7,
  
  P_gut_plasma_MEP		=1.88,	                                      # This is a parameter we added to account for MEP gut blood flow.
  P_liver_plasma_MEP		=1.89,
  P_gonads_plasma_MEP	=1.16,
  P_skin_plasma_MEP		=2.25,
  P_fat_plasma_MEP		=17.42,
  P_restbody_plasma_MEP	=0.75,
  
  fu_DEP			= 0.30,                                               # fu_DEP is the fraction unbound constant for DEP in plasma
  fu_MEP			=0.598,                                               # fu_MEP is the fraction unbound constant for MEP in plasma
  
  kurine_MEP		=0.25,                                              # 1/h
  
  CLLint =CLLint,
  CLGint=CLGint,
  
  k_glucuronidation		= 1000                                        # k_glucuronidation is rate constant of MEP glucuronidation in the liver
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
  Aglucuronidated_MEP =0,
  Agut_lumen = oral_dose
)

DEPMEP <- function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         Cplasma = Aplasma/Vplasma;	  
         Cplasma_MEP = Aplasma_MEP/Vplasma;
         
         Cgut = Agut/Vgut;
         
         RAMG = CLGint * (Cgut*fu_DEP/P_gut_plasma_DEP);
         
         dAgut = Qgut*Cplasma - Qgut*Cgut/P_gut_plasma_DEP - RAMG + kgut_abs*Agut_lumen
         
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
         RGlucuronidation = k_glucuronidation*Cliver_MEP*(fu_MEP/P_liver_plasma_MEP);
         
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
         
         dAglucuronidated_MEP = RGlucuronidation - kurine_MEP*Aglucuronidated_MEP;
         
         dAgut_lumen = -kgut_abs*Agut_lumen
         
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
           dAglucuronidated_MEP,
           dAgut_lumen),dAurine_MEP=dAurine_MEP,RAML=RAML, RAMG=RAMG);
       })
}

times <- seq(0, exposure_duration, by = 0.01)


##---------------------- kurine_MEP----------------------##

#1 fold
##### during exposure
parameters.1<-parameters


parameters.1["kurine_MEP"]<-parameters["kurine_MEP"]*1 

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters


parameters.1["kurine_MEP"]<-parameters["kurine_MEP"]*0.1 

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters


parameters.1["kurine_MEP"]<-parameters["kurine_MEP"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                                     group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
          yvar=out.final_Aplasma$Aplasma_log,
          title="Aplasma",
          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                          yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                          title="Aurine_MEP",
                          ylim=c(0,3.5))
p.urine_MEP



##---------------------- P_skin_plasma_DEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters


parameters.1["P_skin_plasma_DEP"]<-parameters["P_skin_plasma_DEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters


parameters.1["P_skin_plasma_DEP"]<-parameters["P_skin_plasma_DEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters


parameters.1["P_skin_plasma_DEP"]<-parameters["P_skin_plasma_DEP"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP




##---------------------- P_restbody_plasma_DEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters


parameters.1["P_restbody_plasma_DEP"]<-parameters["P_restbody_plasma_DEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters


parameters.1["P_restbody_plasma_DEP"]<-parameters["P_restbody_plasma_DEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_restbody_plasma_DEP"]<-parameters["P_restbody_plasma_DEP"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP



##---------------------- P_fat_plasma_DEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_fat_plasma_DEP"]<-parameters["P_fat_plasma_DEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_fat_plasma_DEP"]<-parameters["P_fat_plasma_DEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_fat_plasma_DEP"]<-parameters["P_fat_plasma_DEP"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP



##---------------------- P_gonads_plasma_DEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_gonads_plasma_DEP"]<-parameters["P_gonads_plasma_DEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_gonads_plasma_DEP"]<-parameters["P_gonads_plasma_DEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["P_gonads_plasma_DEP"]<-parameters["P_gonads_plasma_DEP"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP


##---------------------- k_glucuronidation----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["k_glucuronidation"]<-parameters["k_glucuronidation"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["k_glucuronidation"]<-parameters["k_glucuronidation"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["k_glucuronidation"]<-parameters["k_glucuronidation"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP


##---------------------- CLLint----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLLint"]<-parameters["CLLint"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLLint"]<-parameters["CLLint"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLLint"]<-parameters["CLLint"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP


##---------------------- CLGint----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLGint"]<-parameters["CLGint"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLGint"]<-parameters["CLGint"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["CLGint"]<-parameters["CLGint"]*10

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP




##---------------------- fu_MEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_MEP"]<-parameters["fu_MEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_MEP"]<-parameters["fu_MEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_MEP"]<-parameters["fu_MEP"]*(1/0.598)

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP




##---------------------- fu_DEP----------------------##
#1 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_DEP"]<-parameters["fu_DEP"]*1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_1 <- out.final

cre_adjusted_Auriune_MEP_1 <-cre_adjusted_Auriune_MEP

#0.1 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_DEP"]<-parameters["fu_DEP"]*0.1

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_0.1 <- out.final
cre_adjusted_Auriune_MEP_0.1 <-cre_adjusted_Auriune_MEP

#10 fold
##### during exposure
parameters.1<-parameters

parameters.1["fu_DEP"]<-parameters["fu_DEP"]*(1/0.3)

out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

##### after exposure
# change CI value
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0)

out.dataframe <- as.data.frame(out)
state.time.2<-as.numeric(tail(out.dataframe,n=1))[-1]

state.time.2.vector <- c(  Agut  =state.time.2[1],
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
                           Aglucuronidated_MEP =state.time.2[19],
                           Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

out.final_10 <- out.final
cre_adjusted_Auriune_MEP_10 <-cre_adjusted_Auriune_MEP

## create dataset
out.final_data <- as.data.frame(rbind(out.final_0.1,out.final_1,out.final_10))
out.final_data$group <-as.factor(c(rep(c(0.1,1,10),each=nrow(out.final_0.1))))

cre_adjusted_Auriune_MEP_data <-as.data.frame(c(cre_adjusted_Auriune_MEP_0.1,
                                                cre_adjusted_Auriune_MEP_1,
                                                cre_adjusted_Auriune_MEP_10)) 

cre_adjusted_Auriune_MEP_data$group <-as.factor(c(rep(c(0.1,1,10),
                                                      each=length(cre_adjusted_Auriune_MEP_0.1))))

colnames(cre_adjusted_Auriune_MEP_data) <- c("cre_adjusted_Auriune_MEP","group")
cre_adjusted_Auriune_MEP_data$time <- out.final_data$time

##---------------------- plot ----------------------##
#Aplasma_DEP
out.final_Aplasma <- out.final_data %>% select(time,Aplasma,group)
out.final_Aplasma$Aplasma_log <- log10(out.final_Aplasma$Aplasma/Vplasma)
out.final_Aplasma$group_l <- factor(out.final_Aplasma$group,levels = c(0.1,10,1))

plot_func <- function(dat,yvar,title,ylim){
  p <- ggplot(data=dat, aes(x =time,y =yvar,
                            group=group_l,colour=group))+
    geom_line(size=1)+xlab("Time (h)") +
    ylab("Log Concentration (ug/L)")+
    theme_classic()+
    theme(axis.title.y= element_text(size = 16,colour="black"),
          axis.title.x = element_text(size = 16,colour="black"),
          axis.text.x = element_text(size = 16,colour="black"),
          axis.text.y = element_text(size = 16,colour="black"),
          title=element_text(size=10))+
    ggtitle(title)+
    ylim(ylim)
  return(p)
}
p.plasma_DEP <- plot_func(dat=out.final_Aplasma,
                          yvar=out.final_Aplasma$Aplasma_log,
                          title="Aplasma",
                          ylim=c(-3,0.5))
p.plasma_DEP 


#Agonads_DEP
out.final_Agonads_DEP <- out.final_data %>% select(time,Agonads,group)
out.final_Agonads_DEP$Agonads_DEP_log <- log10(out.final_Agonads_DEP$Agonads/Vgonads)
out.final_Agonads_DEP$group_l <- factor(out.final_Agonads_DEP$group,levels = c(0.1,10,1))

p.gonads_DEP <- plot_func(dat=out.final_Agonads_DEP,
                          yvar=out.final_Agonads_DEP$Agonads_DEP_log,
                          title="Agonads_DEP",
                          ylim=c(-2,1.8))
p.gonads_DEP 


#Aplasma_MEP
out.final_Aplasma_MEP <- out.final_data %>% select(time,Aplasma_MEP,group)
out.final_Aplasma_MEP$Aplasma_MEP_log <- log10(out.final_Aplasma_MEP$Aplasma_MEP/Vplasma)
out.final_Aplasma_MEP$group_l <- factor(out.final_Aplasma_MEP$group,levels = c(0.1,10,1))

p.plasma_MEP <- plot_func(dat=out.final_Aplasma_MEP,
                          yvar=out.final_Aplasma_MEP$Aplasma_MEP_log,
                          title="Aplasma_MEP",
                          ylim=c(-3,-0.5))
p.plasma_MEP 

#Agonads_MEP
out.final_Agonads_MEP <- out.final_data %>% select(time,Agonads_MEP,group)
out.final_Agonads_MEP$Agonads_MEP_log <- log10(out.final_Agonads_MEP$Agonads_MEP/Vgonads)
out.final_Agonads_MEP$group_l <- factor(out.final_Agonads_MEP$group,levels = c(0.1,10,1))

p.gonads_MEP <- plot_func(dat=out.final_Agonads_MEP,
                          yvar=out.final_Agonads_MEP$Agonads_MEP_log,
                          title="Agonads_MEP",
                          ylim=c(-3,-0.5))
p.gonads_MEP

#Aurine_MEP
cre_adjusted_Auriune_MEP_data$group_l <- factor(cre_adjusted_Auriune_MEP_data$group,levels = c(0.1,10,1))

p.urine_MEP <- plot_func(dat=cre_adjusted_Auriune_MEP_data,
                         yvar=log10(cre_adjusted_Auriune_MEP_data$cre_adjusted_Auriune_MEP),
                         title="Aurine_MEP",
                         ylim=c(0,3.5))
p.urine_MEP

