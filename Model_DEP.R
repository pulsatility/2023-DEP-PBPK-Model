
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

Vplasma= Fplasma*BW                       # as a global variable to calculate plasma concentration

QP=5*60			                              # ventilation rate (L/h) (0 for wearing hood)
QC = 15*BW^0.74                           # (L/h) cardiac blood output

FQgut = 0.178                             # fractional GI blood flow
FQliver = 0.053                           # fractional hepatic artery blood flow
FQgonads = 0.0005                       	# fractional gonads blood flow for male
FQskin = 0.058                            # fractional skin blood flow
FQfat = 0.052                         	  # fractional fat blood flow
FQrestbody = 1- FQliver -FQgut-FQfat-FQgonads-FQskin # fractional RB blood flow

Skin_area =1.91                           # total skin area
FXskin =0.91                              # fraction of skin exposed, 0.82 (Weschler hood), 0.91 (Weschler no hood)

kgut_abs <-0                              # absorption rate constant for oral exposure

ave_creatinine_prod_rate = 20             # mg/kg BW/day
creatinine_prod_rate = ave_creatinine_prod_rate *BW/(24*1000)


# paramter values fed to ODE solver
parameters <- c(
  #Exposure Parameters
  CI = CI_Air,
  oral_dose = oral_dose, 
  
  
  #Calculated Physiolocal Parameters
  Vgut= Fgut*BW,
  Vliver= Fliver*BW,
  Vplasma= Fplasma*BW,
  Vgonads= Fgonads*BW,
  Vfat= Ffat*BW,
  Vskin= Fskin*BW,
  Vrestbody= Frestbody*BW,
  
  Qgut= FQgut*QC,
  Qliver= FQliver*QC,
  Qgonads= FQgonads*QC,
  Qskin= FQskin*QC,
  Qfat= FQfat*QC,
  Qrestbody= FQrestbody*QC,
  
  Skin_area_expo = Skin_area*FXskin,
  
  
  #Chemical-specific Parameters
  kskin_absorption = 500,                                          # Skin absorption rate constant
  
  P_gut_plasma_DEP		=5.92,
  P_liver_plasma_DEP	=5.95,
  P_gonads_plasma_DEP	=3.26,
  P_skin_plasma_DEP		=7.29,
  P_fat_plasma_DEP		=64.42,
  P_restbody_plasma_DEP	=1.66, 
  P_plasma_air_DEP		=5e7,
  
  P_gut_plasma_MEP		=1.88,	                                      
  P_liver_plasma_MEP	=1.89,
  P_gonads_plasma_MEP	=1.16,
  P_skin_plasma_MEP		=2.25,
  P_fat_plasma_MEP		=17.42,
  P_restbody_plasma_MEP	=0.75,
  
  fu_DEP			=0.30,                                               # fu_DEP is the fraction unbound constant for DEP in plasma
  fu_MEP			=0.598,                                              # fu_MEP is the fraction unbound constant for MEP in plasma
  
  
  #Biochemical Parameters
  CLGint =991,  
  CLLint =1642,
  
  k_glucuronidation		=1000,                                       # k_glucuronidation is rate constant of MEP glucuronidation in the liver
  kurine_MEP		=0.25                                              # 1/h
  
)


# State variables
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


# ODEs
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
