
# Krais 2018: 300 ug/m^3 exposed for 3 hours wearing clothes.
# Simulation results of the human DEP model for two different DEP exposure scenarios in the Krais study by using the default parameter values.


library(deSolve)
library(tidyverse)
library(ggplot2)
library(ggthemes)

# Import model code containing parameter and ODEs
source("Model_DEP.R")

exposure_duration=3
times <- seq(0, exposure_duration, by = 0.01)

CI_Air <- 0.3
parameters["CI"]<-CI_Air


##### -------------Figure S2A: To simulate Krais hood -------------- ##### 
QP=0			                              # ventilation rate (L/h) 
FXskin =0.04                               # fraction of skin exposed
parameters["Skin_area_expo"]<-Skin_area*FXskin

# Run simulation - during exposure
parameters.1<-parameters
out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

# Run simulation - after exposure
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0) # set air concentration to 0

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
                          Aglucuronidated_MEP =state.time.2[19],
                          Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind.data.frame(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

time <- out.final[,1]
concentration <- log10(cre_adjusted_Auriune_MEP)
dat <- data.frame(time,concentration)
dat$time <- round(dat$time,2)
dat$group <- "orgin"

#Load Krais 2018 datasets - Hood
krais.data.hood.se<-read.table("Krais hood - mean and se.txt",sep = "\t",col.names = c("time","UrineMEP","se_concentration","se_time"),
                               fill = TRUE,strip.white = TRUE,na.strings = "#NA")
krais.data.hood.se$ymin_conc <- log10(krais.data.hood.se$UrineMEP-krais.data.hood.se$se_concentration)
krais.data.hood.se$ymax_conc <- log10(krais.data.hood.se$UrineMEP+krais.data.hood.se$se_concentration)

krais.data.hood.se$xmin_t <- krais.data.hood.se$time-krais.data.hood.se$se_time
krais.data.hood.se$xmax_t <- krais.data.hood.se$time+krais.data.hood.se$se_time

# Plot Figure S2A - Urine MEP
p.krais.hood <-ggplot(krais.data.hood.se)+
  geom_errorbar(aes(x=time,y=log10(UrineMEP),ymin=ymin_conc, ymax=ymax_conc), 
                width=1)+
  geom_errorbarh(aes(x=time,y=log10(UrineMEP),xmin=xmin_t, xmax=xmax_t))+
  geom_point(data =krais.data.hood.se,aes(x=time, y=log10(UrineMEP)) ,shape=1,size=2)+
  geom_line(data=dat,aes(x=time, y=concentration),color="blue",size=1)+
  theme_classic()+ylab("Log10 Urinary MEP/Creatinine (µg/g)")+xlab("Time (h)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  ylim(c(-0.5,2))+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12,colour="black"),
        axis.text.y = element_text(size = 12,colour="black"))

p.krais.hood


##### -------------Figure S2B: To simulate Krais no hood -------------- ##### 
QP=5*60			                              # ventilation rate (L/h) 
FXskin =0.13                               # fraction of skin exposed
parameters["Skin_area_expo"]<-Skin_area*FXskin

# Run simulation - during exposure
parameters.1<-parameters
out <- lsoda(y = state, times = times, func = DEPMEP, parms = parameters.1, jacfunc = NULL)

# Run simulation - after exposure
parameters.2 <- parameters.1
parameters.2["CI"]<-c(0) # set air concentration to 0

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
                          Aglucuronidated_MEP =state.time.2[19],
                          Agut_lumen=state.time.2[20]                  
)

times.2 <- seq(exposure_duration,60,by = 0.01)
out.2 <- lsoda(y = state.time.2.vector , times = times.2, func = DEPMEP, parms = parameters.2, jacfunc=NULL )
out.final <-rbind.data.frame(out,out.2[-1,])
cre_adjusted_Auriune_MEP <- out.final[,"dAurine_MEP"]/creatinine_prod_rate #ug/g

time <- out.final[,1]
concentration <- log10(cre_adjusted_Auriune_MEP)
dat <- data.frame(time,concentration)
dat$time <- round(dat$time,2)
dat$group <- "orgin"


#Load Krais 2018 datasets - no Hood
krais.data.nohood.se<-read.table("Krais no hood - mean and se.txt",sep = "\t",col.names = c("time","UrineMEP","se_concentration","se_time"),
                                 fill = TRUE,strip.white = TRUE,na.strings = "#NA")
krais.data.nohood.se$ymin_conc <- log10(krais.data.nohood.se$UrineMEP-krais.data.nohood.se$se_concentration)
krais.data.nohood.se$ymax_conc <- log10(krais.data.nohood.se$UrineMEP+krais.data.nohood.se$se_concentration)

krais.data.nohood.se$xmin_t <- krais.data.nohood.se$time-krais.data.nohood.se$se_time
krais.data.nohood.se$xmax_t <- krais.data.nohood.se$time+krais.data.nohood.se$se_time

# Plot Figure S2B - Urine MEP
p.krais.nohood <-ggplot(krais.data.nohood.se)+
  geom_errorbar(aes(x=time,y=log10(UrineMEP),ymin=ymin_conc, ymax=ymax_conc), 
                width=1)+
  geom_errorbarh(aes(x=time,y=log10(UrineMEP),xmin=xmin_t, xmax=xmax_t))+
  geom_point(data =krais.data.nohood.se,aes(x=time, y=log10(UrineMEP)) ,shape=1,size=2)+
  geom_line(data=dat,aes(x=time, y=concentration),color="blue",size=1)+
  theme_classic()+ylab("Log10 Urinary MEP/Creatinine (µg/g)")+xlab("Time (h)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  ylim(c(-0.5,3.5))+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12,colour="black"),
        axis.text.y = element_text(size = 12,colour="black"))

p.krais.nohood
