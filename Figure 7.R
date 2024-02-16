
library(FME)         # Package for MCMC simulation and model fitting

##------------------------------------------------ Load data sets------------------------------------------------ ##
load("MCMC.DEP.RData")

MC.H.1.DEP = as.mcmc (MCMC.DEP[[1]]$pars)     # first  chain
MC.H.2.DEP = as.mcmc (MCMC.DEP[[2]]$pars)     # second chain
MC.H.3.DEP = as.mcmc (MCMC.DEP[[3]]$pars)     # third  chain

par(mfrow=c(1,1))
MC.H.1.DEP_dataframe <- as.data.frame(MC.H.1.DEP)
MC.H.2.DEP_dataframe <- as.data.frame(MC.H.2.DEP)
MC.H.3.DEP_dataframe <- as.data.frame(MC.H.3.DEP)


cor_1 <- round(cor(MC.H.1.DEP[,"P_skin_plasma_DEP"][5000:10000],MC.H.1.DEP[,"P_restbody_plasma_DEP"][5000:10000]),2)
cor_2 <- round(cor(MC.H.1.DEP[,"P_skin_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kurine_MEP"][5000:10000]),2)
cor_3 <- round(cor(MC.H.1.DEP[,"P_skin_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000]),2)
cor_4 <- round(cor(MC.H.1.DEP[,"P_restbody_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kurine_MEP"][5000:10000]),2)
cor_5 <- round(cor(MC.H.1.DEP[,"P_restbody_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000]),2)
cor_6 <- round(cor(MC.H.1.DEP[,"kurine_MEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000]),2)

cor.test(MC.H.1.DEP[,"P_skin_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kurine_MEP"][5000:10000])
cor.test(MC.H.1.DEP[,"P_skin_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000])
cor.test(MC.H.1.DEP[,"P_restbody_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kurine_MEP"][5000:10000])
cor.test(MC.H.1.DEP[,"P_restbody_plasma_DEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000])
cor.test(MC.H.1.DEP[,"kurine_MEP"][5000:10000],MC.H.1.DEP[,"kskin_absorption"][5000:10000])


n_seq <- 10

##------------------------------------------------ correlation plot------------------------------------------------ ##

## P_restbody_plasma_DEP
plot(c(MC.H.1.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)]), 
       c(MC.H.1.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
         MC.H.2.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
         MC.H.3.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="P_skin_plasma_DEP",ylab="P_restbody_plasma_DEP",col="dark blue")
text(3, 1.8, paste0("cor="," ",cor_1))
text(3, 1.5, "p< .001")


## P_skin_plasma_DEP_kurine_MEP
plot(c(MC.H.1.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)]), 
     c(MC.H.1.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="P_skin_plasma_DEP",ylab="kurine_MEP",col="dark blue")
text(3, -0.75, paste0("cor="," ",cor_2))
text(3, -0.85, "p< .001")
dev.off()


## P_skin_plasma_DEP_kskin_absorption
plot(c(MC.H.1.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$P_skin_plasma_DEP[seq(5000,10000,n_seq)]), 
     c(MC.H.1.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="P_skin_plasma_DEP",ylab="kskin_absorption",col="dark blue")
text(3, 6.1, paste0("cor="," ",cor_3))
text(3, 6.05, "p< .001")
dev.off()


## P_restbody_plasma_DEP_kurine_MEP
plot(c(MC.H.1.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)]), 
     c(MC.H.1.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="P_restbody_plasma_DEP",ylab="kurine_MEP",col="dark blue")
text(1, -1, paste0("cor="," ",cor_4))
text(1, -0.8, "p< .001")
dev.off()

## P_restbody_plasma_DEP_kskin_absorption
plot(c(MC.H.1.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$P_restbody_plasma_DEP[seq(5000,10000,n_seq)]), 
     c(MC.H.1.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="P_restbody_plasma_DEP",ylab="kskin_absorption",col="dark blue")
text(1, 6.1, paste0("cor="," ",cor_5))
text(1, 6.05, "p= 0.039")
dev.off()

## kurine_MEP_kskin_absorption
plot(c(MC.H.1.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kurine_MEP[seq(5000,10000,n_seq)]), 
     c(MC.H.1.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.2.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)],
       MC.H.3.DEP_dataframe$kskin_absorption[seq(5000,10000,n_seq)]),
     pch=20,cex=0.5,xlab="kurine_MEP",ylab="kskin_absorption",col="dark blue")
text(-0.9, 6.1, paste0("cor="," ",cor_6))
text(-0.9, 6.05, "p< .001")
dev.off()


