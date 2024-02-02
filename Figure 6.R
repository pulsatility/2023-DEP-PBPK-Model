library(ggplot2)   
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggjoy)


### Densities plot of prior and posterior parameters uncertainty distribution  

load("MCMC.DEP.RData")
Human.DEP.MCMC<- MCMC.DEP[[1]]

## Sampling from posterior parameters to generate the posterior distributions 

## Human posteior distributions
M.Human.DEP.pos  <- exp(Human.DEP.MCMC$pars) %>% apply(2,mean)    # Get the mean of each column. 1 indicates row; 2 indicates column. Each column has 50000 values.
SD.Human.DEP.pos <- exp(Human.DEP.MCMC$pars) %>% apply(2,sd)

dis.Human.DEP.pos <- matrix(ncol = 50000, nrow = 9)               # Create an empty matrix
rownames(dis.Human.DEP.pos) <- names(M.Human.DEP.pos)[1:9]        # The first 17 parameters, one in a row

logμ=log(M.Human.DEP.pos^2/sqrt(SD.Human.DEP.pos^2+M.Human.DEP.pos^2))
logsigma=sqrt(log(SD.Human.DEP.pos^2/M.Human.DEP.pos^2+1))

for (i in 1:50000){
  for (j in 1:9){
    dis.Human.DEP.pos[j,i] <- rlnorm(1, meanlog  = logμ[j], sdlog  = logsigma[j])# Resample once each parameter for 50000 iterations to get a smooth distribution
  }
}

dis.Human.DEP.pos<-melt(dis.Human.DEP.pos) # reshape the table, each iteration number now in a colunm, make plotting more convenient
names(dis.Human.DEP.pos)=c("Par","Species","Value") # Species is the column name of iteration number
dis.Human.DEP.pos$distribution = c("Posterior") # Change the interation number in the Species column to Human, make plotting more convenient


## Human prior distributions
mean.human.DEP.prior = exp(theta.MCMC.DEP) # M-value for human
human.DEP.prior.sd    = mean.human.DEP.prior*0.5                 # S-value for human

dis.Human.DEP.prior <- matrix(ncol = 50000, nrow = 9)       # Create an empty matrix
rownames(dis.Human.DEP.prior) <- names(mean.human.DEP.prior)[1:9]        # The first 17 parameters, one in a row

logμ_i=log(mean.human.DEP.prior^2/sqrt(human.DEP.prior.sd^2+mean.human.DEP.prior^2))
logsigma_i=sqrt(log(human.DEP.prior.sd^2/mean.human.DEP.prior^2+1))

for (i in 1:50000){
  for (j in 1:9){
    dis.Human.DEP.prior[j,i] <- rlnorm(1,meanlog  = logμ_i[j], sdlog = logsigma_i[j] ) # Resample once each parameter for 50000 iterations to get a smooth distribution
  }
}

dis.Human.DEP.prior<-melt(dis.Human.DEP.prior) # reshape the table, each iteration number now in a colunm, make plotting more convenient
names(dis.Human.DEP.prior)=c("Par","Species","Value") # Species is the column name of iteration number
dis.Human.DEP.prior$distribution = c("Prior") # Change the interation number in the Species column to Human, make plotting more convenient
dis.Human.DEP <- rbind.data.frame (dis.Human.DEP.prior,dis.Human.DEP.pos)
dis.Human.DEP$log.value <- log(dis.Human.DEP$Value)
dis.Human.DEP$ggjoy = c("A")
dis.Human.DEP$Par = factor(dis.Human.DEP$Par, levels = names(M.Human.DEP.pos)[1:9])


p <- ggplot (dis.Human.DEP, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = distribution)) + # fill, fill one color for each species
  geom_joy (scale = 8, size = 0.25, rel_min_height = 0.01, alpha = 0.4) + # over size of the entire ggjoy plot
  facet_wrap ( ~ Par, nrow = 3, ncol= 3,scale="free_x")+# arrange the plot by parameter name, free_x means that the x scale is not fixed, so automatically updated.
  theme_joy()+ # change the theme of the plot p (incorporate p into the theme)
  theme (
    plot.background         = element_blank(),
    text                    = element_text (family = "Times",face="bold",size = 18),
    panel.background        = element_rect (fill   = "#f7f7f7"),
    panel.grid.major.x      = element_blank(), 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    axis.ticks.x            = element_line (size   = 1.5, colour = "black"), 
    axis.ticks.y            = element_blank(), 
    axis.text.y             = element_blank(),
    strip.text              = element_text (size   = 18),
    legend.title            = element_text (size   = 18, face="bold"),
    legend.justification    = "none",
    legend.position         = "none",
    legend.text             = element_text (size = 18,face="bold")) + 
  labs (x = "", y = "") 

p


