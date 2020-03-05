library(ggplot2)
library(lattice)
library(plyr)
library(dplyr)
library(reshape2)
library(betareg)
library(MCMCglmm)
library(coda)
library(plotrix)
library(car)
library(sandwich)
library(multcomp)
library(gridExtra)
library(tidyverse)
library(tweenr)
library(gganimate)
library(emmeans)
library(lmtest)
library(usdm) ##For VIF
library(fitdistrplus)
library(Hmisc)
library(mice) ##for multiple imputations
library(VIM) ##Imputation plotting

#Data frame manipulation
PDdata <- read.csv("pdam.mp.analysis.csv",head=T)
PDdata <- PDdata[-which(PDdata$Time.Point == 0),c(1:3,5,6,9)]
PDdata <- arrange(PDdata,Treatment)

AAdata <- read.csv("aspera.mp.analysis.csv",head=T)
AAdata <- AAdata[-which(AAdata$Time.Point == 0),c(1:3,6,7,10,11)]
AAdata <- arrange(AAdata,Treatment)

PD.summary <- read.csv("pdamporesummarystats.csv",head=T)
PD.summary <- PD.summary[-which(PD.summary$Time.Point == 0),c(1:3,6:8)]
PD.summary <- arrange(PD.summary,Treatment)
PDdata<-merge(PDdata,PD.summary,by="Sample.ID")
PDdata <- arrange(PDdata,Treatment.x)
PDdata
harddata <- read.csv("HardnessR.csv",head=T)

##Need to impute hardness data before introducing the other parameters
harddata
harddata[c(21:24),] <- NA
harddata[c(21,22),1]<- "Aspera"
harddata[c(23,24),1]<- "Pdam"
harddata[c(21,23),2]<- "C"
harddata[c(22,24),2]<- "H"

aggr(harddata, col=c('navyblue','yellow'),
                    numbers=TRUE, sortVars=TRUE,
                    labels=names(harddata), cex.axis=.7,
                    gap=3, ylab=c("Missing data","Pattern"))

harddata$Int <- interaction(harddata$Species,harddata$Treatment)
spec <- split(harddata,harddata$Int)
spec <- lapply(seq_along(spec), function(x) as.data.frame(spec[[x]])[, 1:5])
  ##gets rid of 'interaction' column that was used to split
  ##Creates list of four data frames, eahc corresponding to a level of the two-factor interaction
str(spec)


testimp<-lapply(seq_along(spec),function(x) mice(as.data.frame(spec[[x]]), m=20,maxit=50,method='pmm'))
  ##Applies the imputation to each data frame in the spec list
str(testimp)

complete(testimp[[1]],c(1:20))
class(testimp[[1]]$imp$Outer) 

con<-rbind(testimp[[1]],testimp[[2]])
heat<-rbind(testimp[[3]],testimp[[4]])
full<-rbind(con,heat)
complete(full,action="long")
  ##This has printed out, effectively, 20 datasets stacked on top of each other
  ##Each imputation is its own dataset, where every sixth value in each factor combination (i.e. Treatment X Species)
 
full.dat <- complete(full,action="long",include=T) ##Need to include the original data so that R knows which

  ##Pocillopora damicornis
full.pd<-full.dat[which(full.dat$Species =="Pdam"),]   
full.pd$MacroP <- PDdata$Macro
full.pd$MicroP <- PDdata$Total.Micro
full.pd$MeanPore <- PDdata$mean
full.pd$MedPore <- PDdata$median
full.pd$Dens <- PDdata$pore_density

   ##Creating 'proportion' [0,1] variables for betareg
full.pd$InnerProp <- full.pd$Inner / 100
full.pd$PoreProp <- full.pd$Pore / 100
full.pd$OuterProp <- full.pd$Outer / 100

test3<-as.mids(full.pd)
complete(test3,"long")
   ##Useful to check it has worked properly
             

imp.fit <- with(test3, betareg(InnerProp ~ MacroP + MicroP + MedPore + Dens))
summary(pool(imp.fit))

imp.fit2 <- with(test3, betareg(OuterProp ~ MacroP + MicroP + MedPore + Dens))
summary(pool(imp.fit2))

imp.fit3 <- with(test3, betareg(PoreProp ~ MacroP + MicroP + MedPore + Dens))
summary(pool(imp.fit3))
