library(ggplot2)
library(gridExtra)
library(lattice)
library(plyr)
library(dplyr)
library(reshape2)
library(betareg)
library(MCMCglmm)
library(coda)
library(plotrix)
library(mltools)
library(gdata)
library(quantreg)
library(lme4)
library(coin) #median_test
library(lmtest)
library(MASS)
library(emmeans)
library(faraway)
library(boot)
library(survreg)
library(fitdistrplus)
library(sm)
library(DescTools)

data <- read.csv("AllMCS_0.5.csv",head=T)
str(data)

##Exploratory plotting
xyplot(Proportion ~ Bin | Treatment + TimePoint, data[which(data$TimePoint == 1),], type=c("p","smooth"))
  ##Looks a bit like a narrower distribution
  ##But with a higher average/median
 
par(mfrow = c(2,2))
plot(density(data$Proportion[which(data$TimePoint==0 & data$Treatment == "C")]))
plot(density(data$Proportion[which(data$TimePoint==0 & data$Treatment == "H")]))
plot(density(data$Proportion[which(data$TimePoint==1 & data$Treatment == "C")]))
plot(density(data$Proportion[which(data$TimePoint==1 & data$Treatment == "H")]))

  ##Skewed but not as much the pore data. This is probably reflective of the nature of carbonates in general, when scanned at this resolution


##As with microporosity analysis, this data structure requires a beta regression
##We need to apply a transformation to deal with zero values
data$Transformed <- ((data$Proportion * 5) + 0.5) / 6

test <- betareg(Transformed ~ Treatment * as.factor(TimePoint) * as.factor(Bin), data, link = 'logit')
summary(test)
plot(test,type="sweighted2")

pred <- predict(test,type="link") #link type relates to linear predictor
resid<- residuals(test,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

qqnorm(resid)
qqline(resid)
hist(resid,breaks=500)
  ##the distribution of weighted residuals looks fine but with a very big zero peak
  ##Look back at the first xyplot
  ##This may be due to the first two bins, 0.25 and 0.75 as output by WebMango
  ##It automatically creates these labelled bin when you select .5 voxels as the radius 
  ##to bin data by. It allows it to bin voxels into 1.25 voxels radius etc. and so maximises the resolution
  ##of your data. It is relates to how you define a sphere in voxelised space
  ##They are consistently zero across all samples

comb %>%
mutate(bin=cut_width(pred,width=0.1,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
  ##This looks very good in terms of zero-centering and consistency of variance
  ##Left bin is a problem but might also reflect the above nuance

dat.clean<-data[-which(data$Bin %in% c(0.25,0.75)),]
test <- betareg(Transformed ~ Treatment * as.factor(TimePoint) * as.factor(Bin), dat.clean, link = 'logit')
summary(test)

pred <- predict(test,type="link") #link type relates to linear predictor
resid<- residuals(test,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

comb %>%
mutate(bin=cut_width(pred,width=0.1,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
  ##Slightly better but the data inherently has a fair number of zeros
  ##because a few large spheres in a few samples
  ##necessitate that there be bins created with zero values in the others

test2 <- betareg(sqrt(Transformed) ~ Treatment * as.factor(TimePoint) * as.factor(Bin), dat.clean, link = 'logit')
  ##Can't apply loq or reciprocal transformations due to [0,1] bounds
pred <- predict(test2,type="link") #link type relates to linear predictor
resid<- residuals(test2,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

comb %>%
mutate(bin=cut_width(pred,width=0.1,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
  ##And a sqrt transformation makes no difference

  ##I'm going to use the original model to avoid additional transformation

lsm <- lsmeans(test, ~Treatment | TimePoint + Bin)
cont<-contrast(lsm, interaction = 'pairwise', adjust='tukey')
cont
str(cont)
cont <- as.data.frame(cont)
sig.cont <- cont[which(cont$p.value <= 0.05),]
sig.cont

summ<-ddply(data,c("Treatment","TimePoint","Bin"),summarise,mean=mean(Proportion))
summ<-cbind(summ[c(1:62),],summ[c(63:124),])
str(summ)
summ2<-summ[which(summ$Bin %in% c("1.25","2.75","3.25","4.75","5.25","5.75","6.25")),-c(6,7)]
summ2 <- summ2[-which(summ2$TimePoint == "0"),]
summ2
  ##So in terms of a comparison, it looks like MCS in Heat is greater for 2.75 and 3.25
  ##But less for the others
  ##Overall it looks like a growth in the 'middle' sig bins and a reduction in the end.
  ##NOTE - these are all relative proportions


write.matrix(as.data.frame(cont),file="MISContrasts.tsv",sep='\t')

MISsummary<-ddply(data,c("Treatment","TimePoint","Bin"),summarise,mean=mean(Proportion),se=sd(Proportion)/sqrt(6))
MISsummary$se[which(MISsummary$TimePoint == 0)] <- MISsummary$se[which(MISsummary$TimePoint == 0)]*sqrt(6)
MISsummary$se[which(MISsummary$TimePoint == 0)] <- MISsummary$se[which(MISsummary$TimePoint == 0)]/sqrt(3)
MISsummary

write.matrix(MISsummary, file="MISBinValues.tsv", sep='\t')


plot <- ddply(data,c("Treatment","TimePoint","Bin"),summarise,mean=mean(Proportion),se=(sd(Proportion)/sqrt(6)))
plot$Bin <- as.factor(plot$Bin)
GG <-ggplot(plot[which(plot$TimePoint == "1"),])+
   geom_point(aes(x=Bin,y=mean,colour=Treatment))+
   geom_smooth(aes(x=Bin,y=mean,colour=Treatment),span=.5)+
   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, x=Bin,colour=Treatment),width=0.2,size=1,alpha=.4)+
   scale_colour_manual(values=c("black","red"))+
   xlab("Sphere Radius (voxels)")+
   ylab("Relative Proportion")+
   theme_classic()
GG + annotate("text",x=c("1.25","2.75","3.25","4.75","5.25","5.75","6.25"),
                     y=c(0.08,0.075,0.195,0.05,0.115,0.055,0.075),label="*",size=8)

