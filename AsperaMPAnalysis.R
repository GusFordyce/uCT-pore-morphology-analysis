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


data <- read.csv ('aspera.mp.analysis.csv',header=T)
str(data)
data$Time.Point <- as.factor(data$Time.Point)

###---------------------------------------------------------------------------------------------------#####
######------ Macro Porosity--------#####
data$MPprop <- data$Macro/100
MPmod <- betareg (MPprop ~ Treatment * Time.Point , data=data, link="logit")
summary(MPmod)
plot(MPmod, which = 1:6, type="pearson")
 ##Looks pretty good
 ##The variance is not v constant but it zero centred
 ##inconsistencies reflect sample size differences
 ##and there is no directionality to the pattern
 ##Potential outliers > 0.22

cooks.distance(MPmod)
 ##Only sample 7

###Contrasts 
lsmMacro <- lsmeans(MPmod, ~Treatment | Time.Point)
Int.ContMacro <- contrast(lsmMacro, interaction = "pairwise",adjust="Tukey")
  ##This mirrors the results for Total. Porosity 

###Excluding outlier
MPmod <- betareg (MPprop ~ Treatment * Time.Point , data=data[-which(data$Sample.ID ==  597),], link="logit")
summary(MPmod)
 ##Better model fit and precision 

###Contrasts 
lsmMacro <- lsmeans(MPmod, ~Treatment | Time.Point)
Int.ContMacro <- contrast(lsmMacro, interaction = "pairwise",adjust="Tukey")
Int.ContMacro
  ##The difference is in fact greater without the outlier
  ##But no biological reason to remove it


###---------------------------------------------------------------------------------------------------#####
######------ Total Micro Porosity--------#####
data$MiPprop <- data$Total.Micro/100
MiPmod <- betareg (MiPprop ~ Treatment * Time.Point , data=data, link="logit")
summary(MiPmod)
plot(MiPmod,which=1:6,type="pearson")
cooks.distance(MiPmod)
 ###No outliers, half-normal plot looks good
 ###Zero-centred variance but possible directionality

pred <- predict(MiPmod,type="link") #link type relates to linear predictor
resid<- residuals(MiPmod,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

comb %>%
mutate(bin=cut_width(pred,width=0.01,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
  ###Because of how values are distributed across time, the variance widen with decreasing residual value

  ###I repeated this with a sqrt() transformation but it made negligible difference to the model fit 
  ###And no difference to the outcome

lsmMicro <- lsmeans(MiPmod, ~Treatment | Time.Point)
Int.ContMicro <- contrast(lsmMicro, interaction = "pairwise",adjust="Tukey")
Int.ContMicro
  ##Sig change in MicroP at T1 but not T0. In contrast to results of Leggat et al. 2019

summdat<-ddply(data,c("Treatment","Time.Point"),summarize,mean=mean(MiPprop),se=sd(MiPprop)/sqrt(6))
summdat$se[c(1,3)] <- summdat$se[c(1,3)]*sqrt(6)
summdat$se[c(1,3)] <- summdat$se[c(1,3)]/sqrt(3)
summdat
  ##MiP lower for Heat at TP1 - unexpected!

###---------------------------------------------------------------------------------------------------#####
######------ Whole Sample SA:Ratio--------#######
colnames(data)
colnames(data)[8] <- "SAVol"
str(data$SAVol)

SAmodLM <- lm(SAVol ~ Treatment * Time.Point, data = data)
summary(SAmodLM)
  ##Crap fit
plot(SAmodLM)
hist(resid(SAmodLM),breaks=15)
cooks.distance(SAmodLM)
   ##QQplot/residual histogram looks ok given the relatively small sample size.
   ##Variance seems pretty constant around zero
   ##Apparent outliers are < 0.22

lsmSAvol <- lsmeans(SAmodLM, ~Treatment | Time.Point)
ContSAvol <- contrast(lsmSAvol, interaction = "pairwise",adjust="Tukey")
ContSAvol

sumsavol<-ddply(data,c("Treatment","Time.Point"),summarise,mean=mean(SAVol),se=sd(SAVol)/sqrt(6))
sumsavol$se[c(1,3)] <- sumsavol$se[c(1,3)]*sqrt(6)
sumsavol$se[c(1,3)] <- sumsavol$se[c(1,3)]/sqrt(3)
sumsavol


###-------------------------------------------------------------------------------------------------###
#################-----------Microporous Phase Distribution-----------#################
str(data)
mpdataall <- data[,-c(4:10)]
str(mpdataall)

alltest<-melt(mpdataall,id.vars=c(1:3),measure.vars=c(4:105))
alltest
alltest <- plyr::rename(alltest,c("variable"="MPbin","value"="no.voxels"))
   ##Clip data
alltest <- alltest[-which(alltest$MPbin %in% c("X1","X0")),]
alltest <- droplevels(alltest)
str(alltest)
   ##create proportion var
propdat <- ddply(alltest, ~Sample.ID,summarise,sample_sum=sum(no.voxels))
alltest <- merge(alltest,propdat,by="Sample.ID") 
   #create propoprtion
alltest$prop <- alltest$no.voxels / alltest$sample_sum
   #identfy zeros/need for transformation
nrow(alltest[which(alltest$prop == 0),])
   ##60 cases - transformation needed as per vignette
alltest$prop_trans <- ((alltest$prop * 5) + 0.5)/6
nrow(alltest[which(alltest$prop_trans == 0),])
   ##Now run model with Time.Point
str(alltest)

all.mod <- betareg(prop_trans ~ Treatment*MPbin*Time.Point, data = alltest, link = "logit")
summary(all.mod)
plot(all.mod, which = 1:6)
a <- cooks.distance(all.mod)
max(a)
##The median deviance is v close to zero
##Nothing has a high cook's distance to suggest it is an outlier

pred <- predict(all.mod,type="link") #link type relates to linear predictor
resid<- residuals(all.mod,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

comb %>%
mutate(bin=cut_width(pred,width=0.05,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
   ##Variance inflates a bit at the lower residual values
   ##But this isn't necessarily a distinct pattern
   ##Likely this is because many bins at the tail ends of the distribution have largely zeroes
   ##So any non-zero values inflate variance alot

###Post hoc
lsm <- lsmeans(all.mod, ~Treatment | MPbin + Time.Point)
Int.Cont <- contrast(lsm, interaction = "pairwise",adjust="Tukey")
str(Int.Cont)
adf.sig <- as.data.frame(Int.Cont)
write.matrix(adf.sig,file="AspMicroPContrasts.tsv",sep='\t')

vals<-ddply(alltest,c("Treatment","MPbin","Time.Point"),summarise,mean=mean(prop),se=sd(prop)/sqrt(6))
vals$se[which(vals$Time.Point==0)] <- vals$se[which(vals$Time.Point==0)]*sqrt(6)
vals$se[which(vals$Time.Point==0)] <- vals$se[which(vals$Time.Point==0)]/sqrt(3)
vals
write.matrix(vals,file="AspMicroPValues.tsv",sep='\t')


str(adf.sig)
adf.sig.end <-adf.sig[which(adf.sig$p.value <= .05 & adf.sig$Time.Point == "1"),]
adf.sig.en

############-------PLOTTING CODE-------#############
##First, summarise the data
plot.data <- ddply(alltest, 
                   c("Treatment","MPbin","Time.Point"),
                   summarise,
                   mean = mean(prop),
                   se = sd(prop)/sqrt(6))
plot.data <- arrange(plot.data, Time.Point)

num_mp <- rep(seq(0.995, 0.005, by= -0.01),4)
plot.data$MPphase <- num_mp
str(plot.data)
   ####################DIVERGING TO NEW SCRIPT FOR PLOTTING SINGLE
   ####################SAMPLE FOR THE PAPER, AS AN EXAMPLE

####BOTH LINES TOGETHER - T1, C + H
TP1CompAspera<-ggplot()+
    geom_point(aes(MPphase,mean),col="black",
               data = plot.data[which(plot.data$Treatment == "C" & plot.data$Time.Point == "1"),])+
    geom_errorbar(aes(ymin=mean-se, ymax = mean+se,x = MPphase),
               data=plot.data[which(plot.data$Treatment == "C" & plot.data$Time.Point == "1"),],
               width = 0,size=1,alpha=.4)+
    geom_point(aes(MPphase,mean),col="forestgreen",
               data = plot.data[which(plot.data$Treatment == "H" & plot.data$Time.Point == "1"),])+
    geom_errorbar(aes(ymin=mean-se, ymax = mean+se,x = MPphase),col = "forestgreen",
               data=plot.data[which(plot.data$Treatment == "H" & plot.data$Time.Point == "1"),],
               width = 0,size=1,alpha=.4)+
    xlab("% Microporosity")+
    ylab("Fraction of Microporosity Voxels")+
    ylim(0,0.045)+
dist_theme

ggsave("TP1-A-Asp-Micro-Comp.pdf",plot=TP1CompAspera,width=8, height = 5)



####BOTH LINES TOGETHER - T0, C + H 
TP0CompAspera<-ggplot()+
    geom_point(aes(MPphase,mean),col="black",
               data = plot.data[which(plot.data$Treatment == "C" & plot.data$Time.Point == "0"),])+
    geom_errorbar(aes(ymin=mean-se, ymax = mean+se,x = MPphase),
               data=plot.data[which(plot.data$Treatment == "C" & plot.data$Time.Point == "0"),],
               width = 0,size=1,alpha=.4)+
    geom_point(aes(MPphase,mean),col="forestgreen",
               data = plot.data[which(plot.data$Treatment == "H" & plot.data$Time.Point == "0"),])+
    geom_errorbar(aes(ymin=mean-se, ymax = mean+se,x = MPphase),col = "forestgreen",
               data=plot.data[which(plot.data$Treatment == "H" & plot.data$Time.Point == "0"),],
               width = 0,size=1,alpha=.4)+
    xlab("% Microporosity")+
    ylab("Fraction of Microporosity Voxels")+
    ylim(0,0.045)+
dist_theme

ggsave("TP0-A-Asp-Micro-Comp.pdf",plot=TP0CompAspera,width=8, height = 5)







