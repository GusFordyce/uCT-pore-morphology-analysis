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

allpores <- read.csv("allporedata.csv",header=T)
str(allpores)
allpores$Timepoint <- as.factor(allpores$Timepoint)

#########################-----------------------Pore Size Distribution (Volume parameter in dataset)----------------##############
##Exploratory plotting
xyplot(Volume ~ Treatment.x | Timepoint, data=allpores, jitter.x=T)
  ##Very skewed 
  ##Try comparing kernel density estimates

h <- sd(allpores$Volume)*(4/3/length(allpores$Volume))^(1/5)
h
  ##Computing Silverman's rule of thumb for bandwidth estimation

sm.density.compare(allpores$Volume[which(allpores$Timepoint==1)],
                   allpores$Treatment.x[which(allpores$Timepoint==1)],
                   xlim=c(0,max(allpores$Volume)),h=0.006956163,model="equal")
  ##Fairly incomprehensible
  ##Lets try make this prettier and clearer

##GGplot theme
clear<- theme(axis.text.x=element_text(colour="black",size=18),
                    axis.text.y=element_text(colour="black",size=18),
                    axis.title.x=element_text(colour="black",size=18),
                    axis.title.y=element_text(colour="black",size=18),
                    axis.ticks=element_line(colour="black"),
                    legend.text=element_text(size=18),
                    legend.title=element_text(size=18),
                    legend.position="right",
                    panel.grid.minor=element_line(colour="transparent"),
                    panel.grid.major=element_line(colour="transparent"),
                    panel.background=element_rect(fill="white",colour="black"),
                    plot.background=element_blank()) 


#####Density estimates for Time Point 1
ConDens <- sm.density(allpores$Volume[which(allpores$Treatment.x == "C" & allpores$Timepoint==1)],display="none",h=0.006956163)
HeatDens <- sm.density(allpores$Volume[which(allpores$Treatment.x == "H" & allpores$Timepoint==1)],display="none",h=0.006956163)
Condf<-as.data.frame(cbind(ConDens[[1]],ConDens[[2]],ConDens[[7]],ConDens[[8]]))
Heatdf <-as.data.frame(cbind(HeatDens[[1]],HeatDens[[2]],HeatDens[[7]],HeatDens[[8]]))
Condf <- plyr::rename(Condf,c("V1" = "EvalPointsCont","V2" = "ProbsDensCont",
                              "V3" = "CIUpperCont","V4" = "CILowerCont"))
Heatdf <- plyr::rename(Heatdf,c("V1" = "EvalPointsHeat","V2" = "ProbsDensHeat",
                              "V3" = "CIUpperHeat","V4" = "CILowerHeat"))
plotdf <- cbind(Condf,Heatdf)
str(plotdf)


#####Density estimates for Time Point 0
ConDens0 <- sm.density(allpores$Volume[which(allpores$Treatment.x == "C" & allpores$Timepoint==0)],display="none",h=0.006956163)
HeatDens0 <- sm.density(allpores$Volume[which(allpores$Treatment.x == "H" & allpores$Timepoint==0)],display="none",h=0.006956163)
Condf0<-as.data.frame(cbind(ConDens0[[1]],ConDens0[[2]],ConDens0[[7]],ConDens0[[8]]))
Heatdf0 <-as.data.frame(cbind(HeatDens0[[1]],HeatDens0[[2]],HeatDens0[[7]],HeatDens0[[8]]))
Condf0 <- plyr::rename(Condf0,c("V1" = "EvalPointsCont","V2" = "ProbsDensCont",
                              "V3" = "CIUpperCont","V4" = "CILowerCont"))
Heatdf0 <- plyr::rename(Heatdf0,c("V1" = "EvalPointsHeat","V2" = "ProbsDensHeat",
                              "V3" = "CIUpperHeat","V4" = "CILowerHeat"))

plotdf0 <- cbind(Condf0,Heatdf0)

##Plotting Code
####FOR TIME POINT 1
ggplot(plotdf)+
     geom_line(aes(x=EvalPointsCont,y=ProbsDensCont),col='black')+
     geom_ribbon(aes(x=EvalPointsCont,ymin=CILowerCont,ymax=CIUpperCont),
                 colour="black",alpha=.2,linetype=0)+
     geom_line(aes(x=EvalPointsHeat,y=ProbsDensHeat),col='red')+
     geom_ribbon(aes(x=EvalPointsHeat,ymin=CILowerHeat,ymax=CIUpperHeat),
                 fill='red',alpha=.2,linetype=0)+
     xlim(-0.005,max(plotdf$EvalPointsCont))+
     xlab(bquote("Poresize" ~ (mm^3)))+
     ylab("Probability Density f(x)")+
     clear

ggplot(plotdf)+
     geom_line(aes(x=EvalPointsCont,y=ProbsDensCont),col='black')+
     geom_ribbon(aes(x=EvalPointsCont,ymin=CILowerCont,ymax=CIUpperCont),
                 colour="black",alpha=.2,linetype=0)+
     geom_line(aes(x=EvalPointsHeat,y=ProbsDensHeat),col='red')+
     geom_ribbon(aes(x=EvalPointsHeat,ymin=CILowerHeat,ymax=CIUpperHeat),
                 fill='red',alpha=.2,linetype=0)+
     xlim(-0.005,0.3)+
     xlab(bquote("Poresize" ~ (mm^3)))+
     ylab("Probability Density f(x)")+
     clear


####FOR TIME POINT 0
ggplot(plotdf0)+
     geom_line(aes(x=EvalPointsCont,y=ProbsDensCont),col='black')+
     geom_ribbon(aes(x=EvalPointsCont,ymin=CILowerCont,ymax=CIUpperCont),
                 colour="black",alpha=.2,linetype=0)+
     geom_line(aes(x=EvalPointsHeat,y=ProbsDensHeat),col='red')+
     geom_ribbon(aes(x=EvalPointsHeat,ymin=CILowerHeat,ymax=CIUpperHeat),
                 fill='red',alpha=.2,linetype=0)+
     xlim(-0.005,max(plotdf$EvalPointsCont))+
     xlab(bquote("Poresize" ~ (mm^3)))+
     ylab("Probability Density f(x)")+
     theme_classic()

ggplot(plotdf0)+
     geom_line(aes(x=EvalPointsCont,y=ProbsDensCont),col='black')+
     geom_ribbon(aes(x=EvalPointsCont,ymin=CILowerCont,ymax=CIUpperCont),
                 colour="black",alpha=.2,linetype=0)+
     geom_line(aes(x=EvalPointsHeat,y=ProbsDensHeat),col='red')+
     geom_ribbon(aes(x=EvalPointsHeat,ymin=CILowerHeat,ymax=CIUpperHeat),
                 fill='red',alpha=.2,linetype=0)+
     xlim(-0.005,0.3)+
     xlab(bquote("Poresize" ~ (mm^3)))+
     ylab("Probability Density f(x)")+
     clear

  ##It is better but nonetheless the extreme skewness of the data towards smaller pores really makes for an ugly graph
  ##And we're not able to dive into where the differences (if any) lie in much detail

summdat<-ddply(allpores,c("Treatment.x","Timepoint"),summarise,mean=mean(Volume),
                                                      sd=sd(Volume),se=sd(Volume)/length(Volume))
  ##Surprisingly looks like the mean volume for H at TP1 would be lower

  ##Unclear on how to model highly skewed data such as this. Poisson and negative binomial aren't designed for
  ##non-integer data which pore volume more certainly is. If I were to implement negative binomial regression:

  ##Using a Cullen-Frey Graph to suggest possible distributions

descdist(allpores$Volume,discrete=F,boot=1000)
  ##Based on this, the candidate distributions are lognormal, gamma and weibull

weib.fit<-fitdist(allpores$Volume,"weibull")
gamma.fit<-fitdist(allpores$Volume,"gamma")
lognorm.fit<-fitdist(allpores$Volume,"lnorm")

  ##Comparative Q-Q plot
qqcomp(list(gamma.fit, weib.fit, lognorm.fit),
       legendtext=c("gamma","weibull", "lnorm") )
  ##Lognorm make v difficult to interpret.
  ##QQ plot compares the theoretical quantiles of a particular distribution (x-axis), to the empirical quantiles
  ##of your data. This function also plots a y=x line. The closer each line sits to the y=x line, the more likely
  ##that your data is of that distribution

qqcomp(list(gamma.fit,weib.fit),
       legendtext=c("gamma","weibull") )
  ##Based on this, gamma is the best distribution for modelling this data, then weibull then lognormal
  ##But still skews away from y=x

plot(lognorm.fit)
plot(weib.fit)
plot(gamma.fit)
  ##Based on these, we can see that gamma is better at describing the tails of the distribution but doesn't do as well
  ##as lognormal at capturing the middle
  ##Other lines of evidence discount lognormal, which means in a gamma vs. weibull, gamma is better

cdfcomp(list(gamma.fit,weib.fit),
       legendtext=c("gamma","weibull") )

ppcomp(list(gamma.fit,weib.fit),
       legendtext=c("gamma","weibull"))

   ##Neither is conclusive - I would lean towards Gamma because it is more likely to capture to tail of this 
   ##heavily skewed data

   ##It seems there is no straightforward way to implement weibull regression with a random effect (sample.ID)
   ##That doesn't explicitly relate to survival analysis
   ##I am selecting Gamma regression to model this relationship - indeed it is designed for heavily skewed, non-negative data


gammatest <- glmer(Volume ~ Treatment.x * as.factor(Timepoint)+(1|Sample.ID), allpores,family=Gamma(link="log"))
summary(gammatest)
   ##Nothing stands out as being horrible
   ##The median deviance residual in the model is near zero which is good

##Check linearity in relationship to predictor on log scale and distribution of variance
res <- residuals(gammatest, type="deviance")
ypred <- predict(gammatest)
plot(ypred, res)
  ##The density of points makes it hard to approxiamte the mean deviance to zero
  ##But it does look like the distribution is still very skewed
fitdistr(allpores$Volume,"gamma")
  ##Indeed when you generate shape and rate, the shape is < 0.25 so extremely skewed
  ##which matches up with the bootstrapped values in the Cullen-Frey graph which show that kurtosis and skewness
  ##are not one-to-one

dat <- as.data.frame(cbind(res,ypred))

dat %>%
mutate(bin=cut_width(ypred,width=0.01,boundary=0)) %>%
ggplot(aes(x=bin,y=res))+geom_boxplot()

  ##The variance is not drastically different based on the IQR
  ##It does increase with larger pore sizes (less negative on the log scale)but this is to be expected given how
  ##few big pores there are relative to small pores


##Confirming consistency of coefficient of variation 
##i.e. checking near-constant shape parameters
boxplot(Volume ~ Timepoint * Treatment.x, data=allpores)
gamma.fitC<-fitdist(allpores$Volume[which(allpores$Timepoint == 1 & allpores$Treatment.x == "C")],"gamma")
1/sqrt(gamma.fitC$estimate[1])
gamma.fitH<-fitdist(allpores$Volume[which(allpores$Timepoint == 1 & allpores$Treatment.x == "H")],"gamma")
1/sqrt(gamma.fitH$estimate[1])
gamma.fitC0<-fitdist(allpores$Volume[which(allpores$Timepoint == 0 & allpores$Treatment.x == "C")],"gamma")
1/sqrt(gamma.fitC0$estimate[1])
gamma.fitH0<-fitdist(allpores$Volume[which(allpores$Timepoint == 0 & allpores$Treatment.x == "H")],"gamma")
1/sqrt(gamma.fitH0$estimate[1])

   ##So the coefficient of variation is consistent across the factorial levels
   ##and the variance in linear predictor vs deviance residuals has no clear or severe pattern
   ##beyond what you might expect for extremely skewed data. 

#Using least squared means constrasts for post-hoc
lsm <- lsmeans(gammatest, ~ Treatment.x | as.factor(Timepoint))
contrast(lsm, interaction  = "pairwise", adjust="tukey")


  ##This is a very extreme distribution that is difficult to model well, as evidenced by the very low shape parameter
  ##The shape parameter is basically the same which means the skewness is not too different
  ## < 1 means the distribution is effectively asymptotic (true for y, not really for x but data structure is so similar)
  ##The rate parameter is about 50% larger for heat treatment, suggesting the the rate of 'decay'
  ##i.e. the steepness of the slope, is greater


##plotting means
poreplot<-ddply(allpores, c("Treatment.x","Timepoint"),summarise, mean=mean(Volume),se=sd(Volume)/length(Volume))
poreplot$Treatment <- ifelse(poreplot$Treatment.x == "C", "Control", "Treatment") 
PorePlotTP0 <- ggplot(poreplot[which(poreplot$Timepoint =="0"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","red"))+
               xlab("")+
               ylab("Mean Pore Volume" ~ (mm^3))+
               ylim(0,0.04)+
               clear


PorePlotTP1 <- ggplot(poreplot[which(poreplot$Timepoint =="1"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","red"))+
               xlab("")+
               ylab("Mean Pore Volume" ~ (mm^3))+
               ylim(0,0.04)+
               clear

##Combined in illustrator

#########################-----------------------Pore Size Distribution (Volume parameter in dataset)----------------##############
##Exploratory plotting
str(allpores)
bwplot(Ratio~Treatment.x|Timepoint,data=allpores,ylim=c(0,0.2))
  ##From this, it looks like maybe the ratio goes up (based on median and IQR) which would indicate a 'roughening' of pore surfaces?

plot(density(allpores$Ratio))
polygon(density(allpores$Ratio),col='red',border='blue')
  ##Based on this probability density plot, I'm going to use gamma again for parity with pore size distribution analysis
allpores$Ratio[which(is.na(allpores$Ratio))]

ratio.mod <- glmer(Ratio ~ Treatment.x * as.factor(Timepoint) + (1|Sample.ID), allpores,family=Gamma(link="log"))
summary(ratio.mod)
  ##Median deviance is v close to zero

boxplot(Ratio ~ Timepoint * Treatment.x, data=allpores)
gamma.fitC<-fitdist(allpores$Ratio[which(allpores$Timepoint == 1 & allpores$Treatment.x == "C")],"gamma")
cv.c1<-1/sqrt(gamma.fitC$estimate[1])
gamma.fitH<-fitdist(allpores$Ratio[which(allpores$Timepoint == 1 & allpores$Treatment.x == "H")],"gamma")
cv.h1<-1/sqrt(gamma.fitH$estimate[1])
gamma.fitC0<-fitdist(allpores$Ratio[which(allpores$Timepoint == 0 & allpores$Treatment.x == "C")],"gamma")
cv.c0<-1/sqrt(gamma.fitC0$estimate[1])
gamma.fitH0<-fitdist(allpores$Ratio[which(allpores$Timepoint == 0 & allpores$Treatment.x == "H")],"gamma")
cv.h0<-1/sqrt(gamma.fitH0$estimate[1])
dat<-as.data.frame(rbind(cv.c0,cv.h0,cv.c1,cv.h1))
mean(dat$shape)
sd(dat$shape)
  ##Coefficient of variation looks pretty consistent and with low sd relative to values
  ##The shape parameters also look way less extreme and skewed.

res <- residuals(ratio.mod, type="deviance")
ypred <- predict(ratio.mod)
plot(ypred, res)
  ##very nice and consistent

dat <- as.data.frame(cbind(res,ypred))

dat %>%
mutate(bin=cut_width(ypred,width=0.01,boundary=0)) %>%
ggplot(aes(x=bin,y=res))+geom_boxplot()
  ##Very nice

lsm <- lsmeans(ratio.mod, ~ Treatment.x | Timepoint)
contrast(lsm,interaction="pairwise",adjust="tukey")
  ##No sig change in pore 'shape'



##Generating sphere SA:Vol data for a plot
#Volume of a sphere is (4/3)*pi*(radius^3)
#S.A. of a sphere is 4 * pi * (radius^2)
#How long is this fake data set?
nlevels(as.factor(allporesSA$Volume))
##41108 different values for volume
nlevels(as.factor(allporesSA$SurfaceArea))
##127522 for S.A
##Try with smaller number first

##Start at zero, where do we end? lets go with 1.5 because of TP1 graph
spherevol <- seq(from=0,to=1.49,by=(1.5/41108))
length(spherevol) #Good enough - 40834

#Now back-calculate what the radii would be for these theoretical spheres
sphereradius <- (3*(spherevol/(4*pi)))^(1/3)
spheresa <- 4*pi*(sphereradius^2)

spheredat<-as.data.frame(cbind(spherevol,spheresa))
class(spheredat)
str(spheredat)
head(spheredat)

cubevol <- spherevol
cubelength <- spherevol ^ (1/3)
cubesa <- 6 * (cubelength ^ 2)
cubedat <- as.data.frame(cbind(cubevol, cubesa))


ggplot()+
      geom_smooth(data=spheredat,
                  aes(x=spherevol,y=spheresa),colour="royalblue3",lty=2,alpha=.2)+
      geom_smooth(data=cubedat,
                  aes(x=cubevol,y=cubesa),colour="purple",lty=2,alpha=.2)+
      xlab("Volume" ~ (mm^3))+
      ylab("Surface Area" ~ (mm^2))+
      clear



##For timepoint end
TP1all<-ggplot()+
      geom_smooth(data=allpores[-which(allpores$Timepoint==0),],
                  aes(x=Volume,y=SurfaceArea,colour=Treatment.x,fill=Treatment.x),alpha=.2)+
      scale_fill_manual(values=c("black","red"))+
      scale_colour_manual(values=c("black","red"))+
      geom_smooth(data=spheredat,
                  aes(x=spherevol,y=spheresa),colour="royalblue3",lty=2,alpha=.2)+
      ylim(0,11)+
      xlim(0,2)+
      xlab("Volume" ~ (mm^3))+
      ylab("Surface Area" ~ (mm^2))+
      clear

##Barely changes
##For timepoint zero
TP0all<-ggplot()+
      geom_smooth(data=allpores[-which(allpores$Timepoint==1),],
                  aes(x=Volume,y=SurfaceArea,colour=Treatment.x,fill=Treatment.x),alpha=.2)+
      scale_fill_manual(values=c("black","red"))+
      scale_colour_manual(values=c("black","red"))+
      geom_smooth(data=spheredat,
                  aes(x=spherevol,y=spheresa),colour="royalblue3",lty=2,alpha=.2)+
      ylim(0,11)+
      xlim(0,2)+
      xlab("Volume" ~ (mm^3))+
      ylab("Surface Area" ~ (mm^2))+
      clear

grid.arrange(TP0all,TP1all,ncol=2)


