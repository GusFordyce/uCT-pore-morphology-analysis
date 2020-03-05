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
library(gdata)
library(pastecs)
library(MASS)
library(sm)
library(fitdistrplus)

##Read in data frame
data <- read.csv ('pdam.mp.analysis.csv',header=T)
str(data)
head(data)

##Each of the columns X1 - X0 refers to a microporosity bin. 
##Values are expressed as number of voxels in each bin.
##We need to redefine Time Point as factor and transform porosity variables into proportions bound
##by [0,1] for the purposes of beta regression


data$Time.Point <- as.factor(data$Time.Point)
data$MPprop <- data$Macro/100
data$MiPprop <- data$Total.Micro/100


#################-----------Macro Porosity-----------#################
MPmod <- betareg (MPprop ~ Treatment * Time.Point , data=data, link="logit")
summary(MPmod)
  #The estimate of the precision parameter in beta regression is akin to 'variance explained' 
  #in linear models. A positive estimate for the precision parameter indicates that 
  #the inclusion of the explanatory variables (i.e. covariates) increases the precision
  ##of the model and reduces variance in predicted values

plot(MPmod)
  #We expect standardised weighted residuals for beta regression to be centred around zero
  #For a continous predictor, we would test for linearity between the dependent and independent variables
  #on the log scale (due to logit link). But as these predictors are categorical, this is not possible.
  #Instead, look at the residuals (deviance) vs linear predictor plot produced by plot(model)
  #Look for any severe deviations from being centred around zero in each factorial group

  ##Residuals, cook's and leverage diagnostics all skewed by the outlier, sample 555
  ##Attempt without for reporting purposes

  ##Pairwise contrasts for the above model us lsmeans from the emmeans package
lsmMP <- lsmeans(MPmod, ~Treatment | Time.Point)
MPcont <- contrast(lsmMP, interaction = "pairwise",adjust="Tukey")
MPcont
  ##No significant difference, as expected given the results of Leggat et al. (2019)
  ##When interpreting these contrasts and back-transforming the estimates (exp(est)), remember these
  ##are calculated off the estimated marginal means in the model. Meaning the estimated % change given
  ##by back transforming will not precisely match the change in the raw data


  ##Due to the outlier (sample ID 555), we re-ran the analysis while excluding this.

MPmodexcl <- betareg (MPprop ~ Treatment * Time.Point , data=data[-which(data$Sample.ID == 555),], link="logit")
summary(MPmodexcl)
  ##Estimates for variables and precision increase, as does pseudo r-squared

lsmMPexcl <- lsmeans(MPmodexcl, ~Treatment | Time.Point)
MPcontexcl <- contrast(lsmMPexcl, interaction = "pairwise",adjust="Tukey")
MPcontexcl
  ##But this does not change the outcome of the pairwise contrasts

##########----PLOTTING----
##Plots are not in the manuscript but can be used to visualise the data
##GGplot theme
Macro_theme<- theme(axis.text.x=element_text(colour="black",size=14),
                    axis.text.y=element_text(colour="black",size=14),
                    axis.title=element_text(colour="black",size=17),
                    axis.ticks=element_line(colour="black"),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=12,face="bold"),
                    legend.position="none",
                    panel.grid.minor=element_line(colour="transparent"),
                    panel.grid.major=element_line(colour="transparent"),
                    panel.background=element_rect(fill="white",colour="black"),
                    plot.background=element_rect(fill="transparent",colour=NA))  

##Creating separate datasets for means with/without outlier
MPexcl555 <- ddply(data[which(data$Sample.ID != 555),],
                      c("Treatment","Time.Point"),
                      summarize,
                      mean=mean(Macro),
                      se=sd(Macro)/6)


MacroPlotTP1 <- ddply(data,
                      c("Treatment","Time.Point"),
                      summarize,
                      mean=mean(Macro),
                      se=sd(Macro)/6)

comb<-rbind(MacroPlotTP1, MPexcl555)
comb <- comb[-c(5,7,8),]
comb$Treatment <- ifelse(comb$mean < 10, "Control (excl. outlier)", ifelse(comb$mean < 12, "Control", ifelse(comb$mean > 20, "Control", "Treatment")))



PDMacroPlotTP0 <- ggplot(comb[which(comb$Time.Point =="0"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","forestgreen"))+
               xlab("")+
               ylab("Macroporosity (%)")+
               ylim(0,25)+
               Macro_theme


PDMacroPlotTP1 <- ggplot(comb[which(comb$Time.Point =="1"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","royalblue3","forestgreen"))+
               xlab("")+
               ylab("Macroporosity (%)")+
               ylim(0,25)+
               Macro_theme

grid.arrange(PDMacroPlotTP0,PDMacroPlotTP1,ncol=2)

###-------------------------------------------------------------------------------------------------###
#################-----------TotalMicro Porosity-----------#################
MiPmod <- betareg (MiPprop ~ Treatment * Time.Point , data=data, link="logit")
summary(MiPmod)
plot(MiPmod, type="deviance")
  #Nothing worrying aside from outlier (sample 555)

lsmMiP <-lsmeans(MiPmod, ~ Treatment | Time.Point)
MiPcont <- contrast(lsmMiP, interaction = "pairwise", adjust = 'Tukey')
MiPcont
  ##In agreement with the original study, we see a significant difference in the means at TP1
  
  ##Repeating analysis but without sample 555
MiPmodexcl <- betareg (MiPprop ~ Treatment * Time.Point , data=data[-which(data$Sample.ID == 555),], link="logit")
summary(MiPmodexcl)
  #A relatively minor increase in phi and pseudo-r squared compared to that in macro porosity

lsmMiPexcl <-lsmeans(MiPmodexcl, ~ Treatment | Time.Point)
MiPcontexcl <- contrast(lsmMiPexcl, interaction = "pairwise", adjust = 'Tukey')
MiPcontexcl
  ##Very minor changes to the estimate and standard error
  ##and the broad results do not change.
 

##########----PLOTTING----
##Plots are not in the manuscript but can be used to visualise the data
##GGplot theme
Micro_theme<- theme(axis.text.x=element_text(colour="black",size=12),
                    axis.text.y=element_text(colour="black",size=14),
                    axis.title.x=element_blank(),
                    axis.title.y=element_text(colour="black",size=17),
                    axis.ticks=element_line(colour="black"),
                    legend.text=element_text(size=12),
                    legend.title=element_blank(),
                    legend.position="none",
                    panel.grid.minor=element_line(colour="transparent"),
                    panel.grid.major=element_line(colour="transparent"),
                    panel.background=element_rect(fill="white",colour="black"),
                    plot.background=element_rect(fill="transparent",colour=NA)) 

#Creating separate data sets
MiPexcl555 <- ddply(data[which(data$Sample.ID != 555),],
                      c("Treatment","Time.Point"),
                      summarize,
                      mean=mean(Total.Micro),
                      se=sd(Total.Micro)/6)


MicroPlot <- ddply(data,
                      c("Treatment","Time.Point"),
                      summarize,
                      mean=mean(Total.Micro),
                      se=sd(Total.Micro)/6)

comb<-rbind(MicroPlot, MiPexcl555)
comb <- comb[-c(5,7,8),]
comb$Treatment <- ifelse(comb$mean < 3, "Control (excl. outlier)", ifelse(comb$mean < 5, "Control", "Treatment"))

#plots

PDMicroPlotTP0 <- ggplot(comb[which(comb$Time.Point =="0"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","forestgreen"))+
               xlab("")+
               ylab("Microporosity (%)")+
               ylim(0,10)+
               Micro_theme


PDMicroPlotTP1 <- ggplot(comb[which(comb$Time.Point =="1"),])+
               geom_bar(aes(x=Treatment,y=mean,fill=Treatment),position=position_dodge(),stat="identity",width=.5,col="black",alpha=.8)+
               geom_errorbar(aes(ymax=mean+se, ymin=mean-se,x=Treatment),width=.2,size=1)+
               scale_fill_manual(values=c("black","royalblue3","forestgreen"))+
               xlab("")+
               ylab("Microporosity (%)")+
               ylim(0,10)+
               Micro_theme

grid.arrange(PDMicroPlotTP0,PDMicroPlotTP1,ncol=2)


###-------------------------------------------------------------------------------------------------###
#################-----------Surface Area to Volume Ratio-----------#################
##SA:Vol is not beta distributed so attempting to fit a linear model first

SAmodLM <- lm(sa.vol ~ Treatment * Time.Point, data = data)
summary(SAmodLM)
plot(SAmodLM)
 ##Two possible candidates for outliers 
 ##Normality questionable as a result (see labeled points)
 ##Those same two outliers are likely responsible for potential deviations from homoscedasiticity
car::ncvTest(SAmodLM) #Breusch-Pagan test
   ##This suggests that there is no significant heteroscedasticity

hist(resid(SAmodLM),breaks=20)
 #Looks to me like one point in particular is skewing normality distirbution

barplot(cooks.distance(SAmodLM))
   ##555 is an outlier in this instance as well
   ##5 and 18 are high but below the threshold

shapiro.test(data$sa.vol)
   ##This suggests that the assumption of normality is not violated and the graphical method is being
   ##muddied by the outlier

lsm <- lsmeans(SAmodLM,~Treatment | Time.Point)
contrast(lsm, interaction = "pairwise",adjust="Tukey")
   ##This highlights no sig difference at either time point


   ##And if we re-run it without the outlier?

SAmodLMexcl <-  lm(sa.vol ~ Treatment * Time.Point, data = data[-which(data$Sample.ID ==  555),])
summary(SAmodLMexcl)
plot(SAmodLMexcl))
   ##R-squared has gone up a fair bit 
   ##Graphical diagnostics look better

contrast(lsmeans(SAmodLMexcl,~Treatment|Time.Point),interaction="pairwise",adjust="tukey")
   ##Controls greater than treatment at TP0


###-------------------------------------------------------------------------------------------------###
#################-----------Pore Density-----------#################
sm.density.compare(data$pore_density[which(data$Time.Point ==1)],
                      data$Treatment[which(data$Time.Point ==1)])
colfill<-c(2:(2+length(levels(data$Treatment)))) 
legend(locator(1), levels(data$Treatment), fill=colfill)
  ##Based on this, it looks like heat treatments might have higher density at TP1

sm.density.compare(data$pore_density[which(data$Time.Point ==0)],
                      data$Treatment[which(data$Time.Point ==0)])
colfill<-c(2:(2+length(levels(data$Treatment)))) 
legend(locator(2), levels(data$Treatment), fill=colfill)
  ##And that relationship switched relative to TP0

descdist(data$pore_density,boot=1000)
  ##This suggests that the data approximate to a uniform/linear distribution

dens <- lm(pore_density ~ Treatment*Time.Point, data)
summary(dens)
plot(dens)
  ##Unexpected large tail due to two points (16/7) which are also disrupting the variance distribution
hist(resid(dens),breaks=20)
  ##Graphical deviations from normality seem to be a result of the two large outliers
car::ncvTest(dens)
  ##Homogeneity of variance assumption met
barplot(cooks.distance(dens))
shapiro.test(data$pore_density)
  ##Normality not violated

##Producing stripchart for Figure S1c
stripchart(poredata$pore_density ~ poredata$Treatment*poredata$Time.Point,
           pch=1,col="orange", xlab="Pore Density (pores mm-3)")

influencePlot(dens)
  ##7 and 16 are considerable outliers
  ##Despite tests above, this is clearly degrading model fit (r squared <.3)
  ##both high densities, at TP1 - one in C, one in H

  ##Try log transformation
logdens <- lm(log(pore_density) ~ Treatment*Time.Point, data)
plot(logdens)
influencePlot(logdens)
summary(logdens)
  ##Better but still a problem

lsm <- lsmeans(logdens, ~ Treatment | Time.Point)
contrast(lsm, interaction = "pairwise", adjust="tukey")
  ##No sig diff.

##TimePoint 0
tp0dens <- lm(pore_density ~ Treatment, poredata[which(poredata$Time.Point=="0"),])
plot(tp0dens)
  ##more var and bigger spread in H
plot(density(resid(tp0dens)))
  ##Looks good
summary(tp0dens)
  ##insig.

##Raw data - no outliers
dens.out <- lm(pore_density ~ Treatment*Time.Point, data[-c(7,16),])
plot(dens.out)
plot(density(resid(dens.out)))
summary(dens.out)
  ##Remarkable improvement in model fit
  ##R-squared jumps from .23 to .8
  ##the estimates don't change considerably

AIC(dens)
AIC(dens.out)
  ##Fit is better without the outliers

lsm <- lsmeans(dens.out, ~ Treatment | Time.Point)
contrast(lsm, interaction = "pairwise", adjust="tukey")
  ##So this suggests that at TP1, control is less dense than H
  ##and at TP0, it is more dense

  ##Looked at the tomogram and data
  ##It seems like it actually does just have loads of little pores and the mean pore size is smaller compared to other
  ##16 is not particularly high compared to the entire population but it is relative to that timepoint/treatment subpopulation

##summarised with outliers
ddply(data,c("Treatment","Time.Point"),summarise, mean=mean(pore_density))

##summarised without outliers
ddply(data[-c(7,16),],c("Time.Point","Treatment"),summarise,mean=mean(pore_density),se=(sd(pore_density)/sqrt(5)))


###-------------------------------------------------------------------------------------------------###
#################-----------Microporous Phase Distribution-----------#################
colnames(data)
mpdata <- data[,-c(4:10,113,114)] #removing the unnecessary columns before melting
colnames(mpdata)

test<-melt(mpdata,id.vars=c(1:3),measure.vars=c(4:105))
head(test)
test
test <- plyr::rename(test,c("variable"="MPbin","value"="no.voxels"))
head(test)

 
 ##diff total no. of voxels in each sample will influence comparisons.
 ##normalising by converting to proportion of total
 ##First create new dataset with total sum of voxels in each sample
denom <- ddply(test, ~Sample.ID,summarise,sample_sum=sum(no.voxels))
 ##Then re-merge this back into the original data set
merged <- merge(test,denom,by="Sample.ID")                         
str(merged)

 ##Create "proportion of total" variable
merged$sample_prop <- merged$no.voxels / merged$sample_sum

 ##Exploratory plotting
par(mfrow=c(2,1))
plot(merged$sample_prop[which(merged$Treatment == "H")]~ merged$MPbin[which(merged$Treatment == "H")])
plot(merged$sample_prop[which(merged$Treatment == "C")]~ merged$MPbin[which(merged$Treatment == "C")])
       ##Distribution for controls looks a lot flatter


###-----MODELLING/STAT TESTING
 #First we need to clip to X0 and X1 bins because these represent solid carbonate and open pore in MP analysis.
 #These exist due to defining thresholds, the MP analysis includes the value at threshold
clipped_merged <- merged[-which(merged$MPbin %in% c("X1","X0")),]
clipped_merged <- droplevels(clipped_merged)
str(clipped_merged)


nrow(clipped_merged[which(clipped_merged$sample_prop == 0),])
 ##73 cases - transformation for model needed (as per vignette for Betareg)
clipped_merged$s_prop_transformed <- ((clipped_merged$sample_prop * 5) + 0.5)/6

##Full model - threeway interaction between MPbin, Treatment and Time point
##Because we are interested in the treatment*time point interaction within each bin
PD.mod <- betareg(s_prop_transformed ~ Treatment*MPbin*Time.Point, data = clipped_merged,link='logit')
plot(PD.mod)
summary(PD.mod)
  ##The funnel shape in the deviance residuals suggests that the model is skewed towards small values
  ##But it is still centred which supports the linearity between dependents and non-dependents
  ##Further exploratory plotting

pred <- predict(PD.mod,type="link") #link type relates to linear predictor
resid<- residuals(PD.mod,type="sweighted2")
comb<-as.data.frame(cbind(pred,resid))

comb %>%
mutate(bin=cut_width(pred,width=0.01,boundary=0)) %>%
ggplot(aes(x=bin,y=resid))+geom_boxplot()
  ##zero-centering is good despite a few deviations. Variance looks better but still not great
  ##The big jump is likely a product of different sample sizes at each time point leading to greater SD/variance at TP0
  ##We also see that the zero cluster in the original plot is largely due to a few point that great that point


lsm <- lsmeans(PD.mod, ~Treatment | MPbin + Time.Point)
Int.Cont <- contrast(lsm, interaction = "pairwise",adjust="Tukey")
 ##.325 up to .185 and .715 and .845




############-------PLOTTING CODE-------#############
##First, summarise the data
plot.data <- ddply(clipped_merged, 
                   c("Treatment","MPbin","Time.Point"),
                   summarise,
                   mean = mean(sample_prop),
                   se = sd(sample_prop)/sqrt(6))
plot.data <- arrange(plot.data, Time.Point)


num_mp <- rep(seq(0.995, 0.005, by= -0.01),4)
plot.data$MPphase <- num_mp

dist_theme<- theme(axis.text.x=element_text(colour="black",size=16),
                    axis.text.y=element_text(colour="black",size=16),
                    axis.title.x=element_text(colour="black",size=20),
                    axis.title.y=element_text(colour="black",size=20),
                    axis.ticks=element_line(colour="black"),
                    legend.text=element_text(size=14),
                    legend.title=element_blank(),
                    legend.position="none",
                    panel.grid.minor=element_line(colour="transparent"),
                    panel.grid.major=element_line(colour="transparent"),
                    panel.background=element_rect(fill="white",colour="black"),
                    plot.background=element_rect(fill="transparent",colour=NA)) 

####BOTH LINES TOGETHER - T1, C + H
TP1Comparison <- ggplot()+
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







####BOTH LINES TOGETHER - T0, C + H 
TP0Comparison<-ggplot()+
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





