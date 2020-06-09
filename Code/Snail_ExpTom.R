####BELOSTOMA FUNCTIONAL RESPONSE ACROSS VARYING HABITAT COMPLEXITY####
###EXPERIMENTAL DESIGN AND ANALYSIS 2018###

#load packages needed
library(dplyr)
library(ggplot2)
library(car)
library(cowplot)
library(frair)

#Set directory that is specific to your machine
snail_data <- read.csv("Data/snail_data.csv") 
bug_size <- read.csv("Data/bug_sizes.csv")
snail_size <- read.csv("Data/snail_size.csv")

#order by Block, elapsed hours and tank (in that order); not a critical step but makes it look nicer
snail_data<-snail_data[with(snail_data,order(Block,ElapsedHrs,Tank)),]

#create unique ID for each block and tank combo (because we re-used tank numbers across blocks), using the interaction function
snail_data$UniqueID<-interaction(snail_data$Block,snail_data$Tank) 

#drop control treatments for now (they were tanks 25 and 26 in Block 1, and tanks 37 and 38 in Block 2)
snail_data<-snail_data[!snail_data$UniqueID%in%c(1.25,1.26,2.37,2.38),]

#separate out (or subset) the final survey of number killed, which happened at 96 elapsed hours
snail_total<- snail_data[which(snail_data$ElapsedHrs==96),]

#calculate proportion killed
snail_total$ProportionKilled<-snail_total$NumberKilled/snail_total$SnailDensity

#create unique IDS for block-tank, and then drop the controls
snail_size$UniqueID<-interaction(snail_size$Block,snail_size$Container) 
bug_size$UniqueID<-interaction(bug_size$Block,bug_size$Container) 
snail_size<-snail_size[!snail_size$UniqueID%in%c(1.25,1.26,2.37,2.38),]
bug_size<-bug_size[!bug_size$UniqueID%in%c(1.25,1.26,2.37,2.38),]

#create average snail sizes for each treatment and block
snail_avg<-snail_size%>%
  group_by(Block,Container)%>%
  dplyr::summarise(Mean_Width=mean(Width_mm,na.rm=T),Mean_Aperture=mean(Aperture_mm,na.rm=T))

#merge all the data sets together
snail_total<-merge(snail_total,bug_size,by.x=c("Block","Tank","UniqueID"),by.y=c("Block","Container","UniqueID"),all=T)
snail_total<-merge(snail_total,snail_avg,by.x=c("Block","Tank"),by.y=c("Block","Container"),all=T)
  
#plot proportion and number killed against density to visually assess functional response type
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
figS2a<-ggplot(data = filter(snail_total,Complexity=="none"),aes(SnailDensity,ProportionKilled))+
  geom_point(position=position_jitter(width=0.05),size=2)+
  theme_cowplot()+
  labs(x="",y="Proportion Killed")+
  binomial_smooth(formula= y ~ poly(x, 2),color="black",se=F)+
  scale_x_continuous(breaks=c(2,4,8,16),labels=c(2,4,8,16))
figS2b<-ggplot(data = filter(snail_total,Complexity=="low"),aes(SnailDensity,ProportionKilled))+
  geom_point(position=position_jitter(width=0.05),size=2)+
  theme_cowplot()+
  labs(x="",y="Proportion Killed")+
  binomial_smooth(formula= y ~ poly(x, 2),color="black",se=F)+
  scale_x_continuous(breaks=c(2,4,8,16),labels=c(2,4,8,16))
figS2c<-ggplot(data = filter(snail_total,Complexity=="high"),aes(SnailDensity,ProportionKilled))+
  geom_point(position=position_jitter(width=0.05),size=2)+
  theme_cowplot()+
  labs(x=expression("Initial"~italic(Helisoma)~"Number"),y="Proportion Killed")+
  binomial_smooth(formula= y ~ poly(x, 2),color="black",se=F)+
  scale_x_continuous(breaks=c(2,4,8,16),labels=c(2,4,8,16))
png("Results/FigS2.png",res=600,width=3.5,height=9,units="in")
plot_grid(figS2a,figS2b,figS2c,ncol=1,labels=LETTERS[1:3])
dev.off()

#Analyze whether bug size and snail size varied by treatment or block
snailsize.mod<-lm(Mean_Aperture~Complexity+SnailDensity+Block,data=snail_total)
bugsize.mod<-lm(Width_mm~Complexity+SnailDensity+Block,data=snail_total)
Anova(snailsize.mod,type=3)
Anova(bugsize.mod,type=3)

#plot snail and bug size by treatment
figS3a<-ggplot(snail_total,aes(SnailDensity,Mean_Aperture,color=Complexity,shape=Complexity))+
  geom_point(size=2)+
  labs(x="",y=expression("Mean"~italic(Helisoma)~"Width (mm)"))+
  scale_color_manual(values=c("#56B4E9","#E69F00","#009E73"),labels=c("High","Low","No"))+
  scale_shape_manual(values=c(19,17,15),labels=c("High","Low","No"))+
  geom_smooth(method="lm",se=F)+
  theme_cowplot()+
  scale_x_continuous(breaks=c(2,4,8,16))

figS3b<-ggplot(snail_total,aes(SnailDensity,Width_mm,color=Complexity,shape=Complexity))+
  geom_point(size=2)+
  labs(x=expression("Initial"~italic(Helisoma)~"Number"),y=expression(italic(Belostoma)~"Length (mm)"))+
  scale_color_manual(values=c("#56B4E9","#E69F00","#009E73"),labels=c("High","Low","No"))+
  scale_shape_manual(values=c(19,17,15),labels=c("High","Low","No"))+
  geom_smooth(method="lm",se=F)+
  scale_x_continuous(breaks=c(2,4,8,16))+
  theme_cowplot()
png("Results/FigS3.png",res=600,height=7,width=5,units="in")
plot_grid(figS3a,figS3b,labels=c("A","B"),align="hv",ncol=1)
dev.off()

#Analyse whether snail size, and block bug size impacted survival
snail_total<-within(snail_total,rm(NumberKilled_Count2))
snail_total_noNA<-na.omit(snail_total)

size.mod<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~Complexity+SnailDensity+Mean_Aperture+Width_mm,data=snail_total_noNA,family=quasibinomial)
Anova(size.mod,type=3,test.statistic = "F")
summary(size.mod)
exp(coef(size.mod)) #to get odds ratios

#plot proportion killed by water bug length (called "Width_mm" here) and snail 
fig3b<-ggplot(snail_total,aes(Width_mm,ProportionKilled))+
  geom_point()+
  binomial_smooth(color="black")+
  labs(x=expression(italic(Belostoma)~"Length (mm)"), y=expression("Proportion"~italic(Helisoma)~"Killed"))+
  theme_cowplot()
fig3a<-ggplot(snail_total,aes(Mean_Aperture,ProportionKilled))+
  geom_point()+
  binomial_smooth(color="black")+
  labs(x=expression("Mean"~italic(Helisoma)~"Width (mm)"), y=expression("Proportion"~italic(Helisoma)~"Killed"))+
  theme_cowplot()
  
png("Results/Fig3.png",res=600,height=7,width=3.5,units="in")
plot_grid(fig3a,fig3b,labels=c("A","B"),align="hv",ncol=1)
dev.off()

# determine Type II vs Type III using frair package (from Pritchard et al. 2017, Methods in Ecol and Evol)
ft.no<-frair_test(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="none",])
ft.low<-frair_test(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="low",])
ft.high<-frair_test(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="high",])

#Alternatively, test flexible functional response by Real (1977)
none_flexpnr<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="none",],response = "flexpnr",
                  start=list(b=1,h=0.5,q=0),fixed=list(T=96/24))
low_flexpnr<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="low",],response = "flexpnr",
                  start=list(b=0.05,h=0.5,q=0.01),fixed=list(T=96/24))
high_flexpnr<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="high",],response = "flexpnr",
                  start=list(b=1,h=0.5,q=0),fixed=list(T=96/24))
summary(none_flexpnr$fit)
summary(low_flexpnr$fit)#doesn't converge
summary(high_flexpnr$fit)

#get bootstrap confidence intervals
noneboot_flxpnr<-frair_boot(none_flexpnr,nboot = 5000)
confint(noneboot_flxpnr,citypes="bca")
highboot_flxpnr<-frair_boot(high_flexpnr,nboot = 5000)
confint(highboot_flxpnr,citypes="bca")

#compare two models where the fit converges
frair_compare(none_flexpnr,high_flexpnr)

#try additional type III model for low complexity treatment
low_hasseliii<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="low",],response = "hassIIInr",
                         start=list(b=1,h=1,c=0.2),fixed=list(T=96/24))
summary(low_hasseliii$fit)#still doesn't converge, so likely a type II response

#Test Roger's Type II models to compare treatments
none_typeii<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="none",],response = "rogersII",
                      start=list(a=1,h=1),fixed=list(T=96/24))
low_typeii<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="low",],response = "rogersII",
                      start=list(a=1,h=1),fixed=list(T=96/24))
high_typeii<-frair_fit(NumberKilled~SnailDensity,data=snail_total[snail_total$Complexity=="high",],response = "rogersII",
                      start=list(a=1,h=1),fixed=list(T=96/24))
summary(none_typeii$fit)
summary(low_typeii$fit)
summary(high_typeii$fit)

#do treatment comparisons based on Juliano (2001)
frair_compare(none_typeii,low_typeii)
frair_compare(none_typeii,high_typeii)
frair_compare(low_typeii,high_typeii)

#get bootstrapped confidence intervals
noneboot<-frair_boot(none_typeii,nboot = 5000)
lowboot<-frair_boot(low_typeii,nboot = 5000)
highboot<-frair_boot(high_typeii,nboot = 5000)

#plot best fitting model for each species
png("Results/Fig1.png",res=600,height=9,width=3.5,units="in")
par(mfrow=c(3,1),mar=c(4,5,1,0))
plot(noneboot_flxpnr,xlab="",ylim=c(0,14),ylab="Number Killed",
     cex.axis=2,cex.lab=2,frame=F,las=1,xaxt="n")
abline(h=0)
axis(1,at=c(2,4,8,16),cex.axis=2,pos=0)
drawpoly(noneboot_flxpnr, col=adjustcolor("#009E73",alpha.f=0.5))
lines(noneboot_flxpnr,lwd=3)
points(noneboot_flxpnr,pch=20, col="#009E73",cex=2)
mtext(text="A",font=2,line=3,side=2,at=14,las=1,cex=1.75)
plot(lowboot,xlab="",ylab="Number Killed",
     cex.axis=2,cex.lab=2,frame=F,las=1,ylim=c(0,14),xaxt="n")
axis(1,at=c(2,4,8,16),cex.axis=2,pos=0)
abline(h=0)
drawpoly(lowboot, col=adjustcolor("#E69F00",alpha.f=0.5))
lines(lowboot,lwd=3,lty=2)
points(lowboot,pch=20, col="#E69F00",cex=2)
mtext(text="B",font=2,line=3,side=2,at=14,las=1,cex=1.75)
plot(highboot_flxpnr,xlab=expression("Initial"~italic(Helisoma)~"Number"),ylab="Number Killed",
     frame=F,las=1,ylim=c(0,14),cex.axis=2,cex.lab=2,xaxt="n")
axis(1,at=c(2,4,8,16),cex.axis=2,pos=0)
abline(h=0)
drawpoly(highboot_flxpnr, col=adjustcolor("#56B4E9",alpha.f=0.5))
lines(highboot_flxpnr,lwd=3,lty=3)
points(highboot_flxpnr,pch=20, col="#56B4E9",cex=2)
mtext(text="C",font=2,line=3,side=2,at=14,las=1,cex=1.75)
dev.off()

#plot treatment comparisons (predicted values only)
png("Results/Fig2.png",res=600,height=7,width=3.5,units="in")
par(mfrow=c(2,1),mar=c(4,5,0,0.1))
plot(NA,xlim=c(2,16),ylim=c(0,14),xlab="",ylab="Number Killed",frame=F,las=1,
     xaxt="n",cex.axis=1.25,cex.lab=1.25)
axis(1,at=c(2,4,8,16),cex.axis=1.25,pos=0)
abline(h=0)
drawpoly(noneboot_flxpnr, col=adjustcolor("#009E73",alpha.f=0.5))
lines(noneboot_flxpnr,lwd=3)
mtext(text="A",font=2,line=3,side=2,at=14,las=1,cex=1.75)
#points(noneboot,pch=21, bg="#009E73")
drawpoly(highboot_flxpnr, col=adjustcolor("#56B4E9",alpha.f=0.5))
lines(highboot_flxpnr,lwd=3,lty=3)

plot(NA,xlim=c(2,16),ylim=c(0,14),xlab=expression("Initial"~italic(Helisoma)~"Number"),
     ylab="Number Killed",frame=F,las=1,xaxt="n",cex.axis=1.25,cex.lab=1.25)
axis(1,at=c(2,4,8,16),cex.axis=1.25,pos=0)
abline(h=0)
drawpoly(noneboot, col=adjustcolor("#009E73",alpha.f=0.5))
lines(noneboot,lwd=3)
#points(noneboot,pch=21, bg="#009E73")
drawpoly(lowboot, col=adjustcolor("#E69F00",alpha.f=0.5))
lines(lowboot,lwd=3,lty=2)
#points(lowboot,pch=21, bg="#E69F00")
drawpoly(highboot, col=adjustcolor("#56B4E9",alpha.f=0.5))
lines(highboot,lwd=3,lty=3)
#points(highboot,pch=21, bg="#56B4E9")
mtext(text="B",font=2,line=3,side=2,at=14,las=1,cex=1.75)
dev.off()