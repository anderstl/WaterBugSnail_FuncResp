
#load packages needed for analysis
library(dplyr)
library(ggplot2)
library(car)
library(cowplot)
library(emdbook)
library(bbmle)
library(deSolve)

#------------------------------------------------------------------------------
# type II FR
#------------------------------------------------------------------------------

eaten.bolker = function(N0, b, h, q=0, Tt, P) {
  a = b*N0^q
  Neaten.est = N0 - lambertW(a*h*N0*exp(-a*(P*Tt-h*N0)))/(a*h)
  return(Neaten.est)
}

nll.bolker = function(Neaten, N0, b, h, q=0, Tt, P){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.bolker(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P)
  nll = -1*sum(dbinom(x = Neaten,
                      size = N0, 
                      prob = y/N0,
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# type III FR
#------------------------------------------------------------------------------

eaten.hassell = function(N0, b, h, Tt, P){
  p = -(b*h*N0^2 + b*P*Tt*N0 + 1) / (b*h*N0)
  q = P*Tt*N0/h
  Neaten.est = -(p/2) - sqrt( (p/2)^2 - q )
  return(Neaten.est)
}

nll.hassell = function(Neaten, N0, b, h, Tt, P){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.hassell(N0=N0, b=b, h=h, Tt=Tt, P=P)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# generalized FR
#------------------------------------------------------------------------------

eq.ode.general = function(t, x, parms){
  with(as.list(parms),{
    dN = -b * x[1]^(1+q) / (1 + b * h * x[1]^(1+q)) * P
    return(list(c(dN)))
  })
}

eaten.ode.general = function(N0, b, h, q, Tt, P, steps=100){
  Neaten.est = vector()
  for(i.eaten in 1:length(N0)){
    Neaten.est[i.eaten] = N0[i.eaten] - lsoda(y = N0[i.eaten],
                                              times = seq(0,Tt[i.eaten],length=steps),
                                              func = eq.ode.general,
                                              parms = c(b=b, h=h, q=q, P=P[i.eaten])
    )[length(seq(0,Tt[i.eaten],length=steps)),2]
  }
  return(Neaten.est)
}

nll.ode.general = function(Neaten, N0, b, h, q, Tt, P, steps=100){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.ode.general(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}


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

#make histograms of snail size by treatment
snail_size<-merge(snail_size,snail_total[,c("Tank","SnailDensity","Complexity")],by.x="Container",by.y="Tank")

#create average snail sizes for each treatment and block
snail_avg<-snail_size%>%
  group_by(Block,Container)%>%
  dplyr::summarise(Mean_Width=mean(Width_mm,na.rm=T),Mean_Aperture=mean(Aperture_mm,na.rm=T),
                   sd_Width=sd(Width_mm,na.rm=T),sd_Aperture=sd(Aperture_mm,na.rm=T))

png("Results/FigS4.png",res=500,height=5,width=7,units="in")
snail_stdev<-ggplot(snail_total,aes(as.factor(SnailDensity),sd_Aperture,fill=Complexity))+
  geom_boxplot()+
  labs(x="Density",y="Aperture Standard Deviation")+
  theme_classic()+
  scale_fill_manual(values=c("#56B4E9","#E69F00","#009E73"),labels=c("High","Low","No"),name="Density")
snail_stdev
dev.off()

#merge all the data sets together
snail_total<-merge(snail_total,bug_size,by.x=c("Block","Tank","UniqueID"),by.y=c("Block","Container","UniqueID"),all=T)
snail_total<-merge(snail_total,snail_avg,by.x=c("Block","Tank"),by.y=c("Block","Container"),all=T)

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
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
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

#polynomial regression
typeI.no<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity,
              data=snail_total,family=binomial,subset=c(Complexity=="none"))
typeII.no<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2),
               data=snail_total,family=binomial,subset=c(Complexity=="none"))
typeIII.no<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2)+I(SnailDensity^3),
                data=snail_total,family=binomial,subset=c(Complexity=="none"))
AIC(typeI.no,typeII.no,typeIII.no)

typeI.low<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity,
               data=snail_total,family=binomial,subset=c(Complexity=="low"))
typeII.low<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2),
                data=snail_total,family=binomial,subset=c(Complexity=="low"))
typeIII.low<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2)+I(SnailDensity^3),
                 data=snail_total,family=binomial,subset=c(Complexity=="low"))
AIC(typeI.low,typeII.low,typeIII.low)

typeI.hi<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity,
              data=snail_total,family=binomial,subset=c(Complexity=="high"))
typeII.hi<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2),
               data=snail_total,family=binomial,subset=c(Complexity=="high"))
typeIII.hi<-glm(cbind(NumberKilled,SnailDensity-NumberKilled)~SnailDensity+I(SnailDensity^2)+I(SnailDensity^3),
                data=snail_total,family=binomial,subset=c(Complexity=="high"))
AIC(typeI.hi,typeII.hi,typeIII.hi)

Tab1<-as.data.frame(round(rbind(summary(typeII.no)$coefficients,
          summary(typeII.low)$coefficients,
          summary(typeII.hi)$coefficients),3))
Tab1$Treatment<-c("None","","","Low","","","High","","")
Tab1$Parameter<-c("Intercept","Linear","Quadratic")
write.csv(Tab1,"Results/NewTable1.csv")

#get maximum likelihood estimates for type 2 and type 3 models
no.t2 <- mle2(nll.bolker,
                 start = list(b = 0.3, h = 1.5),
                 data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="none"],
                             Neaten = snail_total$NumberKilled[snail_total$Complexity=="none"],
                             Tt = 4,
                             P = 1))
no.t3 <- mle2(nll.hassell,
                 start = list(b = 0.5, h =1.5),
                 data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="none"],
                             Neaten = snail_total$NumberKilled[snail_total$Complexity=="none"],
                             Tt = 4,
                             P = 1))
# fit.gen = mle2(minuslogl = nll.ode.general,
#                start = list(b = 1,
#                             h = 1/25,
#                             q = 0),
#                data = list(Neaten = snail_total$NumberKilled[snail_total$Complexity=="none"],
#                            N0 = snail_total$SnailDensity[snail_total$Complexity=="none"],
#                            P = rep(1,nrow(snail_total[snail_total$Complexity=="none",])),
#                            Tt = rep(1,nrow(snail_total[snail_total$Complexity=="none",])))
# )

no.ci<-confint(no.t3,method="quad")


lo.t2 <- mle2(nll.bolker,
              start = list(b = 0.3, h = 1.5),
              data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="low"],
                          Neaten = snail_total$NumberKilled[snail_total$Complexity=="low"],
                          Tt = 4,
                          P = 1))
lo.t3 <- mle2(nll.hassell,
              start = list(b = 0.5, h =1.5),
              data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="low"],
                          Neaten = snail_total$NumberKilled[snail_total$Complexity=="low"],
                          Tt = 4,
                          P = 1))
# lo.g = mle2(minuslogl = nll.ode.general,
#                     start = list(b = 1,
#                                  h = 1/25,
#                                  q = 0),
#                     data = list(Neaten = snail_total$NumberKilled[snail_total$Complexity=="low"],
#                                 N0 = snail_total$SnailDensity[snail_total$Complexity=="low"],
#                                 P = rep(1,nrow(snail_total[snail_total$Complexity=="low",])),
#                                 Tt = rep(1,nrow(snail_total[snail_total$Complexity=="low",])))
#)
low.ci<-confint(lo.t3,method="quad")

hi.t2 <- mle2(nll.bolker,
               start = list(b = 0.3, h = 0.3),
               data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="high"],
                           Neaten = snail_total$NumberKilled[snail_total$Complexity=="high"],
                           Tt = 4,
                           P = 1))
hi.t3 <- mle2(nll.hassell,
               start = list(b = 0.25, h =0.05),
               data = list(N0 = snail_total$SnailDensity[snail_total$Complexity=="high"],
                           Neaten = snail_total$NumberKilled[snail_total$Complexity=="high"],
                           Tt = 4,
                           P = 1))
# #hi.g = mle2(minuslogl = nll.ode.general,
#             start = list(b = 1,
#                          h = 1/25,
#                          q = 0),
#             data = list(Neaten = snail_total$NumberKilled[snail_total$Complexity=="high"],
#                         N0 = snail_total$SnailDensity[snail_total$Complexity=="high"],
#                         P = rep(1,nrow(snail_total[snail_total$Complexity=="high",])),
#                         Tt = rep(1,nrow(snail_total[snail_total$Complexity=="high",])))
#)
hi.ci<-confint(hi.t3,method="quad")

citab<-rbind(no.ci,low.ci,hi.ci)
estimate<-c(coef(no.t3),coef(lo.t3),coef(hi.t3))
Tab3<-as.data.frame(round(cbind(estimate,citab),3))
Tab3$Treatment<-c("None","None","Low","Low","High","High")
Tab3$Parameter<-c(rep(c("b","h"),3))
Tab3
write.csv(Tab3,"Results/NewTable3.csv")

NewFig2<-ggplot(Tab3,aes(Treatment,estimate,color=Parameter,shape=Parameter))+
  geom_point(position=position_dodge(width=0.5))+
  geom_linerange(aes(ymin = `2.5 %`, ymax = `97.5 %`),position=position_dodge(width=0.5))+
  theme_cowplot()+
  labs(y="Estimate")+
  scale_x_discrete(limits=c("None","Low","High"))+
  scale_color_manual(values=c("black","gray"))+
  theme(legend.position=c(0.1,0.9))
png("Results/NewFig2.png",res=500,height=3.5,width=3.5,units="in")
NewFig2
dev.off()

#create vector of new densities
x   <- seq(0,max(snail_total$SnailDensity),length=50)

Tab2<-as.data.frame(rbind(AICc(no.t2,no.t3),AICc(lo.t2,lo.t3),AICc(hi.t2,hi.t3)))
Tab2$AICc<-round(Tab2$AICc,3)
Tab2$dAICc<-c(round(Tab2$AICc[1]-Tab2$AICc[2],3),0,round(Tab2$AICc[3]-Tab2$AICc[4],3),0,round(Tab2$AICc[5]-Tab2$AICc[6],3),0)
Tab2$Treatment<-c("No","","Low","","High","")
Tab2$Model<-c(rep(c("Type II","Type III"),3))
Tab2
write.csv(Tab2,"Results/NewTable2.csv")

#make plot
no.y2  <- eaten.bolker(x,coef(no.t2)[[1]],coef(no.t2)[[2]],0,4,1)
no.y3  <- eaten.hassell(x,coef(no.t3)[[1]],coef(no.t3)[[2]],4,1)
 no.y3upp  <- eaten.hassell(x,Tab3$`97.5 %`[[1]],Tab3$`97.5 %`[[2]],4,1)
 no.y3low  <- eaten.hassell(x,Tab3$`2.5 %`[[1]],Tab3$`2.5 %`[[2]],4,1)
no.pl<-ggplot()+
  geom_point(aes(snail_total$SnailDensity[snail_total$Complexity=="none"],
                  snail_total$NumberKilled[snail_total$Complexity=="none"]))+
  #geom_line(aes(x,no.y2),linetype=2)+
  theme_cowplot()+
  labs(x="",y="Number Killed")+
  scale_y_continuous(breaks=seq(0,16,4),limits=c(0,16))+
  scale_x_continuous(breaks=seq(0,16,4))
  
lo.y2  <- eaten.bolker(x,coef(lo.t2)[[1]],coef(lo.t2)[[2]],0,4,1)
lo.y3  <- eaten.hassell(x,coef(lo.t3)[[1]],coef(lo.t3)[[2]],4,1)
lo.y3upp  <- eaten.hassell(x,Tab3$`97.5 %`[[3]],Tab3$`97.5 %`[[4]],4,1)
lo.y3low  <- eaten.hassell(x,Tab3$`2.5 %`[[3]],Tab3$`2.5 %`[[4]],4,1)
lo.pl<-ggplot()+
  geom_point(aes(snail_total$SnailDensity[snail_total$Complexity=="low"],
                 snail_total$NumberKilled[snail_total$Complexity=="low"]))+
  #geom_line(aes(x,lo.y2),linetype=2)+
  geom_line(aes(x,lo.y3))+
  theme_cowplot()+
  labs(x="",y="Number Killed")+
  scale_y_continuous(breaks=seq(0,16,4),limits=c(0,16))+
  scale_x_continuous(breaks=seq(0,16,4))

hi.y2  <- eaten.bolker(x,coef(hi.t2)[[1]],coef(hi.t2)[[2]],0,4,1)
hi.y3  <- eaten.hassell(x,coef(hi.t3)[[1]],coef(hi.t3)[[2]],4,1)
hi.pl<-ggplot()+
  geom_point(aes(snail_total$SnailDensity[snail_total$Complexity=="high"],
                 snail_total$NumberKilled[snail_total$Complexity=="high"]))+
  #geom_line(aes(x,hi.y2),linetype=2)+
  geom_line(aes(x,hi.y3))+
  theme_cowplot()+
  labs(x=expression("Initial Number of"~italic(Helisoma)),y="Number Killed")+
  scale_y_continuous(breaks=seq(0,16,4),limits=c(0,16))+
  scale_x_continuous(breaks=seq(0,16,4))

png("Results/NewFig1.png",res=600,height=9,width=3.5,units="in")
plot_grid(no.pl,lo.pl,hi.pl,labels=c("A)","B)","C)"),align="hv",ncol=1)
dev.off()

