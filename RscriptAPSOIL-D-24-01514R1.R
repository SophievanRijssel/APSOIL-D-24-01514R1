rm(list=ls())


#####load libraries#####
library(plyr)
library(dplyr)
library("ggplot2")
library(Rmisc)
library(lmerTest)
library(reshape2)
library(car) 
library(phyloseq)
library(readr)
library(ggpattern)
library(metagenomeSeq)
library(metagMisc)
library(vegan)
library(emmeans)
library(multcomp)
library(multcompView)

#####read in data#####

dataEva<-read.csv("xxxxx.csv",header=T,row.names=1,sep=";",text=text,dec=".",stringsAsFactors=FALSE,check.names = FALSE)

#####potential added calculations#####
dataEva$SIR.difference<-dataEva$SIR.end-dataEva$SIR.start
dataEva$net.micr.growth<-dataEva$SIR.end-dataEva$SIR.start
dataEva$POXC.difference<-dataEva$POXC.end-dataEva$POXC.start
dataEva$Organic.Soil.percentage<-factor(dataEva$Organic.Soil.percentage)
dataEva$Substrate<-factor(dataEva$Substrate,levels=c("Cover Crop","Carrot leaves","Alfalfa","Hay","Straw","Control"))
dataEva$replicate<-factor(dataEva$replicate)
dataEva$Pair<-factor(dataEva$Pair)
dataEva$Inoculation<-factor(dataEva$Inoculation)
dataEva$Management<-factor(dataEva$Management)
dataEva$TOTALCO2produced.total<-dataEva$TOTALCO2produced.total/1000000
dataEva$POXC.end.fracC<-dataEva$POXC.end / dataEva$percentC_soil.start / 1000
dataEva$TOTALCO2produced.total.diffb<-dataEva$TOTALCO2produced.total.diffb/1000000
dataEva$Nutrients.start<-dataEva$no3.start+dataEva$nh4.start

#####subsets#####
dataEvasubs<-subset(dataEva,!Substrate=="Control")
dataEvaman<-subset(dataEva,Inoculation=="no")
dataEvabasic<-subset(dataEva,Inoculation=="no"&Substrate=="Control")

#####Table 1 Basic soil info#####
dETable1<-subset(dataEva,Substrate=="Control")
dETable1<-subset(dETable1,Inoculation=="no")
dETable1<-subset(dETable1,select=c(Field,replicate,Management,Year,pH,P.Olsen.start,SIR.start,SOM.start,POXC.start,CN.ratio_soil.start,no3.start,nh4.start,Bacterial.richness,Bacterial.Shannon,Bacterial.PLFAs,Fungal.richness,Fungal.Shannon,Fungal.PLFA))
write.csv(dETable1,"Table1basicsoilinfo.csv")

#####test characteristics diffs #####
summary(aov(percentC.substrate~Substrate,data=dataEva))
summary(aov(percentC.substrate~Substrate,data=dataEva))
summary(aov(percentC.substrate~Substrate,data=dataEva))
summary(aov(percentC.substrate~Substrate,data=dataEva))
summary(aov(percentC.substrate~Substrate,data=dataEva))
summary(aov(percentC.substrate~Substrate,data=dataEva))


#####chemical properties crop residues#####
dataCR<-read.csv("P:\\Chapter 5 Inoculation mesocosm experiment\\Data\\Carbon and Nitrogen\\CNcropresiduesforR.csv",header=T,row.names=1,sep=";",text=text,dec=".",stringsAsFactors=FALSE,check.names = FALSE)
summary(aov(percentC~Sample,data=dataCR))
summary(aov(percentN~Sample,data=dataCR))
summary(aov(CN.ratio~Sample,data=dataCR))
dataCR$CN.ratio
TukeyHSD(aov(percentC~Sample,data=dataCR))
TukeyHSD(aov(percentN~Sample,data=dataCR))
TukeyHSD(aov(CN.ratio~Sample,data=dataCR))
aggregate(percentC~Sample,data=dataCR,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(percentN~Sample,data=dataCR,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(CN.ratio~Sample,data=dataCR,function(x) c(mean = mean(x), sd = sd(x)))


#####Real data#####
aggregate(percentC.substrate~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(percentN.substrate~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(total.C.input.substrate~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(total.N.input.substrate~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Mass.loss~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Mass.loss~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate((TOTALCO2produced.total/1000000)~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(total.gram.C.respired~Substrate,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
ggplot(data=dataEva,aes(x=Management,y=NO3H4uptake.resin, fill=Inoculation))+geom_boxplot()+facet_wrap("Substrate")
#summary absolute changes in C and N
#N > to milligram because too many zero;s
aggregate(total.gram.C.input.soil~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate((total.gram.N.input.soil)*1000~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(total.C.input.substrate~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate((total.N.input.substrate)*1000~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(total.C.output~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
aggregate((g.N.resin)*1000~Substrate*Management*Inoculation,data=dataEva,function(x) c(mean = mean(x), sd = sd(x)))
#make aggregated data frame
dataEvafunctions<-subset(dataEva,select=c("Substrate","Management","Inoculation","Field","Mass.loss","NO3H4uptake.resin","POXC.end","TOTALCO2produced.total","SIR.end"))
#dataEvafunctionssum<-aggregate.data.frame(dataEvafunctions[4:8],by=list(dataEvafunctions[1:3]),FUN=mean)
aggregate(Mass.loss~Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(NO3H4uptake.resin~Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(POXC.end~Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(TOTALCO2produced.total~Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(SIR.end~Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))

#for C and N also per field
Nresin.all.long<-aggregate(NO3H4uptake.resin~Field*Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x)))
Nresin.all.short<-dcast(data=Nresin.all.long,Field*Management*Inoculation~Substrate)
write.csv(Nresin.all.short,"Nresin.datatable.csv")
CumResp.all.long<-aggregate(TOTALCO2produced.total~Field*Substrate*Management*Inoculation,data=dataEvafunctions,function(x) c(mean = mean(x)))
CumResp.all.short<-dcast(data=CumResp.all.long,Field*Management*Inoculation~Substrate)
write.csv(CumResp.all.short,"CumResp.datatable.csv")
LM.CumResp<-lm(TOTALCO2produced.total~Substrate*Field,data=dataEvafunctions)
LM.CumResp.resid<-resid(LM.CumResp)
LM.CumResp.fitted<-fitted(LM.CumResp)
plot(LM.CumResp.resid~LM.CumResp.fitted)


##### ANOVAS #####
#real data
summary(aov(Mass.loss~Substrate*Management*Inoculation,data=dataEvasubs))
summary(aov(NO3H4uptake.resin~Substrate*Management*Inoculation,data=dataEva))
summary(aov(TOTALCO2produced.total~Substrate*Management*Inoculation,data=dataEva))
summary(aov(SIR.end~Substrate*Management*Inoculation,data=dataEva))
summary(aov(POXC.end~Substrate*Management*Inoculation,data=dataEva))

#post hoc plus cld
Mass.loss.aov <- aov(Mass.loss~Substrate*Management*Inoculation,data=dataEvasubs)
Mass.loss.posthoc<-lsmeans (Mass.loss.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(Mass.loss.posthoc)

CO2.aov <- aov(TOTALCO2produced.total~Substrate*Management*Inoculation,data=dataEva)
CO2.posthoc<-lsmeans (CO2.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(CO2.posthoc)

N.resin.aov <- aov(NO3H4uptake.resin~Substrate*Management*Inoculation,data=dataEva)
N.resin.posthoc<-lsmeans (N.resin.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(N.resin.posthoc) #only these are relevant

POXC.aov <- aov(POXC.end~Substrate*Management*Inoculation,data=dataEva)
POXC.posthoc<-lsmeans (POXC.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(POXC.posthoc)

SIR.aov <- aov(SIR.end~Substrate*Management*Inoculation,data=dataEva)
SIR.posthoc<-lsmeans (SIR.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(SIR.posthoc)


#field effects
summary(aov(Mass.loss~Field,data=dataEvasubs))
summary(aov(NO3H4uptake.resin~Field,data=dataEva))
summary(aov(TOTALCO2produced.total~Field,data=dataEva))
summary(aov(SIR.end~Field,data=dataEva))
summary(aov(POXC.end~Field,data=dataEva))

aggregate(SIR.end~Substrate,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(SIR.end~Management*Substrate,data=dataEvafunctions,function(x) c(mean = mean(x), sd = sd(x)))
TukeyHSD(aov(SIR.end~Substrate,data=dataEva))

#####  substrate effect : logresponse to blanc  #####
dElrB<-subset(dataEva,select=c("Code","Pair","Inoculation","Management","replicate","Organic.Soil.percentage","Substrate","C.N.substrate","Mass.loss","NO3H4uptake.resin","TOTALCO2produced.total","SIR.end","POXC.end"))
dESML<-dcast(dElrB,replicate*Management*Inoculation~Substrate,value.var="Mass.loss")
dESN<-dcast(dElrB,replicate*Management*Inoculation~Substrate,value.var="NO3H4uptake.resin")
dESCO2<-dcast(dElrB,replicate*Management*Inoculation~Substrate,value.var="TOTALCO2produced.total")
dESSIR<-dcast(dElrB,replicate*Management*Inoculation~Substrate,value.var="SIR.end")
dESPOXC<-dcast(dElrB,replicate*Management*Inoculation~Substrate,value.var="POXC.end")

dElrB2ML<-subset(dElrB,Substrate=="Straw",select=c("Code","Pair","Inoculation","Management","replicate"))
dElrB2ML$SE.CoverCrop<-dESML$`Cover Crop`
dElrB2ML$SE.CarrotLeaves<-dESML$'Carrot leaves'
dElrB2ML$SE.Alfalfa<-dESML$'Alfalfa'
dElrB2ML$SE.Hay<-dESML$'Hay'
dElrB2ML$SE.Straw<-dESML$'Straw'
dElrB2ML$Measurement<-"Mass loss"

dElrB2N<-subset(dElrB,Substrate=="Straw",select=c("Code","Pair","Inoculation","Management","replicate"))
dElrB2N$SE.CoverCrop<-log(dESN$'Cover Crop'/dESN$'Control')
dElrB2N$SE.CarrotLeaves<-log(dESN$'Carrot leaves'/dESN$'Control')
dElrB2N$SE.Alfalfa<-log(dESN$'Alfalfa'/dESN$'Control')
dElrB2N$SE.Hay<-log(dESN$'Hay'/dESN$'Control')
dElrB2N$SE.Straw<-log(dESN$'Straw'/dESN$'Control')
dElrB2N$Measurement<-"Cumulative N"

dElrB2CO2<-subset(dElrB,Substrate=="Straw",select=c("Code","Pair","Inoculation","Management","replicate"))
dElrB2CO2$SE.CoverCrop<-log(dESCO2$'Cover Crop'/dESCO2$'Control')
dElrB2CO2$SE.CarrotLeaves<-log(dESCO2$'Carrot leaves'/dESCO2$'Control')
dElrB2CO2$SE.Alfalfa<-log(dESCO2$'Alfalfa'/dESCO2$'Control')
dElrB2CO2$SE.Hay<-log(dESCO2$'Hay'/dESCO2$'Control')
dElrB2CO2$SE.Straw<-log(dESCO2$'Straw'/dESCO2$'Control')
dElrB2CO2$Measurement<-"Cumulative C"

dElrB2SIR<-subset(dElrB,Substrate=="Straw",select=c("Code","Pair","Inoculation","Management","replicate"))
dElrB2SIR$SE.CoverCrop<-log(dESSIR$'Cover Crop'/dESSIR$'Control')
dElrB2SIR$SE.CarrotLeaves<-log(dESSIR$'Carrot leaves'/dESSIR$'Control')
dElrB2SIR$SE.Alfalfa<-log(dESSIR$'Alfalfa'/dESSIR$'Control')
dElrB2SIR$SE.Hay<-log(dESSIR$'Hay'/dESSIR$'Control')
dElrB2SIR$SE.Straw<-log(dESSIR$'Straw'/dESSIR$'Control')
dElrB2SIR$Measurement<-"SIR.end"

dElrB2POXC<-subset(dElrB,Substrate=="Straw",select=c("Code","Pair","Inoculation","Management","replicate"))
dElrB2POXC$SE.CoverCrop<-log(dESPOXC$'Cover Crop'/dESPOXC$'Control')
dElrB2POXC$SE.CarrotLeaves<-log(dESPOXC$'Carrot leaves'/dESPOXC$'Control')
dElrB2POXC$SE.Alfalfa<-log(dESPOXC$'Alfalfa'/dESPOXC$'Control')
dElrB2POXC$SE.Hay<-log(dESPOXC$'Hay'/dESPOXC$'Control')
dElrB2POXC$SE.Straw<-log(dESPOXC$'Straw'/dESPOXC$'Control')
dElrB2POXC$Measurement<-"POXC.end"

dElrB3<-rbind(dElrB2N,dElrB2CO2,dElrB2SIR,dElrB2POXC)
colnames(dElrB3)<-c("Code","Pair","Inoculation","Management","replicate","Cover crop","Carrot leaves","Alfalfa","Hay","Straw","Measurement")
dElrB4<-melt(dElrB3)
colnames(dElrB4) <- c("Code","Pair","Inoculation","Management","replicate","Measurement","Substrate","Outcome")
dEOSP<-subset(dataEva,select=c("Code","Organic.Soil.percentage"))
dElrB4<-left_join(dElrB4,dEOSP)

ggplot(dElrB4,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate,scales="free")+
  theme_bw(base_size=15)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")

ggplot(dElrB4,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate,scales="free")+
  theme_bw(base_size=15)+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")

dElrB4$Substrate<-factor(dElrB4$Substrate)

dElrB4N.resin<-subset(dElrB4,Measurement=="Cumulative N")
dElrB4CO2<-subset(dElrB4,Measurement=="Cumulative C")
dElrB4SIR<-subset(dElrB4,Measurement=="SIR.end")
dElrB4POXC<-subset(dElrB4,Measurement=="POXC.end")

#general AOVs
summary(aov(data=dElrB4CO2,Outcome~Substrate*Management*Inoculation))
summary(aov(data=dElrB4N.resin,Outcome~Substrate*Management*Inoculation))
summary(aov(data=dElrB4POXC,Outcome~Substrate*Management*Inoculation))
summary(aov(data=dElrB4SIR,Outcome~Substrate*Management*Inoculation))

#tukey posthoc + cld


dElrB4CO2.aov <- aov(data=dElrB4CO2,Outcome~Substrate*Management*Inoculation)
dElrB4CO2.posthoc<-lsmeans (dElrB4CO2.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(dElrB4CO2.posthoc)

dElrB4N.resin.aov <- aov(data=dElrB4N.resin,Outcome~Substrate*Management*Inoculation)
dElrB4N.resin.posthoc<-lsmeans (dElrB4N.resin.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(dElrB4N.resin.posthoc) #only these are relevant

dElrB4POXC.aov <- aov(data=dElrB4POXC,Outcome~Substrate*Management*Inoculation)
dElrB4POXC.posthoc<-lsmeans (dElrB4POXC.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(dElrB4POXC.posthoc)

dElrB4SIR.aov <- aov(data=dElrB4SIR,Outcome~Substrate*Management*Inoculation)
dElrB4SIR.posthoc<-lsmeans (dElrB4SIR.aov, pairwise ~ Substrate*Management*Inoculation,adjust="tukey")
cld(dElrB4SIR.posthoc)

#create overview for general substrate effects 
TukeyHSD(aov(data=dElrB4CO2,Outcome~Substrate))
dElrB4CO2.aov <- aov(data=dElrB4CO2,Outcome~Substrate)
dElrB4CO2.posthoc<-lsmeans (dElrB4CO2.aov, pairwise ~ Substrate,adjust="tukey")
cld(dElrB4CO2.posthoc)

dElrB4POXC.aov <- aov(data=dElrB4POXC,Outcome~Substrate)
dElrB4POXC.posthoc<-lsmeans (dElrB4POXC.aov, pairwise ~ Substrate,adjust="tukey")
cld(dElrB4POXC.posthoc)

dElrB4SIR.aov <- aov(data=dElrB4SIR,Outcome~Substrate)
dElrB4SIR.posthoc<-lsmeans (dElrB4SIR.aov, pairwise ~ Substrate,adjust="tukey")
cld(dElrB4SIR.posthoc)

#different from 0
dElrB4N.resinCoverCrop<-subset(dElrB4N.resin,Substrate=="Cover crop")
dElrB4N.resinCarrotLeaves<-subset(dElrB4N.resin,Substrate=="Carrot leaves")
dElrB4N.resinAlfalfa<-subset(dElrB4N.resin,Substrate=="Alfalfa")
dElrB4N.resinHay<-subset(dElrB4N.resin,Substrate=="Hay")
dElrB4N.resinStraw<-subset(dElrB4N.resin,Substrate=="Straw")
t.test(dElrB4N.resinCoverCrop$Outcome)
t.test(dElrB4N.resinCarrotLeaves$Outcome)
t.test(dElrB4N.resinAlfalfa$Outcome)
t.test(dElrB4N.resinHay$Outcome)
t.test(dElrB4N.resinStraw$Outcome)

# 
# summary(aov(data=dElrB4N.resinCoverCrop,Outcome~Management*Inoculation))
# summary(aov(data=dElrB4N.resinCarrotLeaves,Outcome~Management*Inoculation))
# summary(aov(data=dElrB4N.resinAlfalfa,Outcome~Management*Inoculation))
# summary(aov(data=dElrB4N.resinHay,Outcome~Management*Inoculation))
# summary(aov(data=dElrB4N.resinStraw,Outcome~Management*Inoculation))
# 
# TukeyHSD(aov(data=dElrB4N.resinCoverCrop,Outcome~Organic.Soil.percentage))
# TukeyHSD(aov(data=dElrB4N.resinCarrotLeaves,Outcome~Organic.Soil.percentage))
# TukeyHSD(aov(data=dElrB4N.resinAlfalfa,Outcome~Organic.Soil.percentage))
# TukeyHSD(aov(data=dElrB4N.resinHay,Outcome~Organic.Soil.percentage))
# TukeyHSD(aov(data=dElrB4N.resinStraw,Outcome~Organic.Soil.percentage))
# summary(aov(data=dElrB4N.resinStraw,Outcome~Management*Inoculation))



dElrB4CO2CoverCrop<-subset(dElrB4CO2,Substrate=="Cover crop")
dElrB4CO2CarrotLeaves<-subset(dElrB4CO2,Substrate=="Carrot leaves")
dElrB4CO2Alfalfa<-subset(dElrB4CO2,Substrate=="Alfalfa")
dElrB4CO2Hay<-subset(dElrB4CO2,Substrate=="Hay")
dElrB4CO2Straw<-subset(dElrB4CO2,Substrate=="Straw")

t.test(dElrB4CO2CoverCrop$Outcome)
t.test(dElrB4CO2CarrotLeaves$Outcome)
t.test(dElrB4CO2Alfalfa$Outcome)
t.test(dElrB4CO2Hay$Outcome)
t.test(dElrB4CO2Straw$Outcome)

summary(aov(data=dElrB4SIR,Outcome~Substrate*Management*Inoculation))
TukeyHSD(aov(data=dElrB4SIR,Outcome~Substrate))
TukeyHSD(aov(data=dElrB4SIR,Outcome~Organic.Soil.percentage))
dElrB4SIRCoverCrop<-subset(dElrB4SIR,Substrate=="Cover crop")
dElrB4SIRCarrotLeaves<-subset(dElrB4SIR,Substrate=="Carrot leaves")
dElrB4SIRAlfalfa<-subset(dElrB4SIR,Substrate=="Alfalfa")
dElrB4SIRHay<-subset(dElrB4SIR,Substrate=="Hay")
dElrB4SIRStraw<-subset(dElrB4SIR,Substrate=="Straw")
TukeyHSD(aov(data=dElrB4SIRCoverCrop,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4SIRCarrotLeaves,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4SIRAlfalfa,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4SIRHay,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4SIRStraw,Outcome~Organic.Soil.percentage))
t.test(dElrB4SIRCoverCrop$Outcome)
t.test(dElrB4SIRCarrotLeaves$Outcome)
t.test(dElrB4SIRAlfalfa$Outcome)
t.test(dElrB4SIRHay$Outcome)
t.test(dElrB4SIRStraw$Outcome)



TukeyHSD(aov(data=dElrB4POXC,Outcome~Substrate))
TukeyHSD(aov(data=dElrB4POXC,Outcome~Organic.Soil.percentage))
dElrB4POXCCoverCrop<-subset(dElrB4POXC,Substrate=="Cover crop")
dElrB4POXCCarrotLeaves<-subset(dElrB4POXC,Substrate=="Carrot leaves")
dElrB4POXCAlfalfa<-subset(dElrB4POXC,Substrate=="Alfalfa")
dElrB4POXCHay<-subset(dElrB4POXC,Substrate=="Hay")
dElrB4POXCStraw<-subset(dElrB4POXC,Substrate=="Straw")
TukeyHSD(aov(data=dElrB4POXCCoverCrop,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4POXCCarrotLeaves,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4POXCAlfalfa,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4POXCHay,Outcome~Organic.Soil.percentage))
TukeyHSD(aov(data=dElrB4POXCStraw,Outcome~Organic.Soil.percentage))
t.test(dElrB4POXCCoverCrop$Outcome)
t.test(dElrB4POXCCarrotLeaves$Outcome)
t.test(dElrB4POXCAlfalfa$Outcome)
t.test(dElrB4POXCHay$Outcome)
t.test(dElrB4POXCStraw$Outcome)

ggplot(dElrB4,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate,scales="free")+
  theme_bw(base_size=20)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")


dElrB4CO2<-subset(dElrB4,Measurement=="TOTALCO2produced.total")
dElrB4SIR<-subset(dElrB4,Measurement=="SIR.end")
dElrB4POXC<-subset(dElrB4,Measurement=="POXC.end")

ggplot(dElrB4N.resin,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate)+
  theme_bw(base_size=20)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")

ggplot(dElrB4CO2,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate)+
  theme_bw(base_size=20)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")

ggplot(dElrB4SIR,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate)+
  theme_bw(base_size=20)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")


ggplot(dElrB4POXC,aes(y=Outcome,x=Management,fill=Management,pattern=Inoculation))+
  geom_boxplot_pattern(color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(Measurement~Substrate)+
  theme_bw(base_size=20)+
  labs(y="Effect size of substrate \n(log (Substrate/Blanc))")






data1 <- summarySE(dataEvasubs, measurevar="Mass.loss", groupvars=c("Management","Inoculation","Substrate"), na.rm=TRUE)
ggplot(data=data1, aes(x=Management, y=Mass.loss,fill=Management,pattern=Inoculation)) +
  geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  facet_grid(.~Substrate)+
  geom_errorbar(aes(ymin=Mass.loss, ymax=Mass.loss+se),
                width=.2, size =.3,                   # Width of the error bars
                position=position_dodge(.9)) 



ggplot(data=data1, aes(x=Management, y=Mass.loss,fill=Management,pattern=Inoculation)) +
  geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),
                   color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_fill_manual(values=c("white","grey")) +
  scale_pattern_manual(values = c(no = "none", yes = "stripe")) +
  facet_grid(.~Substrate)+
  geom_errorbar(aes(ymin=Mass.loss, ymax=Mass.loss+se),
                width=.2, size =.3,                   # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

