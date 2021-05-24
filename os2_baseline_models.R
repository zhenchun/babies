#Packages
library(dlnm)
library(lme4)
library(lmerTest)
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(stringr)
library(dplyr)
library(readr)
library(robustHD)

path<-"C:/Users/zy125/Box Sync/Postdoc/os2/data"

setwd(path)

getwd()

Fulldata=read.csv("OS2.csv")
View(Fulldata)
Baseline=subset(Fulldata,Post==0)

###section A####
##########################################################################################
####divide distance to road length into 4 groups########################################## 
####group four as the reference and then use other three groups to compare with group 4###
##########################################################################################


##AA is a baseline database include all variables for baseline level analysis
AA=data.frame()
AA[1:238,1]=Baseline$ANAP1_Cr
AA[1:238,2]=Baseline$ANAP2_Cr
AA[1:238,3]=Baseline$AFLU2_Cr
AA[1:238,4]=Baseline$APHE9_Cr
AA[1:238,5]=Baseline$APYR1_Cr
AA[1:238,6]=Baseline$TAPAHs_Cr
AA[1:238,7]=Baseline$COPD
AA[1:238,8]=Baseline$IHD
AA[1:238,9]=Baseline$Age
AA[1:238,10]=Baseline$Male
AA[1:238,11]=Baseline$BMI
AA[1:238,12]=Baseline[,1]
AA[1:238,13]=Baseline$all_roads_length_100m
AA[1:238,14]=Baseline$major_roads_length_100m
AA[1:238,15]=Baseline$minor_roads_length_100m
AA[1:238,16]=Baseline$Q1
AA[1:238,17]=Baseline$Q2
AA[1:238,18]=Baseline$Q3
AA[1:238,19]=Baseline$Q4
AA[1:238,20]=Baseline$Near100
AA[1:238,21]=Baseline$Far100
AA[1:238,22]=Baseline$Group
AA[1:238,23]=Baseline$Gender
AA[1:238,24]=Baseline$NO2_2013

AA=na.omit(AA)
colnames(AA)<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id","all_roads_length_100m",
                "major_roads_length_100m","minor_roads_length_100m","Q1", "Q2","Q3","Q4", 
                "Near100","Far100", "Group", "Gender", "NO2_2013")

bio_names<- c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

###section A#########################################################
####test the distance to roads effect on concentrations##

dist.coefs<-list()


for (i in 1:6) {
  
  lmerA=lmer(log10(AA[,i])~AA[,16]+AA[,17]+AA[,18]+AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  dist.coefs[[i]]<-data.frame(coef(summary(lmerA)))# extract coefficients
  dist.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(dist.coefs[[i]]$t.value)))
  dist.coefs[[i]]$LCI=dist.coefs[[i]]$Estimate-1.96*dist.coefs[[i]]$Std..Error
  dist.coefs[[i]]$UCI=dist.coefs[[i]]$Estimate+1.96*dist.coefs[[i]]$Std..Error
  
  dist.coefs[[i]]$FC=10^(dist.coefs[[i]]$Estimate)
  dist.coefs[[i]]$LFC=10^(dist.coefs[[i]]$LCI)
  dist.coefs[[i]]$UFC=10^(dist.coefs[[i]]$UCI)
  dist.coefs[[i]]$PC=((10^(dist.coefs[[i]]$Estimate))-1)*100
  dist.coefs[[i]]$LPC=((10^(dist.coefs[[i]]$LCI))-1)*100
  dist.coefs[[i]]$UPC=((10^(dist.coefs[[i]]$UCI))-1)*100
}

names(dist.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(dist.coefs[[i]], paste0("ModelA_", bio_names[i],".csv"))
}

#############################################################################
###section B#################################################################
####test the in or outside 100m of a major road effect on concentrations#####
#############################################################################
M100.coefs<-list()


for (i in 1:6) {
  
  lmerB=lmer(log10(AA[,i])~AA[,20]+AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  M100.coefs[[i]]<-data.frame(coef(summary(lmerB)))# extract coefficients
  M100.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(M100.coefs[[i]]$t.value)))#
  M100.coefs[[i]]$LCI=M100.coefs[[i]]$Estimate-1.96*M100.coefs[[i]]$Std..Error
  M100.coefs[[i]]$UCI=M100.coefs[[i]]$Estimate+1.96*M100.coefs[[i]]$Std..Error
  
  M100.coefs[[i]]$FC=10^(M100.coefs[[i]]$Estimate)
  M100.coefs[[i]]$LFC=10^(M100.coefs[[i]]$LCI)
  M100.coefs[[i]]$UFC=10^(M100.coefs[[i]]$UCI)
  M100.coefs[[i]]$PC=((10^(M100.coefs[[i]]$Estimate))-1)*100
  M100.coefs[[i]]$LPC=((10^(M100.coefs[[i]]$LCI))-1)*100
  M100.coefs[[i]]$UPC=((10^(M100.coefs[[i]]$UCI))-1)*100
}

names(M100.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(M100.coefs[[i]], paste0("ModelB_", bio_names[i],".csv"))
}

