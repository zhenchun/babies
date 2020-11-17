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
  dist.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(dist.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
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
  M100.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(M100.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
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

#######################################################################
###section C############################################################
####test the roadlength within 100m buffers effect on concentrations####
#######################################################################


####major roads
major.length.coefs<-list()


for (i in 1:6) {
  
  lmerCMA=lmer(log10(AA[,i])~AA[,14]+AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  major.length.coefs[[i]]<-data.frame(coef(summary(lmerCMA)))# extract coefficients
  major.length.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(major.length.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  major.length.coefs[[i]]$LCI=major.length.coefs[[i]]$Estimate-1.96*major.length.coefs[[i]]$Std..Error
  major.length.coefs[[i]]$UCI=major.length.coefs[[i]]$Estimate+1.96*major.length.coefs[[i]]$Std..Error
  
  major.length.coefs[[i]]$FC=10^(major.length.coefs[[i]]$Estimate)
  major.length.coefs[[i]]$LFC=10^(major.length.coefs[[i]]$LCI)
  major.length.coefs[[i]]$UFC=10^(major.length.coefs[[i]]$UCI)
  major.length.coefs[[i]]$PC=((10^(major.length.coefs[[i]]$Estimate))-1)*100
  major.length.coefs[[i]]$LPC=((10^(major.length.coefs[[i]]$LCI))-1)*100
  major.length.coefs[[i]]$UPC=((10^(major.length.coefs[[i]]$UCI))-1)*100
  major.length.coefs[[i]]$names=bio_names[i]
  major.length.coefs[[i]]$type="major.road"
}

major.length<-data.frame()

for (i in 1:6){
  major.length<-rbind(major.length, major.length.coefs[[i]][2,])
}


####minor roads

minor.length.coefs<-list()


for (i in 1:6) {
  
  lmerCMI=lmer(log10(AA[,i])~AA[,15]+AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  minor.length.coefs[[i]]<-data.frame(coef(summary(lmerCMI)))# extract coefficients
  minor.length.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(minor.length.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  minor.length.coefs[[i]]$LCI=minor.length.coefs[[i]]$Estimate-1.96*minor.length.coefs[[i]]$Std..Error
  minor.length.coefs[[i]]$UCI=minor.length.coefs[[i]]$Estimate+1.96*minor.length.coefs[[i]]$Std..Error
  
  minor.length.coefs[[i]]$FC=10^(minor.length.coefs[[i]]$Estimate)
  minor.length.coefs[[i]]$LFC=10^(minor.length.coefs[[i]]$LCI)
  minor.length.coefs[[i]]$UFC=10^(minor.length.coefs[[i]]$UCI)
  minor.length.coefs[[i]]$PC=((10^(minor.length.coefs[[i]]$Estimate))-1)*100
  minor.length.coefs[[i]]$LPC=((10^(minor.length.coefs[[i]]$LCI))-1)*100
  minor.length.coefs[[i]]$UPC=((10^(minor.length.coefs[[i]]$UCI))-1)*100
  minor.length.coefs[[i]]$names=bio_names[i]
  minor.length.coefs[[i]]$type="minor.road"
}

names(minor.length.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")


minor.length<-data.frame()

for (i in 1:6){
  minor.length<-rbind(minor.length, minor.length.coefs[[i]][2,])
}


####all roads


all.length.coefs<-list()


for (i in 1:6) {
  
  lmerCAL=lmer(log10(AA[,i])~AA[,13]+AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  all.length.coefs[[i]]<-data.frame(coef(summary(lmerCAL)))# extract coefficients
  all.length.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(all.length.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  all.length.coefs[[i]]$LCI=all.length.coefs[[i]]$Estimate-1.96*all.length.coefs[[i]]$Std..Error
  all.length.coefs[[i]]$UCI=all.length.coefs[[i]]$Estimate+1.96*all.length.coefs[[i]]$Std..Error
  
  all.length.coefs[[i]]$FC=10^(all.length.coefs[[i]]$Estimate)
  all.length.coefs[[i]]$LFC=10^(all.length.coefs[[i]]$LCI)
  all.length.coefs[[i]]$UFC=10^(all.length.coefs[[i]]$UCI)
  all.length.coefs[[i]]$PC=((10^(all.length.coefs[[i]]$Estimate))-1)*100
  all.length.coefs[[i]]$LPC=((10^(all.length.coefs[[i]]$LCI))-1)*100
  all.length.coefs[[i]]$UPC=((10^(all.length.coefs[[i]]$UCI))-1)*100
  all.length.coefs[[i]]$names=bio_names[i]
  all.length.coefs[[i]]$type="all.road"
}

names(all.length.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")



all.length<-data.frame()

for (i in 1:6){
  all.length<-rbind(all.length, all.length.coefs[[i]][2,])
}


road_length<-rbind(major.length,minor.length,all.length)####combine different roads types
road_length$type<-factor(road_length$type, levels=c("major.road","minor.road","all.road"))
road_length$names<-factor(road_length$names, c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs"))

write.csv(road_length, "road_length.csv")

#######code to plot length######################################### 

pd <- position_dodge(width = 1)

p1<-ggplot(road_length, aes(y=PC*10, x=names))+
  geom_pointrange(aes(ymax=UPC*10, ymin=LPC*10, color=type), position = pd, size=0.8)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5),linetype="dotted", size=1)+
  geom_hline(yintercept = 0)+theme_classic()+
  scale_color_manual(values=c("black","deepskyblue","coral"),labels = c("Major Roads","Minor Roads","All Roads"))+
  labs(x="biomarker names", y="Percentage Change (%, 95%CI)")+
  scale_x_discrete(labels=c("1-ANAP", "2-ANAP", "2-AFLU", "9-APHE", "1-APYR", "TAPAHs"))+
  geom_text(x=2.67, y=7.75, label="*", size=5)+
  theme(legend.title=element_blank(),legend.position="bottom",axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.text=element_text(size=14), 
        axis.title.x=element_blank(), axis.title.y = element_text(size=16),legend.key.size=unit(4, "points"))

##############################################################################################
####test the roadlength within 100m buffers effect on concentrations modified by sex##########
##############################################################################################
##############################################################################################

major.length.sex.coefs<-list() 

for (i in 1:6) {
     
  lmerCMAS=lmer(log10(AA[,i])~AA[,14]+AA[,7]+AA[,8]+AA[,9]+AA[,11]+(AA[,14]*AA[,10])+(1|AA[,12]))
  major.length.sex.coefs[[i]]<-data.frame(coef(summary(lmerCMAS)))# extract coefficients
  major.length.sex.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(major.length.sex.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  major.length.sex.coefs[[i]]$LCI=major.length.sex.coefs[[i]]$Estimate-1.96*major.length.sex.coefs[[i]]$Std..Error
  major.length.sex.coefs[[i]]$UCI=major.length.sex.coefs[[i]]$Estimate+1.96*major.length.sex.coefs[[i]]$Std..Error
       
  major.length.sex.coefs[[i]]$FC=10^(major.length.sex.coefs[[i]]$Estimate)
  major.length.sex.coefs[[i]]$LFC=10^(major.length.sex.coefs[[i]]$LCI)
  major.length.sex.coefs[[i]]$UFC=10^(major.length.sex.coefs[[i]]$UCI)
  major.length.sex.coefs[[i]]$PC=((10^(major.length.sex.coefs[[i]]$Estimate))-1)*100
  major.length.sex.coefs[[i]]$LPC=((10^(major.length.sex.coefs[[i]]$LCI))-1)*100
  major.length.sex.coefs[[i]]$UPC=((10^(major.length.sex.coefs[[i]]$UCI))-1)*100
  major.length.sex.coefs[[i]]$names=bio_names[i]
  major.length.sex.coefs[[i]]$type="major.road"
              
  }

major.length.sex<-data.frame()

for (i in 1:6){
  major.length.sex<-rbind(major.length.sex, major.length.sex.coefs[[i]][8,])
}


minor.length.sex.coefs<-list() 

for (i in 1:6) {
  
  lmerCMIS=lmer(log10(AA[,i])~AA[,15]+AA[,7]+AA[,8]+AA[,9]+AA[,11]+AA[,15]*AA[,10]+(1|AA[,12]))
  minor.length.sex.coefs[[i]]<-data.frame(coef(summary(lmerCMIS)))# extract coefficients
  minor.length.sex.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(minor.length.sex.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  minor.length.sex.coefs[[i]]$LCI=minor.length.sex.coefs[[i]]$Estimate-1.96*minor.length.sex.coefs[[i]]$Std..Error
  minor.length.sex.coefs[[i]]$UCI=minor.length.sex.coefs[[i]]$Estimate+1.96*minor.length.sex.coefs[[i]]$Std..Error
  
  minor.length.sex.coefs[[i]]$FC=10^(minor.length.sex.coefs[[i]]$Estimate)
  minor.length.sex.coefs[[i]]$LFC=10^(minor.length.sex.coefs[[i]]$LCI)
  minor.length.sex.coefs[[i]]$UFC=10^(minor.length.sex.coefs[[i]]$UCI)
  minor.length.sex.coefs[[i]]$PC=((10^(minor.length.sex.coefs[[i]]$Estimate))-1)*100
  minor.length.sex.coefs[[i]]$LPC=((10^(minor.length.sex.coefs[[i]]$LCI))-1)*100
  minor.length.sex.coefs[[i]]$UPC=((10^(minor.length.sex.coefs[[i]]$UCI))-1)*100
  minor.length.sex.coefs[[i]]$names=bio_names[i]
  minor.length.sex.coefs[[i]]$type="minor.road"
  
  }


minor.length.sex<-data.frame()

for (i in 1:6){
  minor.length.sex<-rbind(minor.length.sex, minor.length.sex.coefs[[i]][8,])
}


all.length.sex.coefs<-list() 

for (i in 1:6) {
  
  lmerCAS=lmer(log10(AA[,i])~AA[,13]+AA[,7]+AA[,8]+AA[,9]+AA[,11]+(AA[,13]*AA[,10])+(1|AA[,12]))
  all.length.sex.coefs[[i]]<-data.frame(coef(summary(lmerCAS)))# extract coefficients
  all.length.sex.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(all.length.sex.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  all.length.sex.coefs[[i]]$LCI=all.length.sex.coefs[[i]]$Estimate-1.96*all.length.sex.coefs[[i]]$Std..Error
  all.length.sex.coefs[[i]]$UCI=all.length.sex.coefs[[i]]$Estimate+1.96*all.length.sex.coefs[[i]]$Std..Error
  
  all.length.sex.coefs[[i]]$FC=10^(all.length.sex.coefs[[i]]$Estimate)
  all.length.sex.coefs[[i]]$LFC=10^(all.length.sex.coefs[[i]]$LCI)
  all.length.sex.coefs[[i]]$UFC=10^(all.length.sex.coefs[[i]]$UCI)
  all.length.sex.coefs[[i]]$PC=((10^(all.length.sex.coefs[[i]]$Estimate))-1)*100
  all.length.sex.coefs[[i]]$LPC=((10^(all.length.sex.coefs[[i]]$LCI))-1)*100
  all.length.sex.coefs[[i]]$UPC=((10^(all.length.sex.coefs[[i]]$UCI))-1)*100
  all.length.sex.coefs[[i]]$names=bio_names[i]
  all.length.sex.coefs[[i]]$type="all.road"
}


all.length.sex<-data.frame()

for (i in 1:6){
  all.length.sex<-rbind(all.length.sex, all.length.sex.coefs[[i]][8,])
}


road_length_sex<-rbind(major.length.sex,minor.length.sex,all.length.sex)####combine different roads types
road_length_sex$type<-factor(road_length_sex$type, levels=c("major.road","minor.road","all.road"))
road_length_sex$names<-factor(road_length_sex$names, c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs"))

write.csv(road_length_sex, "road_length_sex.csv")


###############################################################################
#############base mixed effect models with random intercepts of subjects.###### 
#####Fixed effects includes COPD, IHD, age, sex, and BMI. #####################
###############################################################################
base.coefs<-list()
 
   
for (i in 1:6) {
      
  lmerBA=lmer(log10(AA[,i])~AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+(1|AA[,12]))
  base.coefs[[i]]<-data.frame(coef(summary(lmerBA)))# extract coefficients
  base.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(base.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  base.coefs[[i]]$LCI=base.coefs[[i]]$Estimate-1.96*base.coefs[[i]]$Std..Error
  base.coefs[[i]]$UCI=base.coefs[[i]]$Estimate+1.96*base.coefs[[i]]$Std..Error
        
  base.coefs[[i]]$FC=10^(base.coefs[[i]]$Estimate)
  base.coefs[[i]]$LFC=10^(base.coefs[[i]]$LCI)
  base.coefs[[i]]$UFC=10^(base.coefs[[i]]$UCI)
  base.coefs[[i]]$PC=((10^(base.coefs[[i]]$Estimate))-1)*100
  base.coefs[[i]]$LPC=((10^(base.coefs[[i]]$LCI))-1)*100
  base.coefs[[i]]$UPC=((10^(base.coefs[[i]]$UCI))-1)*100
  base.coefs[[i]]$names=bio_names[i]
         
}


base<-data.frame()

for (i in 1:6){
  base<-rbind(base, base.coefs[[i]][8,])
}
 
base_COPD<-data.frame()
 
for (i in 1:6){
       base_COPD<-rbind(base_COPD, base.coefs[[i]][2,])
     }
 
base_IHD<-data.frame()
 
for (i in 1:6){
       base_IHD<-rbind(base_IHD, base.coefs[[i]][3,])
   }

   
##################################section C#######################################################################
###################################test the long-term no2 concentrations effects on biomarkers' concentrations####
##################################################################################################################
base.no2.coefs<-list()
    
      
for (i in 1:6) {
         
  lmerBAN=lmer(log10(AA[,i])~AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+AA[,24]+(1|AA[,12]))
  base.no2.coefs[[i]]<-data.frame(coef(summary(lmerBAN)))# extract coefficients
  base.no2.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(base.no2.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  base.no2.coefs[[i]]$LCI=base.no2.coefs[[i]]$Estimate-1.96*base.no2.coefs[[i]]$Std..Error
  base.no2.coefs[[i]]$UCI=base.no2.coefs[[i]]$Estimate+1.96*base.no2.coefs[[i]]$Std..Error
      
  base.no2.coefs[[i]]$FC=10^(base.no2.coefs[[i]]$Estimate)
  base.no2.coefs[[i]]$LFC=10^(base.no2.coefs[[i]]$LCI)
  base.no2.coefs[[i]]$UFC=10^(base.no2.coefs[[i]]$UCI)
  base.no2.coefs[[i]]$PC=((10^(base.no2.coefs[[i]]$Estimate))-1)*100
  base.no2.coefs[[i]]$LPC=((10^(base.no2.coefs[[i]]$LCI))-1)*100
  base.no2.coefs[[i]]$UPC=((10^(base.no2.coefs[[i]]$UCI))-1)*100
  base.no2.coefs[[i]]$names=bio_names[i]
          
}


base.no2<-data.frame()

for (i in 1:6){
  base.no2<-rbind(base.no2, base.no2.coefs[[i]][7,])
}

base.no2$PC.IQR<-IQR(AA$NO2)*base.no2$PC
base.no2$LPC.IQR<-IQR(AA$NO2)*base.no2$LPC
base.no2$UPC.IQR<-IQR(AA$NO2)*base.no2$UPC

base.no2$names<-factor(base.no2$names, c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs"))

write.csv(base.no2, "base.no2.csv")
 
base.no2_COPD<-data.frame()

for (i in 1:6){
base.no2_COPD<-rbind(base.no2_COPD, base.no2.coefs[[i]][2,])
       }

base.no2_IHD<-data.frame()

for (i in 1:6){
base.no2_IHD<-rbind(base.no2_IHD, base.no2.coefs[[i]][3,])
       }

write.csv(base.no2_COPD, "base.no2_COPD.csv") 

write.csv(base.no2_IHD, "base.no2_IHD.csv") 

p2<- ggplot(base.no2, aes(y=PC.IQR, x=names))+
     geom_pointrange(aes(ymax=UPC.IQR, ymin=LPC.IQR), position = pd, size=0.8)+
     geom_hline(yintercept = 0, linetype="dotted")+theme_classic()+
     labs(x="biomarker names", y="Percentage Change (%, 95%CI)")+
     scale_x_discrete(labels=c("1-ANAP", "2-ANAP", "2-AFLU", "9-APHE", "1-APYR", "TAPAHs"))+
     scale_y_continuous(limits = c(-60, 140),breaks=seq(-60,140, by=20))+
       theme(axis.text= element_text(size=14), 
                         axis.title.x=element_blank(), axis.title.y = element_text(size=16))



AA$Group<-factor(AA$Group, c("HEALTHY","COPD","IHD"))

 p3<-ggplot(AA, aes(x=Group, y=ANAP1_Cr, fill=Group))+geom_boxplot()+
   scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+ 
   theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.position="none")+ggtitle("A.1-ANAP")+
   stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
   annotate("rect", xmin = 1, xmax = 2, ymin = 110, ymax =110, alpha=1,colour = "black")+
   annotate("rect", xmin = 1, xmax = 1, ymin = 90, ymax =110, alpha=1,colour = "black")+
   annotate("rect", xmin = 2, xmax = 2, ymin = 90, ymax =110, alpha=1,colour = "black")+
   geom_text(aes(x=1.5, y=150, label="*", size=8))+
   annotate("rect", xmin = 1, xmax = 3, ymin = 350, ymax =350, alpha=1,colour = "black")+
   annotate("rect", xmin = 1, xmax = 1, ymin = 250, ymax =350, alpha=1,colour = "black")+
   annotate("rect", xmin = 3, xmax = 3, ymin = 250, ymax =350, alpha=1,colour = "black")+
   geom_text(aes(x=2, y=400, label="***", size=8))
 
 
  p4<-ggplot(AA, aes(x=Group, y=ANAP2_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none")+ggtitle("B.2-ANAP")+
    stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
    annotate("rect", xmin = 1, xmax = 2, ymin = 200, ymax =200, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin=150 , ymax =200, alpha=1,colour = "black")+
    annotate("rect", xmin = 2, xmax = 2, ymin = 150, ymax =200, alpha=1,colour = "black")+
    geom_text(aes(x=1.5, y=250, label="*", size=8))+
    annotate("rect", xmin = 1, xmax = 3, ymin = 550, ymax =550, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin = 400, ymax =550, alpha=1,colour = "black")+
    annotate("rect", xmin = 3, xmax = 3, ymin = 400, ymax =550, alpha=1,colour = "black")+
    geom_text(aes(x=2, y=600, label="***", size=8))
  
  
  p5<-ggplot(AA, aes(x=Group, y=AFLU2_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+ 
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none")+ggtitle("C.2-AFLU")+
    stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)
  
  
  p6<-ggplot(AA, aes(x=Group, y=APHE9_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000,10000),labels=c(0.001, 0.01,0.1,1,10,100,1000,10000),limits=c(0.001,10000))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none")+ggtitle("D.9-APHE")+
    stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
    annotate("rect", xmin = 1, xmax = 2, ymin = 1000, ymax =1000, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin = 700, ymax =1000, alpha=1,colour = "black")+
    annotate("rect", xmin = 2, xmax = 2, ymin = 700, ymax =1000, alpha=1,colour = "black")+
    geom_text(aes(x=1.5, y=1100, label="*", size=8))
  
  
  
  p7<-ggplot(AA, aes(x=Group, y=APYR1_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100),labels=c(0.001, 0.01,0.1,1,10,100),limits=c(0.001,100))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none")+ggtitle("E.1-APYR")+
    stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
    annotate("rect", xmin = 1, xmax = 3, ymin = 10, ymax =10, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin = 7, ymax =10, alpha=1,colour = "black")+
    annotate("rect", xmin = 3, xmax = 3, ymin = 7, ymax =10, alpha=1,colour = "black")+
    geom_text(aes(x=2, y=13, label="***", size=8))
  
  p8<-ggplot(AA, aes(x=Group, y=TAPAHs_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),labels=c(0.01,0.1,1,10,100,1000, 10000),limits=c(0.01,10000))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none")+ggtitle("F.\u03A3TAPAHs")+
    stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
    annotate("rect", xmin = 1, xmax = 2, ymin = 1000, ymax =1000, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin = 800, ymax =1000, alpha=1,colour = "black")+
    annotate("rect", xmin = 2, xmax = 2, ymin = 800, ymax =1000, alpha=1,colour = "black")+
    geom_text(aes(x=1.5, y=1100, label="**", size=8))+
    annotate("rect", xmin = 1, xmax = 3, ymin = 3500, ymax =3500, alpha=1,colour = "black")+
    annotate("rect", xmin = 1, xmax = 1, ymin = 2500, ymax =3500, alpha=1,colour = "black")+
    annotate("rect", xmin = 3, xmax = 3, ymin = 2500, ymax =3500, alpha=1,colour = "black")+
    geom_text(aes(x=2, y=4000, label="***", size=8))
  
  ####p9 is not in the final figure
  p9<-ggplot(AA, aes(x=Group, y=TAPAHs_Cr, fill=Group))+geom_boxplot()+
    scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),labels=c(0.01,0.1,1,10,100,1000, 10000),limits=c(0.01,10000))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                     axis.title.y=element_blank(), legend.position = "bottom")+ggtitle("F.\u03A3TAPAHs")
    
  
  #p9 is only used for extract legend
  
  
  ###function to extract the legend
  extract_legend <- function(my_ggp) {
        step1 <- ggplot_gtable(ggplot_build(my_ggp))
        step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
        step3 <- step1$grobs[[step2]]
        return(step3)
   }
  shared_legend <- extract_legend(p9)
  
  grid.arrange(arrangeGrob(p3, p4,p5, p6, p7, p8,nrow=2, 
      left=textGrob("Urinary concentrations (\u03BCg/g creatinine)", 
                    gp = gpar(fontsize = 12, fontface = 'bold'),rot = 90, vjust = 1)), shared_legend,heights = c(10, 1))
  

  
  
##################################################################################################################################
#################################models adjusted by external exposure of distance to road(divided by 100m) and NO2################ 
##################################################################################################################################  
    
base.no2.100.coefs<-list()
  
  
for (i in 1:6) {
    
    lmerBAN=lmer(log10(AA[,i])~AA[,7]+AA[,8]+AA[,9]+AA[,10]+AA[,11]+AA[,20]+AA[,24]+(1|AA[,12]))
    base.no2.100.coefs[[i]]<-data.frame(coef(summary(lmerBAN)))# extract coefficients
    base.no2.100.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(base.no2.100.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
    base.no2.100.coefs[[i]]$LCI=base.no2.100.coefs[[i]]$Estimate-1.96*base.no2.100.coefs[[i]]$Std..Error
    base.no2.100.coefs[[i]]$UCI=base.no2.100.coefs[[i]]$Estimate+1.96*base.no2.100.coefs[[i]]$Std..Error
    
    base.no2.100.coefs[[i]]$FC=10^(base.no2.100.coefs[[i]]$Estimate)
    base.no2.100.coefs[[i]]$LFC=10^(base.no2.100.coefs[[i]]$LCI)
    base.no2.100.coefs[[i]]$UFC=10^(base.no2.100.coefs[[i]]$UCI)
    base.no2.100.coefs[[i]]$PC=((10^(base.no2.100.coefs[[i]]$Estimate))-1)*100
    base.no2.100.coefs[[i]]$LPC=((10^(base.no2.100.coefs[[i]]$LCI))-1)*100
    base.no2.100.coefs[[i]]$UPC=((10^(base.no2.100.coefs[[i]]$UCI))-1)*100
    base.no2.100.coefs[[i]]$names=bio_names[i]
    
  }
  
  
base.no2.100_100<-data.frame()
  
for (i in 1:6){
    base.no2.100_100<-rbind(base.no2.100_100, base.no2.100.coefs[[i]][7,])
  }

write.csv(base.no2.100_100, "base.no2.100_100.csv")  

base.no2.100_no2<-data.frame()

for (i in 1:6){
    base.no2.100_no2<-rbind(base.no2.100_no2, base.no2.100.coefs[[i]][8,])
  }

write.csv(base.no2.100_no2, "base.no2.100_no2.csv")  

  
base.no2.100_COPD<-data.frame()
  
for (i in 1:6){
    base.no2.100_COPD<-rbind(base.no2.100_COPD, base.no2.100.coefs[[i]][2,])
  }
  

write.csv(base.no2.100_COPD, "base.no2.100_COPD.csv")  


base.no2.100_IHD<-data.frame()
  
for (i in 1:6){
    base.no2.100_IHD<-rbind(base.no2.100_IHD, base.no2.100.coefs[[i]][3,])
}


write.csv(base.no2.100_IHD, "base.no2.100_IHD.csv")


######################################################################################plot#######
#################################################################################################


p10<-ggplot(AA, aes(x=Group, y=ANAP1_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("A.1-ANAP")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
  annotate("rect", xmin = 1, xmax = 3, ymin = 350, ymax =350, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin = 250, ymax =350, alpha=1,colour = "black")+
  annotate("rect", xmin = 3, xmax = 3, ymin = 250, ymax =350, alpha=1,colour = "black")+
  geom_text(aes(x=2, y=400, label="***", size=8))


p11<-ggplot(AA, aes(x=Group, y=ANAP2_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("B.2-ANAP")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
  annotate("rect", xmin = 1, xmax = 2, ymin = 200, ymax =200, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin=150 , ymax =200, alpha=1,colour = "black")+
  annotate("rect", xmin = 2, xmax = 2, ymin = 150, ymax =200, alpha=1,colour = "black")+
  geom_text(aes(x=1.5, y=250, label="*", size=8))+
  annotate("rect", xmin = 1, xmax = 3, ymin = 550, ymax =550, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin = 400, ymax =550, alpha=1,colour = "black")+
  annotate("rect", xmin = 3, xmax = 3, ymin = 400, ymax =550, alpha=1,colour = "black")+
  geom_text(aes(x=2, y=600, label="***", size=8))


p12<-ggplot(AA, aes(x=Group, y=AFLU2_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000),labels=c(0.001, 0.01,0.1,1,10,100,1000),limits=c(0.001,1000))+ 
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("C.2-AFLU")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)


p13<-ggplot(AA, aes(x=Group, y=APHE9_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100,1000,10000),labels=c(0.001, 0.01,0.1,1,10,100,1000,10000),limits=c(0.001,10000))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("D.9-APHE")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)


p14<-ggplot(AA, aes(x=Group, y=APYR1_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.001, 0.01,0.1,1,10,100),labels=c(0.001, 0.01,0.1,1,10,100),limits=c(0.001,100))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("E.1-APYR")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
  annotate("rect", xmin = 1, xmax = 3, ymin = 10, ymax =10, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin = 7, ymax =10, alpha=1,colour = "black")+
  annotate("rect", xmin = 3, xmax = 3, ymin = 7, ymax =10, alpha=1,colour = "black")+
  geom_text(aes(x=2, y=13, label="***", size=8))

p15<-ggplot(AA, aes(x=Group, y=TAPAHs_Cr, fill=Group))+geom_boxplot()+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000,10000),labels=c(0.01,0.1,1,10,100,1000, 10000),limits=c(0.01,10000))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none")+ggtitle("F.\u03A3TAPAHs")+
  stat_summary(fun= mean, geom="point", color="red", shape=15,size=2)+
  annotate("rect", xmin = 1, xmax = 2, ymin = 1000, ymax =1000, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin = 800, ymax =1000, alpha=1,colour = "black")+
  annotate("rect", xmin = 2, xmax = 2, ymin = 800, ymax =1000, alpha=1,colour = "black")+
  geom_text(aes(x=1.5, y=1100, label="**", size=8))+
  annotate("rect", xmin = 1, xmax = 3, ymin = 3500, ymax =3500, alpha=1,colour = "black")+
  annotate("rect", xmin = 1, xmax = 1, ymin = 2500, ymax =3500, alpha=1,colour = "black")+
  annotate("rect", xmin = 3, xmax = 3, ymin = 2500, ymax =3500, alpha=1,colour = "black")+
  geom_text(aes(x=2, y=4000, label="***", size=8))




###function to extract the legend
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
shared_legend <- extract_legend(p9)

grid.arrange(arrangeGrob(p10, p11,p12, p13, p14, p15,nrow=2, 
                         left=textGrob("Urinary concentrations (\u03BCg/g creatinine)", 
                                       gp = gpar(fontsize = 12, fontface = 'bold'),rot = 90, vjust = 1)), shared_legend,heights = c(10, 1))


#######################################################################################################
######################################add diesel data#################################################
######################################################################################################

os2_LAIE_VKM_polyline <- read_csv("C:/Users/zy125/Box Sync/Postdoc/os2/os2_LAIE_VKM_polyline.csv")

n<-str_remove_all(colnames(os2_LAIE_VKM_polyline), "join_")

n[1]<-"id1"

colnames(os2_LAIE_VKM_polyline)<-n


Fulldata_polyline<-merge(Fulldata, os2_LAIE_VKM_polyline, by="id1", all.x=TRUE)

############################################################################################################

baseline_vkm=subset(Fulldata_polyline,Post==0)

BB=data.frame()
BB[1:238,1]=baseline_vkm$ANAP1_Cr
BB[1:238,2]=baseline_vkm$ANAP2_Cr
BB[1:238,3]=baseline_vkm$AFLU2_Cr
BB[1:238,4]=baseline_vkm$APHE9_Cr
BB[1:238,5]=baseline_vkm$APYR1_Cr
BB[1:238,6]=baseline_vkm$TAPAHs_Cr
BB[1:238,7]=baseline_vkm$COPD
BB[1:238,8]=baseline_vkm$IHD
BB[1:238,9]=baseline_vkm$Age
BB[1:238,10]=baseline_vkm$Male
BB[1:238,11]=baseline_vkm$BMI
BB[1:238,12]=baseline_vkm[,1]
BB[1:238,13]=baseline_vkm$AADTDcar13 
BB[1:238,14]=baseline_vkm$VKMDCar13 



BB=na.omit(BB)
colnames(BB)<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id", "AADTDcar13","VKMDCar13")
#################################################################################################
################################################################################################
#######################AADTD####################################################################

aadt.coefs<-list()


for (i in 1:6) {
  
  lmerD=lmer(log10(BB[,i])~log10(BB[,13])+BB[,7]+BB[,8]+BB[,9]+BB[,10]+BB[,11]+(1|BB[,12]))
  aadt.coefs[[i]]<-data.frame(coef(summary(lmerD)))# extract coefficients
  aadt.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(aadt.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  aadt.coefs[[i]]$LCI=aadt.coefs[[i]]$Estimate-1.96*aadt.coefs[[i]]$Std..Error
  aadt.coefs[[i]]$UCI=aadt.coefs[[i]]$Estimate+1.96*aadt.coefs[[i]]$Std..Error
  
  aadt.coefs[[i]]$FC=10^(aadt.coefs[[i]]$Estimate)
  aadt.coefs[[i]]$LFC=10^(aadt.coefs[[i]]$LCI)
  aadt.coefs[[i]]$UFC=10^(aadt.coefs[[i]]$UCI)
  aadt.coefs[[i]]$PC=((10^(aadt.coefs[[i]]$Estimate))-1)*100
  aadt.coefs[[i]]$LPC=((10^(aadt.coefs[[i]]$LCI))-1)*100
  aadt.coefs[[i]]$UPC=((10^(aadt.coefs[[i]]$UCI))-1)*100
}

names(aadt.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(aadt.coefs[[i]], paste0("ModelD_", bio_names[i],".csv"))
}



aadt<-data.frame()

for (i in 1:6){
  aadt<-rbind(aadt, aadt.coefs[[i]][2,])
}


#################################################################################################
################################################################################################
#######################VKMD####################################################################

vkmd.coefs<-list()


for (i in 1:6) {
  
  lmerE=lmer(log10(BB[,i])~BB[,14]+BB[,7]+BB[,8]+BB[,9]+BB[,10]+BB[,11]+(1|BB[,12]))
  vkmd.coefs[[i]]<-data.frame(coef(summary(lmerE)))# extract coefficients
  vkmd.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(vkmd.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  vkmd.coefs[[i]]$LCI=vkmd.coefs[[i]]$Estimate-1.96*vkmd.coefs[[i]]$Std..Error
  vkmd.coefs[[i]]$UCI=vkmd.coefs[[i]]$Estimate+1.96*vkmd.coefs[[i]]$Std..Error
  
  vkmd.coefs[[i]]$FC=10^(vkmd.coefs[[i]]$Estimate)
  vkmd.coefs[[i]]$LFC=10^(vkmd.coefs[[i]]$LCI)
  vkmd.coefs[[i]]$UFC=10^(vkmd.coefs[[i]]$UCI)
  vkmd.coefs[[i]]$PC=((10^(vkmd.coefs[[i]]$Estimate))-1)*100
  vkmd.coefs[[i]]$LPC=((10^(vkmd.coefs[[i]]$LCI))-1)*100
  vkmd.coefs[[i]]$UPC=((10^(vkmd.coefs[[i]]$UCI))-1)*100
}

names(vkmd.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(vkmd.coefs[[i]], paste0("ModelD_", bio_names[i],".csv"))
}



vkmd<-data.frame()

for (i in 1:6){
  vkmd<-rbind(vkmd, vkmd.coefs[[i]][2,])
}

####################################################################################################
##############################################

os2_LAIE_VKM_polygon <- read_csv("C:/Users/zy125/Box Sync/Postdoc/os2/os2_LAIE_VKM_polygon.csv")



colnames(os2_LAIE_VKM_polygon)[1]<-"id1"


Fulldata_polygon<-merge(Fulldata, os2_LAIE_VKM_polygon, by="id1", all.x=TRUE)



baseline_gon=subset(Fulldata_polygon,Post==0)

CC=data.frame()
CC[1:238,1]=baseline_gon$ANAP1_Cr
CC[1:238,2]=baseline_gon$ANAP2_Cr
CC[1:238,3]=baseline_gon$AFLU2_Cr
CC[1:238,4]=baseline_gon$APHE9_Cr
CC[1:238,5]=baseline_gon$APYR1_Cr
CC[1:238,6]=baseline_gon$TAPAHs_Cr
CC[1:238,7]=baseline_gon$COPD
CC[1:238,8]=baseline_gon$IHD
CC[1:238,9]=baseline_gon$Age
CC[1:238,10]=baseline_gon$Male
CC[1:238,11]=baseline_gon$BMI
CC[1:238,12]=baseline_gon[,1]
CC[1:238,13]=baseline_gon$VKMDCar13
CC[1:238,14]=baseline_gon$VKMBus13 



CC=na.omit(CC)
colnames(CC)<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id", "VKMDCar13","VKMBus13")

####################################################################################################
##############################################

plgon.coefs<-list()


for (i in 1:6) {
  
  lmerF=lmer(log10(CC[,i])~CC[,13]+CC[,7]+CC[,8]+CC[,9]+CC[,10]+CC[,11]+(1|CC[,12]))
  plgon.coefs[[i]]<-data.frame(coef(summary(lmerF)))# extract coefficients
  plgon.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(plgon.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  plgon.coefs[[i]]$LCI=plgon.coefs[[i]]$Estimate-1.96*plgon.coefs[[i]]$Std..Error
  plgon.coefs[[i]]$UCI=plgon.coefs[[i]]$Estimate+1.96*plgon.coefs[[i]]$Std..Error
  
  plgon.coefs[[i]]$FC=10^(plgon.coefs[[i]]$Estimate)
  plgon.coefs[[i]]$LFC=10^(plgon.coefs[[i]]$LCI)
  plgon.coefs[[i]]$UFC=10^(plgon.coefs[[i]]$UCI)
  plgon.coefs[[i]]$PC=((10^(plgon.coefs[[i]]$Estimate))-1)*100
  plgon.coefs[[i]]$LPC=((10^(plgon.coefs[[i]]$LCI))-1)*100
  plgon.coefs[[i]]$UPC=((10^(plgon.coefs[[i]]$UCI))-1)*100
}

names(plgon.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(plgon.coefs[[i]], paste0("ModelF_", bio_names[i],".csv"))
}



plgon<-data.frame()

for (i in 1:6){
  plgon<-rbind(plgon, plgon.coefs[[i]][2,])
}

################################################################################################################
plgon.bus.coefs<-list()


for (i in 1:6) {
  
  lmerG=lmer(log10(CC[,i])~log10(CC[,13])+CC[,7]+CC[,8]+CC[,9]+CC[,10]+CC[,11]+(1|CC[,12]))
  plgon.bus.coefs[[i]]<-data.frame(coef(summary(lmerG)))# extract coefficients
  plgon.bus.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(plgon.bus.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  plgon.bus.coefs[[i]]$LCI=plgon.bus.coefs[[i]]$Estimate-1.96*plgon.bus.coefs[[i]]$Std..Error
  plgon.bus.coefs[[i]]$UCI=plgon.bus.coefs[[i]]$Estimate+1.96*plgon.bus.coefs[[i]]$Std..Error
  
  plgon.bus.coefs[[i]]$FC=10^(plgon.bus.coefs[[i]]$Estimate)
  plgon.bus.coefs[[i]]$LFC=10^(plgon.bus.coefs[[i]]$LCI)
  plgon.bus.coefs[[i]]$UFC=10^(plgon.bus.coefs[[i]]$UCI)
  plgon.bus.coefs[[i]]$PC=((10^(plgon.bus.coefs[[i]]$Estimate))-1)*100
  plgon.bus.coefs[[i]]$LPC=((10^(plgon.bus.coefs[[i]]$LCI))-1)*100
  plgon.bus.coefs[[i]]$UPC=((10^(plgon.bus.coefs[[i]]$UCI))-1)*100
}

names(plgon.bus.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(plgon.bus.coefs[[i]], paste0("ModelF_", bio_names[i],".csv"))
}



plgon.bus<-data.frame()

for (i in 1:6){
  plgon.bus<-rbind(plgon.bus, plgon.coefs[[i]][2,])
}



#######################################################################################################



os2_LAIE_VKM_polyline_major <- read_csv("C:/Users/zy125/Box Sync/Postdoc/os2/os2_LAIE_VKM_polyline_major.csv")

m<-str_remove_all(colnames(os2_LAIE_VKM_polyline_major), "join_")

m[1]<-"id1"

colnames(os2_LAIE_VKM_polyline_major)<-m


Fulldata_polyline_major<-merge(Fulldata, os2_LAIE_VKM_polyline_major, by="id1", all.x=TRUE)

############################################################################################################

baseline_vkm_major=subset(Fulldata_polyline_major,Post==0)

DD=data.frame()
DD[1:238,1]=baseline_vkm_major$ANAP1_Cr
DD[1:238,2]=baseline_vkm_major$ANAP2_Cr
DD[1:238,3]=baseline_vkm_major$AFLU2_Cr
DD[1:238,4]=baseline_vkm_major$APHE9_Cr
DD[1:238,5]=baseline_vkm_major$APYR1_Cr
DD[1:238,6]=baseline_vkm_major$TAPAHs_Cr
DD[1:238,7]=baseline_vkm_major$COPD
DD[1:238,8]=baseline_vkm_major$IHD
DD[1:238,9]=baseline_vkm_major$Age
DD[1:238,10]=baseline_vkm_major$Male
DD[1:238,11]=baseline_vkm_major$BMI
DD[1:238,12]=baseline_vkm_major[,1]
DD[1:238,13]=baseline_vkm_major$AADTDcar13 
DD[1:238,14]=baseline_vkm_major$VKMDCar13 
DD[1:238,15]=baseline_vkm_major$distance
DD[1:238,16]=baseline_vkm_major$VKMDCar13/baseline_vkm_major$distance


DD=na.omit(DD)
colnames(DD)<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id", "AADTDcar13","VKMDCar13", "d.NN.major","vkm.inv.dis"
                )


#################################################################################################
################################################################################################
#######################AADTD####################################################################

aadt.major.coefs<-list()


for (i in 1:6) {
  
  lmerH=lmer(log10(DD[,i])~DD[,13]+DD[,7]+DD[,8]+DD[,9]+DD[,10]+DD[,11]+(1|DD[,12]))
  aadt.major.coefs[[i]]<-data.frame(coef(summary(lmerH)))# extract coefficients
  aadt.major.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(aadt.major.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  aadt.major.coefs[[i]]$LCI=aadt.major.coefs[[i]]$Estimate-1.96*aadt.major.coefs[[i]]$Std..Error
  aadt.major.coefs[[i]]$UCI=aadt.major.coefs[[i]]$Estimate+1.96*aadt.major.coefs[[i]]$Std..Error
  
  aadt.major.coefs[[i]]$FC=10^(aadt.major.coefs[[i]]$Estimate)
  aadt.major.coefs[[i]]$LFC=10^(aadt.major.coefs[[i]]$LCI)
  aadt.major.coefs[[i]]$UFC=10^(aadt.major.coefs[[i]]$UCI)
  aadt.major.coefs[[i]]$PC=((10^(aadt.major.coefs[[i]]$Estimate))-1)*100
  aadt.major.coefs[[i]]$LPC=((10^(aadt.major.coefs[[i]]$LCI))-1)*100
  aadt.major.coefs[[i]]$UPC=((10^(aadt.major.coefs[[i]]$UCI))-1)*100
}

names(aadt.major.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(aadt.major.coefs[[i]], paste0("ModelH_", bio_names[i],".csv"))
}



aadt.major<-data.frame()

for (i in 1:6){
  aadt.major<-rbind(aadt.major, aadt.major.coefs[[i]][2,])
}

#################################################################################################
################################################################################################
#######################VKMD####################################################################

vkmd.major.coefs<-list()


for (i in 1:6) {
  
  lmerI=lmer(log10(DD[,i])~DD[,16]+DD[,7]+DD[,8]+DD[,9]+DD[,10]+DD[,11]+(1|DD[,12]))
  vkmd.major.coefs[[i]]<-data.frame(coef(summary(lmerI)))# extract coefficients
  vkmd.major.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(vkmd.major.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  vkmd.major.coefs[[i]]$LCI=vkmd.major.coefs[[i]]$Estimate-1.96*vkmd.major.coefs[[i]]$Std..Error
  vkmd.major.coefs[[i]]$UCI=vkmd.major.coefs[[i]]$Estimate+1.96*vkmd.major.coefs[[i]]$Std..Error
  
  vkmd.major.coefs[[i]]$FC=10^(vkmd.major.coefs[[i]]$Estimate)
  vkmd.major.coefs[[i]]$LFC=10^(vkmd.major.coefs[[i]]$LCI)
  vkmd.major.coefs[[i]]$UFC=10^(vkmd.major.coefs[[i]]$UCI)
  vkmd.major.coefs[[i]]$PC=((10^(vkmd.major.coefs[[i]]$Estimate))-1)*100
  vkmd.major.coefs[[i]]$LPC=((10^(vkmd.major.coefs[[i]]$LCI))-1)*100
  vkmd.major.coefs[[i]]$UPC=((10^(vkmd.major.coefs[[i]]$UCI))-1)*100
}

names(vkmd.major.coefs)<-c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

for (i in 1:6){
  write.csv(vkmd.major.coefs[[i]], paste0("ModelI_", bio_names[i],".csv"))
}



vkmd.major<-data.frame()

for (i in 1:6){
  vkmd.major<-rbind(vkmd.major, vkmd.major.coefs[[i]][2,])
}




##############################################################################################################
#####################################9/11/2020 after discussing with yan, add some new data and analysis#####
##################################################################################################################


os2_LAEI_VKM_polyline_buffer100m <- read_csv("C:/Users/zy125/Box Sync/Postdoc/os2/os2_LAEI_VKM_polyline_buffer100m.csv")


colnames(os2_LAEI_VKM_polyline_buffer100m)

os2_buffer100m<-os2_LAEI_VKM_polyline_buffer100m[, c(1,51:76, 103, 109)]

####diesel
os2_buffer100m$aadt.tot.d<-os2_buffer100m$AADTDcar13+os2_buffer100m$AADTDLgv13+os2_buffer100m$AADTCoac_1+os2_buffer100m$AADTRigi_1+os2_buffer100m$AADTArti_1
os2_buffer100m$aadt.tot.d.length<-os2_buffer100m$aadt.tot.d*os2_buffer100m$length_indi

os2_buffer100m$vkm.tot.d<-os2_buffer100m$VKMDCar13+os2_buffer100m$VKMDLgv13+os2_buffer100m$VKMCoach13+os2_buffer100m$VKMRigid13+os2_buffer100m$VKMArtic13

#####petrol

os2_buffer100m$aadt.tot.p<-os2_buffer100m$AADTPcar13+os2_buffer100m$AADTPLgv13
os2_buffer100m$aadt.tot.p.length<-os2_buffer100m$aadt.tot.p*os2_buffer100m$length_indi

os2_buffer100m$vkm.tot.p<-os2_buffer100m$VKMPCar13+os2_buffer100m$VKMPLgv13
#####six variables

os2_buffer100m$AADTMoto_1_len<-os2_buffer100m$AADTMoto_1*os2_buffer100m$length_indi
os2_buffer100m$AADTTaxi13_len<-os2_buffer100m$AADTTaxi13*os2_buffer100m$length_indi
os2_buffer100m$AADTLtBu_1_len<-os2_buffer100m$AADTLtBu_1*os2_buffer100m$length_indi
os2_buffer100m$AADTCoac_1_len<-os2_buffer100m$AADTCoac_1*os2_buffer100m$length_indi
os2_buffer100m$AADTRigi_1_len<-os2_buffer100m$AADTRigi_1*os2_buffer100m$length_indi
os2_buffer100m$AADTArti_1_len<-os2_buffer100m$AADTArti_1*os2_buffer100m$length_indi
os2_buffer100m$AADTTOTA_1_len<-os2_buffer100m$AADTTOTA_1*os2_buffer100m$length_indi







colnames(os2_buffer100m)[28]<-"id"

os2_b100m<-os2_buffer100m %>% group_by(id) %>% 
  summarize(aadt.d=sum(aadt.tot.d),aadt.p=sum(aadt.tot.p), 
            aadt.length.d=sum(aadt.tot.d.length),aadt.length.p=sum(aadt.tot.p.length),
            aadt.taxi=sum(AADTTaxi13), aadt.moto=sum(AADTMoto_1),
            aadt.bus=sum(AADTLtBu_1), aadt.d.heavy=sum(AADTDLgv13+AADTCoac_1+AADTRigi_1+AADTArti_1),
            aadt.e=sum(AADTECar13+AADT_ELg_1),aadt.total=sum(AADTTOTA_1),
            vkm.d=sum(vkm.tot.d), vkm.p=sum(vkm.tot.p),vkm.bus=sum(VKMBus13),
            vkm.taxi=sum(VKMTaxi13), vkm.d.heavy=sum(VKMDLgv13+VKMCoach13+VKMRigid13+VKMArtic13),
            vkm.e=sum(VKMECar13+VKMELgv13),vkm.total=sum(VKMTOTAL13)
            )


colnames(DD)[12]<-"id"

DD<-merge(DD, os2_b100m, all.x = TRUE)

DD[is.na(DD)]<-0

DD<-DD[, c(2,3,4,5,6,7,8,9,10,11,12,1,13:33)]

colnames(DD)[1:16]<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id", "nn.aadt.d","nn.vkm.d", "d.NN.major","vkm.inv.dis")





os2_b100m_len<-os2_buffer100m %>% group_by(id) %>% 
  summarize(aadt.moto_len=sum(AADTMoto_1_len), aadt.taxi_len=sum(AADTTaxi13_len),
            aadt.bus_len=sum(AADTLtBu_1_len),aadt.coach_len=sum(AADTCoac_1_len),
            aadt.rigid_len=sum(AADTRigi_1_len), aadt.artic_len=sum(AADTArti_1_len),
            aadt.total_len=sum(AADTTOTA_1_len)
             )


EE<-merge(DD, os2_b100m_len, all.x = TRUE)

EE[is.na(EE)]<-0

EE<-EE[, c(2,3,4,5,6,7,8,9,10,11,12,1,13:40)]

#################################################################################################
################################################################################################
#######################AADTE####################################################################
buffer.var<-list()
result<-data.frame()

for(j in 17:33){

  indi.buffer.coefs<-list()


   for (i in 1:6) {
  
  lmerY=lmer(log10(DD[,i])~DD[,j]+DD[,7]+DD[,8]+DD[,9]+DD[,10]+DD[,11]+(1|DD[,12]))
  indi.buffer.coefs[[i]]<-data.frame(coef(summary(lmerY)))# extract coefficients
  indi.buffer.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(indi.buffer.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
  indi.buffer.coefs[[i]]$LCI=indi.buffer.coefs[[i]]$Estimate-1.96*indi.buffer.coefs[[i]]$Std..Error
  indi.buffer.coefs[[i]]$UCI=indi.buffer.coefs[[i]]$Estimate+1.96*indi.buffer.coefs[[i]]$Std..Error
  
  indi.buffer.coefs[[i]]$FC=10^(indi.buffer.coefs[[i]]$Estimate)
  indi.buffer.coefs[[i]]$LFC=10^(indi.buffer.coefs[[i]]$LCI)
  indi.buffer.coefs[[i]]$UFC=10^(indi.buffer.coefs[[i]]$UCI)
  indi.buffer.coefs[[i]]$PC=((10^(indi.buffer.coefs[[i]]$Estimate))-1)*100
  indi.buffer.coefs[[i]]$LPC=((10^(indi.buffer.coefs[[i]]$LCI))-1)*100
  indi.buffer.coefs[[i]]$UPC=((10^(indi.buffer.coefs[[i]]$UCI))-1)*100
    }


   indi.buffer<-data.frame()

for (p in 1:6){

    indi.buffer<-rbind(indi.buffer, indi.buffer.coefs[[p]][2,])
   
    }
 buffer.var[[j-16]]<-indi.buffer
 names(buffer.var[[j-16]])<-colnames(DD[j])

 for (k in 1:6){
  
   result[j-16, 1+(k-1)*4]=buffer.var[[j-16]][k, 9]
   result[j-16, 2+(k-1)*4]=buffer.var[[j-16]][k, 10]
   result[j-16, 3+(k-1)*4]=buffer.var[[j-16]][k, 11]
   result[j-16, 4+(k-1)*4]=buffer.var[[j-16]][k, 6]
  
}



}

nam<-data.frame()

for (u in 1:6){
  nam<-append(nam, c(paste0("FC","(",bio_names[u],")"),
                     paste0("LFC","(",bio_names[u],")"),paste0("UFC","(",bio_names[u],")"),
                     paste0("p.value","(",bio_names[u],")")))
  
}
colnames(result)<-nam

result$nam<-colnames(DD[17:33])


for (i in 1:17){
  write.csv(buffer.var[[i]], paste0("Buffer_Model_",colnames(DD[i+16]),".csv"))
}


################################################################################################
#####six length variables

FF=data.frame()

FF<-EE[, c(1:12,34:40)]
FF[,13]<-standardize(EE[,34], centerFun=mean, scaleFun = sd)
FF[,14]<-standardize(EE[,35], centerFun=mean, scaleFun = sd)
FF[,15]<-standardize(EE[,36], centerFun=mean, scaleFun = sd)
FF[,16]<-standardize(EE[,37], centerFun=mean, scaleFun = sd)
FF[,17]<-standardize(EE[,38], centerFun=mean, scaleFun = sd)
FF[,18]<-standardize(EE[,39], centerFun=mean, scaleFun = sd)
FF[,19]<-standardize(EE[,40], centerFun=mean, scaleFun = sd)



len.var<-list()
result_len<-data.frame()

for(j in 13:19){
  
  indi.len.coefs<-list()
  
  
  for (i in 1:6) {
    
    lmerZ=lmer(log10(FF[,i])~FF[,j]+FF[,7]+FF[,8]+FF[,9]+FF[,10]+FF[,11]+(1|FF[,12]))
    indi.len.coefs[[i]]<-data.frame(coef(summary(lmerZ)))# extract coefficients
    indi.len.coefs[[i]]$p.z <- 2 * (1 - pnorm(abs(indi.len.coefs[[i]]$t.value)))# use normal distribution to approximate p-value
    indi.len.coefs[[i]]$LCI=indi.len.coefs[[i]]$Estimate-1.96*indi.len.coefs[[i]]$Std..Error
    indi.len.coefs[[i]]$UCI=indi.len.coefs[[i]]$Estimate+1.96*indi.len.coefs[[i]]$Std..Error
    
    indi.len.coefs[[i]]$FC=10^(indi.len.coefs[[i]]$Estimate)
    indi.len.coefs[[i]]$LFC=10^(indi.len.coefs[[i]]$LCI)
    indi.len.coefs[[i]]$UFC=10^(indi.len.coefs[[i]]$UCI)
    indi.len.coefs[[i]]$PC=((10^(indi.len.coefs[[i]]$Estimate))-1)*100
    indi.len.coefs[[i]]$LPC=((10^(indi.len.coefs[[i]]$LCI))-1)*100
    indi.len.coefs[[i]]$UPC=((10^(indi.len.coefs[[i]]$UCI))-1)*100
  }
  
  
  indi.len<-data.frame()
  
  for (p in 1:6){
    
    indi.len<-rbind(indi.len, indi.len.coefs[[p]][2,])
    
  }
  len.var[[j-12]]<-indi.len
  names(len.var[[j-12]])<-colnames(FF[j])
  
  for (k in 1:6){
    
    result_len[j-12, 1+(k-1)*4]=len.var[[j-12]][k, 9]
    result_len[j-12, 2+(k-1)*4]=len.var[[j-12]][k, 10]
    result_len[j-12, 3+(k-1)*4]=len.var[[j-12]][k, 11]
    result_len[j-12, 4+(k-1)*4]=len.var[[j-12]][k, 6]
    
  }
  
  
  
}

nam.len<-data.frame()

for (u in 1:6){
  nam.len<-append(nam.len, c(paste0("FC","(",bio_names[u],")"),
                     paste0("LFC","(",bio_names[u],")"),paste0("UFC","(",bio_names[u],")"),
                     paste0("p.value","(",bio_names[u],")")))
  
}
colnames(result_len)<-nam.len

result_len$nam.len<-colnames(FF[13:19])


for (i in 1:6){
  write.csv(len.var[[i]], paste0("Len_Model_",colnames(FF[i+12]),".csv"))
}


###############################################################

library(readxl)
length <- read_excel("length.xlsx")
colnames(length)[7]<-"significance"

length$significance<-as.logical(length$significance)

length$biomarker<-factor(length$biomarker, c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs"))

ggplot(length, aes(y=FC, x=biomarker,color=significance))+
  geom_pointrange(aes(ymax=UFC, ymin=LFC, shape=type), position = pd, size=0.8)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5),linetype="dotted", size=1)+
  geom_hline(yintercept = 1)+theme_classic()+theme(legend.position="bottom")+ 
  guides(guide_legend(nrow = 1))+scale_colour_manual(values=c("black", "red"))+labs(y="Fold Change")










