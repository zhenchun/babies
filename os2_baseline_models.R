

#Packages
library(dlnm)
library(lme4)
library(lmerTest)
require(ggplot2)

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
AA=na.omit(AA)
colnames(AA)<-c("ANAP1_Cr", "ANAP2_Cr","AFLU2_Cr","APHE9_Cr", "APYR1_Cr","TAPAHs_Cr",
                "COPD", "IHD","Age","Male", "BMI","id","all_roads_length_100m",
                "major_roads_length_100m","minor_roads_length_100m","Q1", "Q2","Q3","Q4", 
                "Near100","Far100")

bio_names<- c("ANAP1", "ANAP2","AFLU2","APHE9", "APYR1","TAPAHs")

###section A#########################################################
####test the roadlength within 100m buffers effect on concentrations##

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


###section B#########################################################
####test the in or outside 100m of a major road effect on concentrations##

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



###section C#########################################################
####test the roadlength within 100m buffers effect on concentrations##


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


road_length<-rbind(major.length,minor.length,all.length)




ggplot()


pd <- position_dodge(width = 1)

star.label<-data.frame(Group="AFLU2", Value=2.5)

p1<-ggplot(road_length, aes(y=PC*10, x=names))+
       geom_pointrange(aes(ymax=UPC*10, ymin=LPC*10, color=type), position = pd, size=0.8)+
       geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5),linetype="dotted", size=1)+
       geom_hline(yintercept = 0)+theme_classic()+
       scale_color_manual(values=c("black","deepskyblue","coral"),labels = c("Major Roads","Minor Roads","All Roads"))+
       labs(x="biomarker names", y="Percentage Change (%, 95%CI)")+
       scale_x_discrete(labels=c("1-ANAP", "2-ANAP", "2-AFLU", "9-APHE", "1-APYR", "TAPAHs"))+
       geom_text(x=2.67, y=7.75, label="*", size=5)+
       theme(legend.title=element_blank(),legend.position="bottom",  
         axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.text=element_text(size=14), 
         axis.title.x=element_blank(), axis.title.y = element_text(size=16),legend.key.size=unit(4, "points"))




road_length$names<-factor(road_length$names, levels=c("major.road","minor.road","all.road"))








