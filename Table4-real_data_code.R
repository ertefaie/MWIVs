
# prison vs no prison code. 

rm(list=ls())
library(ivpack);library(MASS);library(ivmodel);library(Matrix);library(dummies);library(foreign)
library(lme4)
#dat<-read.csv("C:\\SentencingProject\\upenn\\analysis_master_file.csv", header=T)
#dat<-read.dta("C:\\SentencingProject\\upenn\\analysis_master_file2.dta")

#save(dat,file="C:\\SentencingProject\\upenn\\Sentencing_data2.RData")

# define subgroups in line 26.
Michigan.Sentencing.SubG.CIs<-function(delta) {
  

load("C:\\SentencingProject\\upenn\\Sentencing_data.RData")


table(dat$race_cat)
dat<-dat[dat$race_cat=="Black"|dat$race_cat=="White",]
dat<-dat[as.numeric(dat$educ_cat)!=2,]

table(dat$age_cat_ssd)/dim(dat)[1]

set.seed(301)
################### DEFINE A SUPGROUP  ########## 
note<-"dat$r_sex==Male & dat$race_cat!=White&dat$priorfel_dummy==1&as.numeric(dat$age_cat_ssd)<4"

dat<-dat[dat$r_sex=="Male"&dat$race_cat=="White"&as.numeric(dat$age_cat_ssd)<3,]

################################################
totfel<-rep(1,nrow(dat))
temp<-data.frame(totfel,dat$r_judge_id)
temp<-aggregate(totfel~dat$r_judge_id,temp,FUN=sum)
colnames(temp)[1]<-"r_judge_id"
temp<-temp[temp[2]>4,]
dat<-merge(dat,temp,by="r_judge_id")


n_obs<-nrow(dat)
n_obs
#########################
length(unique(dat$r_judge_id))
Z<-dummy(dat$r_judge_id)

dat$crime_cont_subs<-as.numeric(dat$crime_cat=="Controlled Substance")
dat$crime_person<-as.numeric(dat$crime_cat=="Person")
dat$crime_property<-as.numeric(dat$crime_cat=="Property")
dat$crime_pub_order<-as.numeric(dat$crime_cat=="Public Order")
dat$crime_pub_safe<-as.numeric(dat$crime_cat=="Public Safety")
dat$race_b<-as.numeric(dat$race_cat=="Black")
dat$marital_ds<-as.numeric(dat$marital_cat=="Divorced or Separated")
dat$marital_single<-as.numeric(dat$marital_cat=="Single")
dat$marital_marr<-as.numeric(dat$marital_cat=="Married or Common Law Marriage")
dat$r_senyear3<-as.numeric(dat$r_senyear=="2003")
dat$r_senyear4<-as.numeric(dat$r_senyear=="2004")
dat$r_senyear5<-as.numeric(dat$r_senyear=="2005")
cat.county<-dummy(dat$r_county_new)#[,-1]
dat$cat.county<-cat.county
totcases<-rep(1,length(dat$r_county_new))
dat$cat.county<-dat$cat.county[,apply(dat$cat.county,2,mean)!=0]


temp<-data.frame(cbind(dat$r_judge_id,dat$r_county_new))
colnames(temp)<-c("r_judge_id","r_county_new")
temp<-temp[order(temp[,2]),]

temp$diff<-c(1,diff(temp[,2]))
head(temp)
#NotExcluded<-temp[as.integer(temp[,3])==FALSE|temp[,3]==0,1]
Excluded<-temp[temp[,3]!=0,1]

Z<-Z[,colnames(Z)%in%Excluded==FALSE]

temp<-cat.county

covariates<-cbind(dat$age_cat_ssd,dat$priorfel_dummy, dat$r_drug, dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,dat$marital_ds,dat$marital_single,dat$marital_marr,edu=dat$educ_cat, priormis=dat$priormis_dummy,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,narrest=dat$arrest_cat,dat$r_senyear3,dat$r_senyear4,dat$r_senyear5,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
#covariates<-cbind(dat$priorfel_dummy, dat$r_drug, dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,dat$marital_ds,dat$marital_single,dat$marital_marr,edu=dat$educ_cat, priormis=dat$priormis_dummy,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,narrest=dat$arrest_cat,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)

#covariates<-covariates[,-7]
#covariates<-covariates[,-28]

#covariates<-cbind(dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,age=dat$age_demean,age2=dat$age_demeansq,sex=dat$r_sex,dat$marital_ds,dat$marital_single,dat$marital_marr,race=dat$race_b,edu=dat$educ_cat, priormis=dat$priormis_dummy,priorfel=dat$priorfelony_1,drug=dat$r_drug,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,mental=dat$r_mental_h,narrest=dat$arrest_cat,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
table(as.numeric(dat$educ_cat))

apply(dat$cat.county,2,mean)

table(dat$marital_cat)

sen_type_probs<-aggregate(cbind(Prison=dat$sent_type=="Prison",JailProb=dat$sent_type=="Jail with Probation",Jail=dat$sent_type=="Jail Only", Prob=dat$sent_type=="Probation")~dat$r_judge_id,FUN=mean)
sen_type_probs<-sen_type_probs[order(sen_type_probs$Prison),]

###C:\Users\ashkanertefaie\Dropbox\Dylan\Jeff files\Code-Jeff\Sentencing_data_code

pdf("file.pdf",width=15,height=8,paper='special')
plot(1:150, sen_type_probs$Prison[1:150],ylim=c(0,0.8),ylab="Proportion",xlab="Judges",cex.axis=1.5,cex.lab=1.5)
points(sen_type_probs$JailProb[1:150],pch=3,cex=.8)
points(sen_type_probs$Jail[1:150],pch=4,cex=.8)
points(sen_type_probs$Prob[1:150],pch=17,cex=.8)
legend("topright",c("Prison","Jail+ Probation","Jail","Probation"),pch=c(1,3,4,17),cex=1.5,bty="n")
dev.off()


txt2<-as.numeric(dat$sent_type=="Prison")
txt<-scale(txt2,scale=FALSE)
table(dat$sent_type)
X<-txt

y<-scale(log(dat$mean_earn24_ssd+1),scale=FALSE)

####Orthogonal PRojections onto space X#######

temp<-cbind(y,txt)

temp2<- t(temp)- t(temp)%*%covariates%*%solve(t(covariates)%*% covariates)%*%t(covariates)
temp3<- t(Z)- t(Z)%*%covariates%*%solve(t(covariates)%*% covariates)%*%t(covariates)

y<-temp2[1,]
txt<-temp2[2,]
X<-txt<-scale(txt,scale=FALSE)
Z<-t(temp3)
#Z<-Z[,-1]
Z<-apply(Z,2,scale)


####Orthogonal PRojections onto space X#######

temp<-cbind(y,txt)
temp2<- resid(lm(temp~covariates))
temp3<-resid(lm(Z~covariates))

y<-temp2[,1]
txt<-temp2[,2]
X<-txt<-scale(txt,scale=FALSE)
Z<-(temp3)#Z<-Z[,-1]
Z<-apply(Z,2,scale)
K<-ncol(Z)
sdelt<- delta# -0.0  # sensitivity parameter



n_obs<-nrow(dat)


dat.temp<-data.frame(cbind(y,txt,Z))
###############2SLS###########################
anova(lm(txt~Z))
fm2<-ivreg(y~covariates+ txt | covariates+Z,data= dat.temp )
summary(fm2,vcov=sandwich);

summary(lm(y~covariates+txt,data= dat.temp))
##########################################

gamma3<-gamma2<-t(coef(lm(X ~Z-1)))
delt.vec<-seq(-5,5,by=0.05)
#delt.vec<-seq(-0.06,0.00,by=0.01)
#delt.vec<-seq(-30.3,-20,by=0.05)
pval.AR<-pval.K<-pval.J <-NULL
#PP<-solve(t(Z)%*%Z)%*%t(Z)
#delt.vec<-0.7
for(j in 1:length(delt.vec)){
  delta<- delt.vec[j]
  u<-y-X*delta-Z%*% (gamma2[1,])*sdelt
  #  u<-scale(u,scale=FALSE)
  
  sev<-as.vector((t(u)%*%X)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
  svv<-((t(X)%*%X)/(n_obs-K))-((t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
  see<-as.vector((t(u)%*%u)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u)/(n_obs-K))
  
  
  ##################Kleibergen 2002#######################
  pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-u%*%sev/see)
  Ytilde<- Z %*% pi.tilde
  #P_Ytilda<-Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)
  LM2<-t(u)%*%Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)%*%u/(see + svv[1,1]* sdelt^2 -2*(sdelt)*sev[1])
  ncp<-0
  if(LM2< qchisq(1-0.05,1,ncp= 0,lower.tail=TRUE)){print(c(delta,"KL CI"))}
  pval.K[j]<-1-pf(LM2,1,n_obs-K-1,ncp= 0,lower.tail=TRUE)
  
  ###################Kleibergen 2006####################
  pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-u%*%sev/see)
  D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
  middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(see + svv[1,1]* sdelt^2 -2*(sdelt)*sev[1])
  
  J_KLM<-t(u)%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%(u)
  ncpJ<-0
  
  pval.J[j]<-1-pchisq(J_KLM,K-1,ncp= 0,lower.tail=TRUE)
   
  
  ##################AR#######################
  AR<-t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u/(see + svv[1,1]* sdelt^2 -2*(sdelt)*sev[1])/K
  ncpAR<-0
  AR
  
  if(AR< qf(1-0.05,K,n_obs-K,ncp=0,lower.tail=TRUE)){print(c(delta,"AR CI"))}
  
  pval.AR[j]<-1-pf(AR,K,n_obs-K-1,ncp= 0,lower.tail=TRUE)
  
  print(c(AR,LM2))
  if(j%%10==0) print(j)
}

cbind(pval.AR,pval.K,pval.J,delt.vec,CI.AR=pval.AR>0.05,CI.K=pval.K>0.05,CI.KJ14=(pval.K>0.04)*(pval.J>0.01) )

}
ptm<-proc.time()
Michigan.Sentencing.SubG.CIs(delta=0)
proc.time()-ptm


#fm2<-ivreg(y~ txt | Z,data= dat.temp )
#tail(summary(fm2,vcov=sandwich)$coef,1);
#tail(confint((fm2)),1)

 
#tail(summary(lm(y~covariates+txt,data= dat.temp))$coef,1)
#tail(confint(lm(y~covariates+txt,data= dat.temp)),1)


