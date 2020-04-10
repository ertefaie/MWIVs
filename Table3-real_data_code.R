# prison vs no prison code. 

rm(list=ls())
library(ivpack);library(MASS);library(ivmodel);library(Matrix);library(dummies);library(foreign)
library(lme4)
#dat<-read.csv("C:\\SentencingProject\\upenn\\analysis_master_file.csv", header=T)
#dat<-read.dta("C:\\SentencingProject\\upenn\\analysis_master_file2.dta")

#save(dat,file="C:\\SentencingProject\\upenn\\Sentencing_data2.RData")


Michigan.Sentencing.CIs<-function(delta) {
  
  load("C:\\SentencingProject\\upenn\\Sentencing_data.RData")
  
  
  
  ############# DATA PREPARATION ################
  dat<-dat[dat$race_cat=="Black"|dat$race_cat=="White",]
  dat<-dat[as.numeric(dat$educ_cat)!=2,]
  
  
  set.seed(301)
 
  n_obs<-nrow(dat)
  #########################
  length(unique(dat$r_judge_id))
  Z<-dummy(dat$r_judge_id)
  #Z<-scale(Z)
  temp<-Z[1,]%*%solve(t(Z)%*%Z/n_obs)%*%Z[1,]/sqrt(n_obs)
  
  
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
  cat.county<-dummy(dat$r_county_new)[,-1]
  dat$cat.county<-cat.county
  totcases<-rep(1,length(dat$r_county_new))
  
  temp<-data.frame(cbind(dat$r_judge_id,dat$r_county_new))
  colnames(temp)<-c("r_judge_id","r_county_new")
  temp<-temp[order(temp[,2]),]
  temp$diff<-c(1,diff(temp[,2]))
  head(temp)
  #NotExcluded<-temp[as.integer(temp[,3])==FALSE|temp[,3]==0,1]
  Excluded<-temp[temp[,3]!=0,1]
  
  Z<-Z[,colnames(Z)%in%Excluded==FALSE]
  
  temp<-cat.county
  
  covariates<-cbind(dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,age=dat$age_cat_ssd, sex=dat$r_sex,dat$marital_ds,dat$marital_single,dat$marital_marr,race=dat$race_b,edu=dat$educ_cat, priormis=dat$priormis_dummy,priorfel=dat$priorfel_dummy,drug=dat$r_drug,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,mental=dat$r_mental_h,narrest=dat$arrest_cat,dat$r_senyear3,dat$r_senyear4,dat$r_senyear5,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
  
  #covariates<-cbind(dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,age=dat$age_demean,age2=dat$age_demeansq,sex=dat$r_sex,dat$marital_ds,dat$marital_single,dat$marital_marr,race=dat$race_b,edu=dat$educ_cat, priormis=dat$priormis_dummy,priorfel=dat$priorfelony_1,drug=dat$r_drug,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,mental=dat$r_mental_h,narrest=dat$arrest_cat,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
  
  txt2<-as.numeric(dat$sent_type=="Prison")  # treatment indicator
  txt<-scale(txt2,scale=FALSE)
  X<-txt
  y<-scale(log(dat$mean_earn24_ssd+1),scale=FALSE)  # outcome
  
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
  sdelt<- delta# -0.0 
  
  #######removing extreme judges#####
  
  harshness<-t(coef(lm(X ~Z-1)))
  q<-quantile(harshness,prob=c(.1,.9),na.rm=TRUE)
  
  remove10<-na.omit(colnames(Z)[harshness<q[1]])
  remove90<-na.omit(colnames(Z)[harshness>q[2]])
  remove<-c(remove10,remove90)
  
  
  judge.index<-substr(remove,start=11,stop=15)
  
  dat<-dat[dat$r_judge_id%in%judge.index==FALSE,]
  Z<-dummy(dat$r_judge_id)
  temp<-data.frame(cbind(dat$r_judge_id,dat$r_county_new))
  colnames(temp)<-c("r_judge_id","r_county_new")
  temp<-temp[order(temp[,2]),]
  temp$diff<-c(1,diff(temp[,2]))
  head(temp)
  #NotExcluded<-temp[as.integer(temp[,3])==FALSE|temp[,3]==0,1]
  Excluded<-temp[temp[,3]!=0,1]
  judge.index.Z<-substr(colnames(Z),start=11,stop=15)
  
  Z<-Z[,judge.index.Z%in%Excluded==FALSE]
  
  covariates<-cbind(dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,age=dat$age_cat_ssd, sex=dat$r_sex,dat$marital_ds,dat$marital_single,dat$marital_marr,race=dat$race_b,edu=dat$educ_cat, priormis=dat$priormis_dummy,priorfel=dat$priorfel_dummy,drug=dat$r_drug,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,mental=dat$r_mental_h,narrest=dat$arrest_cat,dat$r_senyear3,dat$r_senyear4,dat$r_senyear5,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
  
  #covariates<-cbind(dat$cat.county,dat$crime_cont_sub,dat$crime_person,dat$crime_property,dat$crime_pub_safe,dat$crime_pub_order,age=dat$age_demean,age2=dat$age_demeansq,sex=dat$r_sex,dat$marital_ds,dat$marital_single,dat$marital_marr,race=dat$race_b,edu=dat$educ_cat, priormis=dat$priormis_dummy,priorfel=dat$priorfelony_1,drug=dat$r_drug,alcohol=dat$r_alcohol,employ23=dat$pre_employ23,employ12=dat$pre_employ12, wages23=dat$pre_wages23,wages12=dat$pre_wages12,mental=dat$r_mental_h,narrest=dat$arrest_cat,dat$marijuana_use,dat$opoids_use,dat$other_drug_use)
  
  table(dat$marital_cat)
  
  txt2<-as.numeric(dat$sent_type=="Prison")
  txt<-scale(txt2,scale=FALSE)
  X<-txt
  
  y<-scale(log(dat$mean_earn24_ssd+1),scale=FALSE)
  
  ####Orthogonal PRojections onto space X #######
  
  temp<-cbind(y,txt)
  temp2<- resid(lm(temp~covariates))
  temp3<-resid(lm(Z~covariates))

  y<-temp2[,1]
  txt<-temp2[,2]
  X<-txt<-scale(txt,scale=FALSE)
  Z<-(temp3)#Z<-Z[,-1]
  Z<-apply(Z,2,scale)
  K<-ncol(Z)
  sdelt<- delta# -1.2 
  
  
  gamma3<-gamma2<-t(coef(lm(X ~Z-1)))
  delt.vec<-seq(-0.45,-0.35,by=0.01)
  #delt.vec<-seq(-0.06,0.00,by=0.01)
  #delt.vec<-seq(-2.2,-1.3,by=0.02)
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
  
  cbind(delt.vec, pval.AR,pval.K,pval.J,CI.AR=pval.AR>0.05,CI.K=pval.K>0.05,CI.KJ14=(pval.K>0.04)*(pval.J>0.01) )
  
}


Michigan.Sentencing.CIs(delta=0)


