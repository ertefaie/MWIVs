
 

#rm(list=ls())
library(ivpack);library(MASS);library(ivmodel);library(Matrix);library(dummies)

#dat<-read.csv("C:\\SentencingProject\\upenn\\analysis_master_file.csv", header=T)
K<-2 # number of instruments
n_obs<-500 # sample size
irep<-1000 # number of datasets
strength<-c(2,rep(0,K-1)) # strength of the IV
sdelt<- 0# sensitivity parameter




mean.true.delta<- 1.0
delt.vec<-seq(-1.2,1.2,by=0.1) # true parameter value.
delt.vec<-seq(0.1,1.3,by=0.05) # true parameter value.
#delt.vec<-seq(0.9,1.0,by=0.05) # true parameter value.
cont<-contK<-contKJ14 <-0
pval.AR<-pval.K<-pval.KJ14 <-pval.KJ32 <-pval.J<-CI.AR<-CI.K<-CI.KJ14<-CI.KJ41 <-CI.KJ32<-matrix(0,nrow=irep,ncol=length(delt.vec))
herto_effect<-NULL
#########################
for( i in 1:irep){

Z1<-rbinom(n_obs ,1,0.5) 
Z2<-rbinom(n_obs ,1,0.5) 
Z3<-rbinom(n_obs ,1,0.5) 
Z<-cbind(Z1,Z2,Z3)
dim(Z);head(Z)

# Scaling is critical. It may fail without scaling
Z<-apply(Z,2,scale)
Z1<-Z[,1]
Z2<-Z[,2]
Z3<-Z[,3]


Sigma<-matrix(0.5,ncol=2,nrow=2)
diag(Sigma)<-1
Sigma[2,2]<-1
mx<-mvrnorm(n_obs,c(0,0),Sigma)

  logitprob<-Z1*strength[1]*1+1.0*Z2*strength[1]-0.0*Z3*strength[1]+mx[,1] 
#  txt2<-txt<-rbinom(n_obs,1,exp(logitprob)/(1+exp(logitprob)))
txt2<-txt<-logitprob
txt<-scale(txt,scale=FALSE)
X<-txt
  
#  true.delta<-rnorm(n_obs, mx[,1],.1)+ mean.true.delta
#  true.delta<- mx[,1]*1.0+ mean.true.delta+rnorm(n_obs)
#  true.delta<- Z1*mx[,1]*1.0+ mean.true.delta+rnorm(n_obs)

  true.delta<- (Z1*0-Z2*1.2)*mx[,1]+ mean.true.delta+rnorm(n_obs,0,0.5) # TRue weighted average txt effect is 0.70
  
#  true.delta<- txt+ mean.true.delta+rnorm(n_obs,0,0.5)
#  true.delta<- 1*txt*(Z1-1*Z2)+ mean.true.delta+rnorm(n_obs)
#  true.delta<- 0.5*txt*(1*Z1-1*Z2)+ mean.true.delta+rnorm(n_obs)
  cor( true.delta,mx[,1])
  y<-true.delta*txt+mx[,2]
  
  y<-scale(y,scale=FALSE)

y.star=y-mean(true.delta)*txt
y.star=y-0.75*txt
anova(lm(y.star ~ Z))


dat.temp<-data.frame(cbind(y,txt2,Z))
###############2SLS###########################
temp<-fitted(lm(txt~Z))
herto_effect[i]<-mean(y*(temp-mean(txt)))/mean(temp*(temp-mean(txt)))
#fm2<-ivreg(y~ txt2 | Z,data= dat.temp )
#summary(fm2,vcov=sandwich);

summary(lm(y~txt2,data= dat.temp))
##########################################

gamma3<-gamma2<-t(coef(lm(X ~Z-1)))
#pval.AR<-pval.K<-pval.KJ14 <-pval.KJ32 <-pval.J <-NULL
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
#if(LM2< qchisq(1-0.05,1,ncp= 0,lower.tail=TRUE)){print(c(delta,"KL CI"))}
pval.K[i,j]<-1-pf(LM2,1,n_obs-K-1,ncp= 0,lower.tail=TRUE)

###################Kleibergen 2006####################
pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-u%*%sev/see)
D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(see + svv[1,1]* sdelt^2 -2*(sdelt)*sev[1])

J_KLM<-t(u)%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%(u)
ncpJ<-0

pval.J[i,j]<-1-pchisq(J_KLM,K-1,ncp= 0,lower.tail=TRUE)

pval.KJ32[i,j]<-1-pchisq(J_KLM,K-1,ncp= 0,lower.tail=TRUE)* pchisq(LM2,1,ncp= 0,lower.tail=TRUE)

pval.KJ14[i,j]<-1-pchisq(J_KLM,K-1,ncp= 0,lower.tail=TRUE)* pchisq(LM2,1,ncp= 0,lower.tail=TRUE)

##################AR#######################
AR<-t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u/(see + svv[1,1]* sdelt^2 -2*(sdelt)*sev[1])/K
ncpAR<-0
AR

#if(AR< qf(1-0.05,K,n_obs-K,ncp=0,lower.tail=TRUE)){print(c(delta,"AR CI"))}

pval.AR[i,j]<-1-pf(AR,K,n_obs-K-1,ncp= 0,lower.tail=TRUE)

#print(c(AR,LM2))
#if(j%%10==0) print(j)
}
CI.AR[i,]<-(pval.AR[i,]>0.05)*1
CI.K[i,]<-(pval.K[i,]>0.05)*1
CI.KJ14[i,]<-(pval.K[i,]>0.04)*(pval.J[i,]>0.01)
CI.KJ41[i,]<-(pval.K[i,]>0.01)*(pval.J[i,]>0.04)
CI.KJ32[i,]<-(pval.K[i,]>0.02)*(pval.J[i,]>0.03)
if(sum(pval.AR[i,]>0.05)>0){cont<-cont+1}
if(sum(pval.K[i,]>0.05)>0){contK<-contK+1}
if(sum(pval.K[i,]>0.04)>0 & sum(pval.J[i,]>0.01)>0){contKJ14<-contKJ14+1}
if(i%%10==0) print(i)

}
cont/irep
contK/irep


cbind(pval.AR=apply(pval.AR,2,mean),pval.K=apply(pval.K,2,mean),pval.J=apply(pval.J,2,mean),beta=delt.vec,CI.AR=apply(CI.AR,2,mean),CI.K =apply(CI.K,2,mean),CI.KJ14 =apply(CI.KJ14,2,mean),CI.KJ41 =apply(CI.KJ41,2,mean),CI.KJ32 =apply(CI.KJ32,2,mean) )

tab<-data.frame(cbind(pval.AR=apply(pval.AR,2,mean),pval.K=apply(pval.K,2,mean),pval.J=apply(pval.J,2,mean),beta=delt.vec,CI.AR=apply(CI.AR,2,mean),CI.K =apply(CI.K,2,mean),CI.KJ14 =apply(CI.KJ14,2,mean) ))
round(tab,3)
mean(herto_effect)
sd(herto_effect)

xtable(cbind(tab[,4],tab[,-4]))

library(xtable)
xtable(round(cbind(tab[,4],tab[,-4],tab2[,-4]),3))


