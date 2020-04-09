
  

# This code generates Fig 6 in the manuscript. It compares the power of AR and different KJ statistics when the outcome is continoues and there is one treatment option 
# K is the dimension of the vector of instruments
# n_obs is the sample size
# gamma sets the strength of the instruments
# delta is the sensitivity parameter
# nrep is the number of simulated datasets



rm(list=ls())
library(ivpack);library(MASS);library(ivmodel);library(Matrix)

PowerPlot<-function(K,n_obs, gamma, delta){

#K<-20*1 # number of instruments
#n_obs<-500
strength<-c(1,rep(0.0,K-1))*gamma
#sdelt<- 0.2
sdelt<-delta
library(ggplot2);library(gridExtra)
library(lme4); library(optmatch)
expit<-function(x){
	exp(x)/(1+exp(x))
	}

bet<-bet.ols<-sd.bet.ols<-sd.bet<-csd.bet<-sd.bet2 <-delt.full <-NULL
#set.seed(301)
k06type1.er<- AR.type1.er<-type1.er<-type1.er.naive<-AR.unad.type1.er<-0
nrep<-10
true.delta.vec<- c(seq(-2,.95,by=.05),seq(.95,1.05,by=.01),seq(1.05,4,by=.05))
true.delta.vec<- seq(-2,4,by=.015)
true.gam<-"unknown"

th.power<-th.power.K <-th.power.KJ14 <-th.power.KJ32 <-matrix(0,ncol=nrep,nrow=length(true.delta.vec))



for(i in 1:length(true.delta.vec)){
true.delta<- true.delta.vec[i]

for(irep in 1:nrep){


u<-rnorm(n_obs) # unmeasured confounder
x1<-rnorm(n_obs)
x2<-rnorm(n_obs)
z1<-rnorm(n_obs)

Sigma<-matrix(0.0,ncol=K,nrow=K)
diag(Sigma)<-1
Z<-mvrnorm(n_obs,c(rep(0,K)),Sigma) # instruments
dim(Z);head(Z)
#  logitprob<-0+ 1*x1 + 1*x2+Z%*%strength +u

Sigma<-matrix(0.5,ncol=2,nrow=2)
diag(Sigma)<-c(1,1)
mx<-mvrnorm(n_obs,c(0,0),Sigma)
tsig2y<-Sigma[2,2]
tsig2d<-Sigma[1,1]

  logitprob<-Z%*%strength+mx[,1] 
  txt<-logitprob

#  strength2<-runif(K,0,0.5)
#  logitprob<-Z%*%strength2+mx[,1] 
#  txt<-logitprob



txt<-scale(txt,scale=FALSE)


X<-cbind(txt) # design matrix of the outcome model
  
  y<-true.delta*txt+mx[,2]
  

##########################################
delta<- 1
###################POWER CALCULATION###########
lambda<-true.delta-delta
gamma2<-coef(lm(txt~Z-1))#-2*coef(summary(lm(txt~Z-1)))[,2]
gamma<-strength
if(true.gam=="known"){gamma2<-gamma}

u<-y-X%*%delta-sdelt* Z%*% gamma2
###################POWER CALCULATION###########
##################Kleibergen 2002#######################
sev<-as.vector((t(u)%*%X)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
svv<-as.vector((t(X)%*%X)/(n_obs-K))-as.vector((t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
see<-as.vector((t(u)%*%u)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u)/(n_obs-K))

P_Ytilda_fun<-function(gam){
pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-(y-X%*%delta-sdelt* Z%*% gam)*sev/see)
Ytilde<- Z %*% pi.tilde
P_Ytilda<-Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)
P_Ytilda
}
ncp<-0
ncpalt.distK<-((lambda*gamma-1*sdelt* gamma)%*%t(Z)%*% P_Ytilda_fun(gamma)%*%Z%*%(lambda*gamma-1*sdelt* gamma))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)




LM2<-t(u)%*%(P_Ytilda_fun(gamma2))%*%u/(see + svv* sdelt^2 -2*(sdelt)*sev)
LM.naive<-t(u)%*%(P_Ytilda_fun(gamma2))%*%u/see#*K/n_obs
if(LM2> qchisq(1-0.05,1,ncp= ncp,lower.tail=TRUE)){type1.er<-type1.er+1}
if(LM.naive> qchisq(1-0.05,1,ncp= ncp,lower.tail=TRUE)){type1.er.naive<-type1.er.naive+1}
FF<-qf(1-0.05,1,n_obs-K-1,ncp= ncp,lower.tail=TRUE)
th.power.K[i,irep]<-1-pf(FF,1,n_obs-K-1,ncp= ncpalt.distK,lower.tail=TRUE)


###################Kleibergen 2006####################
see2<-as.vector(t(u)%*%u/(n_obs-K))
pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-(y-X%*%delta-sdelt* Z%*% gamma2)*sev/see)
D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(see + svv* sdelt^2 -2*(sdelt)*sev)
J_KLM<-t(u)%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%(u)

ncpJ<-0

#pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-(y-X%*%delta-sdelt* Z%*% gamma)*sev/see)
#D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
#middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(1 + 1* sdelt^2 -2*(sdelt)*0.5)
ncpalt.distKJ<-(((lambda*gamma-1*sdelt* gamma)%*%t(Z))%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%t((lambda*gamma-1*sdelt* gamma)%*%t(Z)))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)*(see + svv* sdelt^2 -2*(sdelt)*sev)


FF1<-qchisq(1-0.03,K-1,ncp= ncpJ,lower.tail=TRUE)
FF2<-qchisq(1-0.02,1,ncp= ncp,lower.tail=TRUE)
th.power.KJ32[i,irep]<-1-pchisq(FF1,K-1,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,1,ncp= ncpalt.distK,lower.tail=TRUE)

FF1<-qchisq(1-0.01,K-1,ncp= ncpJ,lower.tail=TRUE)
FF2<-qchisq(1-0.04,1,ncp= ncp,lower.tail=TRUE)
th.power.KJ14[i,irep]<-1-pchisq(FF1,K-1,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,1,ncp= ncpalt.distK,lower.tail=TRUE)




if(J_KLM<qchisq(1-0.03,K-1,ncp= ncpJ,lower.tail=TRUE)){ if( LM2>qchisq(1-0.02,1,ncp=0,lower.tail=TRUE)){k06type1.er<-k06type1.er+1} }else{k06type1.er<-k06type1.er+1}

##################AR#######################
AR<-t(u)%*%(Z%*%solve(t(Z)%*%Z)%*%t(Z)/K)%*%u/(see + svv* sdelt^2 -2*(sdelt)*sev)
if(AR> qf(1-0.05,K,n_obs-K-1,ncp=0,lower.tail=TRUE)){AR.type1.er<-AR.type1.er+1}

ncp<-0
ncpalt.dist<-(t(lambda*gamma-1*sdelt* gamma)%*%t(Z)%*%Z%*%(lambda*gamma-1*sdelt* gamma))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)


FF<-qf(1-0.05,K,n_obs-K-1,ncp= ncp,lower.tail=TRUE)
th.power[i,irep]<-1-pf(FF,K,n_obs-K-1,ncp= ncpalt.dist,lower.tail=TRUE)


  }
if(i%%1==0){print(i)}
}


return(list(beta= true.delta.vec,KJ14=th.power.KJ14,KJ32=th.power.KJ32,KJ05=th.power.K, AR=th.power))

}

ptm<-proc.time()
PP<-PowerPlot(K=50,n_obs=500,gamma= 1, delta =0.2)
proc.time()-ptm


plot(PP$beta, apply(PP$KJ14,1,mean),type="n",ylab="Power", xlab=expression(beta),ylim=c(0,1),xlim=c(-2,2.5))
smoothingSpline = smooth.spline(PP$beta, apply(PP$KJ14,1,mean), spar=0.3)
lines(smoothingSpline,col=4,lty=2,lwd=2)
smoothingSpline = smooth.spline(PP$beta, apply(PP$KJ32,1,mean), spar=0.25)
lines(smoothingSpline,col=3)
smoothingSpline = smooth.spline(PP$beta, apply(PP$KJ05,1,mean), spar=0.2)
lines(smoothingSpline,col=2)
smoothingSpline = smooth.spline(PP$beta, apply(PP$AR,1,mean), spar=0.4)
lines(smoothingSpline,col=1)
legend("topright",c("AR","KJ05","KJ32","KJ14"),col=c(1,2,3,4),lty=c(1,1,1,2),cex=0.8,bty="n", text.font=1)


