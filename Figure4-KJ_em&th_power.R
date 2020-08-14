

# This code generates Fig 5 in the manuscript. It compares the theoretical amd the empirical power of the KJ statistic with a continoues outcome and one treatment options.
# Line 33 specifie the grid on beta1 and beta2. 
# number of replicates is specified in line 32. Fig 5 uses nrep=500.
# K is the dimension of the vector of instruments
# n_obs is the sample size
# gamma sets the strength of the instruments
# delta is the sensitivity parameter
# nrep is the number of simulated datasets


rm(list=ls())
library(ivpack);
library(MASS);library(ivmodel);library(Matrix); #library(expm)


PowerPlotEmTh<-function(K,n_obs, gamma, delta){


#K<-50*1 # number of instruments
#n_obs<-500
strength<-c(1,rep(0.0,K-1))*gamma
#sdelt<- 0.2
sdelt<-delta
expit<-function(x){
	exp(x)/(1+exp(x))
	}

tsls.bet<-bet<-bet.ols<-sd.bet.ols<-sd.bet<-csd.bet<-sd.bet2 <-delt.full <-th.power<-th.power.K <-AR.unad.type1 <-th.power.KJ<-NULL
#set.seed(300)
k06type1.er<- AR.type1.er<-type1.er<-type1.er.naive<-AR.unad.type1.er<-0
nrep<-500
true.delta.vec<-seq(-2,3,by=.02)
#true.delta.vec<-seq(-2,3,by=.2)
true.gam<-"unknown"

th.power.K <-th.power.KJ32 <-th.power.KJ14 <-matrix(0,ncol=nrep,nrow=length(true.delta.vec))
em.power.KJ32<-em.power.KJ14<-NULL
for(i in 1:length(true.delta.vec)){
true.delta<- true.delta.vec[i]

k0632type1.er<-k0614type1.er<-type1.er<-0
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

txt<-scale(txt,scale=FALSE)


X<-cbind(txt) # design matrix of the outcome model
  
  y<-true.delta*txt+mx[,2]

##########################################
#P<-Z%*%solve(t(Z)%*%Z)%*%t(Z)
#M<-diag(dim(P)[1])-P
delta<- 1
###################POWER CALCULATION###########
lambda<-true.delta-delta
gamma2<-coef(lm(txt~Z-1))#-2*coef(summary(lm(txt~Z-1)))[,2]
gamma<-strength
if(true.gam=="known"){gamma2<-gamma}

u<-y-X%*%delta-sdelt* Z%*% gamma2

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

ncpalt.distK<-((lambda* gamma-1*sdelt* gamma)%*%t(Z)%*% P_Ytilda_fun(gamma2)%*%Z%*%(lambda* gamma-1*sdelt* gamma))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)

#ncpalt.distK<-((lambda* gamma-1*sdelt* gamma)%*%t(Z)%*% P_Ytilda_fun(gamma)%*%Z%*%(lambda* gamma-1*sdelt* gamma))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)



LM2<-t(u)%*%(P_Ytilda_fun(gamma2))%*%u/(see + svv* sdelt^2 -2*(sdelt)*sev)


if(LM2> qchisq(1-0.05,1,ncp=0,lower.tail=TRUE)){type1.er<-type1.er+1}
FF<-qf(1-0.05,1,n_obs-K-1,ncp= 0,lower.tail=TRUE)
th.power.K[i,irep]<-1-pf(FF,1,n_obs-K-1,ncp= ncpalt.distK,lower.tail=TRUE)


###################Kleibergen 2006####################
pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-(y-X%*%delta-sdelt* Z%*% gamma2)*sev/see)
D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(see + svv* sdelt^2 -2*(sdelt)*sev)
J_KLM<-t(u)%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%(u)
ncpJ<-0

ncpalt.distKJ<-(((lambda* gamma2-1*sdelt* gamma2)%*%t(Z))%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%t((lambda* gamma2-1*sdelt* gamma2)%*%t(Z)))/(1 + 1 *(lambda-sdelt)^2+2*(lambda-sdelt)* 0.5)*(see + svv* sdelt^2 -2*(sdelt)*sev)


FF1<-qchisq(1-0.01,K-1,ncp= ncpJ,lower.tail=TRUE)
FF2<-qchisq(1-0.04,1,ncp= 0,lower.tail=TRUE)
th.power.KJ14[i,irep]<-1-pchisq(FF1,K-1,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,1,ncp= ncpalt.distK,lower.tail=TRUE)

FF1<-qchisq(1-0.03,K-1,ncp= ncpJ,lower.tail=TRUE)
FF2<-qchisq(1-0.02,1,ncp= 0,lower.tail=TRUE)
th.power.KJ32[i,irep]<-1-pchisq(FF1,K-1,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,1,ncp= ncpalt.distK,lower.tail=TRUE)


if(J_KLM<qchisq(1-0.01,K-1,ncp= 0,lower.tail=TRUE)){ if( LM2>qchisq(1-0.04,1,ncp=0,lower.tail=TRUE)){k0614type1.er<-k0614type1.er+1} }else{k0614type1.er<-k0614type1.er+1}


if(J_KLM<qchisq(1-0.03,K-1,ncp= 0,lower.tail=TRUE)){ if( LM2>qchisq(1-0.02,1,ncp=0,lower.tail=TRUE)){k0632type1.er<-k0632type1.er+1} }else{k0632type1.er<-k0632type1.er+1}



if(irep%%10==0){print(irep)}

}

em.power.KJ14[i]<-k0614type1.er/nrep
em.power.KJ32[i]<-k0632type1.er/nrep

if(i%%1==0){print(i)}

}

type1.er/nrep;apply(th.power.K,1,mean)

(em.power.KJ14);apply(th.power.KJ14,1,mean)
(em.power.KJ32);apply(th.power.KJ32,1,mean)


return(list(beta= true.delta.vec,KJ14=th.power.KJ14,KJ32=th.power.KJ32,KJ05=th.power.K,EMKJ32= em.power.KJ32,EMKJ14= em.power.KJ14))

}


ptm<-proc.time()
PP<-PowerPlotEmTh(K=50,n_obs=500,gamma= 0.5, delta =0.2)
proc.time()-ptm


plot(PP$beta,apply(PP$KJ32,1,mean) ,type="n",ylab="Power", xlab=expression(beta),ylim=c(0,1))
smoothingSpline = smooth.spline(PP$beta,apply(PP$KJ32,1,mean), spar=0.3)
lines(smoothingSpline, lty=1)
smoothingSpline = smooth.spline(PP$beta, PP$EMKJ32, spar=0.3)
lines(smoothingSpline, lty=2)
legend("bottomright",c("th-KJ32","em-KJ32"),lty=1:2,cex=0.6,bty="n")



#save(list=ls(),file="Sens-MWIVs-em&th-KJ3214-Sp5-deltp2-n500-K50-1txt.RData")

