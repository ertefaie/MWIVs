


rm(list=ls())
library(ivpack);library(MASS);library(ivmodel);library(Matrix);library(truncnorm)


PowerPlot2D<-function(K,n_obs, gamma, delta){

K<-20*1 # number of instruments
n_obs<-500
gamma<-1.0
strength<-c(1,rep(0.0,K-1))*gamma
strength2<-c(0,1,rep(0.0,K-2))*gamma
sdelt<- delta<-0.0
#sdelt<-0.1
library(ggplot2);library(gridExtra)
library(lme4); library(optmatch)
expit<-function(x){
	exp(x)/(1+exp(x))
	}

bet<-bet.ols<-sd.bet.ols<-sd.bet<-csd.bet<-sd.bet2 <-delt.full <-NULL
#set.seed(301)
k06type1.er<- AR.type1.er<-type1.er<-type1.er.naive<-AR.unad.type1.er<-0
nrep<-20
#temp1<- seq(0,3,by=.1)
#temp2<- seq(0,3,by=.1)

# these are for sdelt=0.2 
temp1<- seq(-0.2,1.5,by=0.05)
temp2<- seq(-0.6,0.8,by=0.05)#seq(0.8,1.2,by=0.005)

# these are for sdelt=0.0 
temp1<- seq(-1,1,by=0.05)
temp2<- seq(-1,1,by=0.05)#seq(0.8,1.2,by=0.005)


true.delta.vec<-cbind(rep(temp1,each=length(temp2)),rep(temp2,times=length(temp1)))

#true.delta.vec<-cbind(1.4,1.4)



true.gam<-"unknown"

th.power<-th.power.K <-th.power.KJ14 <-th.power.KJ32 <-matrix(0,ncol=nrep,nrow=dim(true.delta.vec)[1])
em.power.KJ32<-em.power.KJ14<-em.power.K<-em.power.AR <-NULL



for(i in 1:dim(true.delta.vec)[1]){
k0632type1.er<-k0614type1.er<-k06type1.er<-type1.er<-AR.type1.er<-0
true.delta<- true.delta.vec[i,]

#source('~/true_parm_zeroinf.R', chdir = TRUE)
###################################################

true.delta<- true.delta.vec[i,]
n_obs1<-400000
u<-rnorm(n_obs1) # unmeasured confounder
x1<-rnorm(n_obs1)
x2<-rnorm(n_obs1)
z1<-rnorm(n_obs1)

Sigma2<-matrix(0.0,ncol=K,nrow=K)
diag(Sigma2)<-1
Z<-mvrnorm(n_obs1,c(rep(0,K)),Sigma2) # instruments
dim(Z);head(Z)
#  logitprob<-0+ 1*x1 + 1*x2+Z%*%strength +u

Sigma2<-matrix(0.0,ncol=3,nrow=3)
diag(Sigma2)<-c(1,1,1)
Sigma2[1,2]<-Sigma2[2,1]<-0.5
mx<-mvrnorm(n_obs1,c(0,0,0),Sigma2)
tsig2y<-Sigma2[3,3]
tsig2d1<-Sigma2[1,1]
tsig2d2<-Sigma2[2,2]

  D1<-Z%*%strength+mx[,1]
  D2<-Z%*%strength2+mx[,2]
  X<-cbind(D1,D2) # design matrix of the outcome model
  
  y<-25+true.delta[1]*D1+true.delta[2]*D2+mx[,3]
  q<-quantile(y,prob=c(0,0.75))
  y[y> q[2]]<- q[2]
#  hist(y)
      y<-y.trans<-log((y+0.5)/(q[2]+0.5-y))
#  hist(y)
  y<-scale(y,scale=FALSE)
  true_parm_zeroinf<-coef(lm(y~X))[-1]





#####################################################
for(irep in 1:nrep){
true.delta<- true.delta.vec[i,]
#n_obs<-500000
u<-rnorm(n_obs) # unmeasured confounder
x1<-rnorm(n_obs)
x2<-rnorm(n_obs)
z1<-rnorm(n_obs)

Sigma<-matrix(0.0,ncol=K,nrow=K)
diag(Sigma)<-1
Z<-mvrnorm(n_obs,c(rep(0,K)),Sigma) # instruments
dim(Z);head(Z)
#  logitprob<-0+ 1*x1 + 1*x2+Z%*%strength +u

Sigma<-matrix(0.5,ncol=3,nrow=3)
diag(Sigma)<-c(1,1,1)
Sigma[1,2]<-Sigma[2,1]<-0.5
mx<-mvrnorm(n_obs,c(0,0,0),Sigma)
tsig2y<-Sigma[3,3]
tsig2d1<-Sigma[1,1]
tsig2d2<-Sigma[2,2]

  D1<-Z%*%strength+mx[,1]
  D2<-Z%*%strength2+mx[,2]
  X<-cbind(D1,D2) # design matrix of the outcome model
  
  y<-25+true.delta[1]*D1+true.delta[2]*D2+mx[,3]
  q<-quantile(y,prob=c(0,0.75))
  y[y> q[2]]<- q[2]
#  hist(y)
      y<-y.trans<-log((y+0.5)/(q[2]+0.5-y))
#  hist(y)
  y<-scale(y,scale=FALSE)
#  lm(y~X)
  ut<-y-(true_parm_zeroinf[1]*D1+ true_parm_zeroinf[2]*D2)
  sevt<-as.vector((t(ut)%*%X)/(n_obs-K))-as.vector((t(ut)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
svvt<-((t(X)%*%X)/(n_obs-K))-((t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
seet<-as.vector((t(ut)%*%ut)/(n_obs-K))

#  S<-(resid(lm(X~Z[,1:2])))
#  sevt <-c(cov(S[,1],resid(lm(y~X))),cov(S[,2],resid(lm(y~X))) )
#  seet <-var(resid(lm(y~X)))

##########################################
true.delta <-true_parm_zeroinf
delta<- c(0,0)
###################POWER CALCULATION###########
lambda<- true.delta-delta
gamma3<-gamma2<-t(coef(lm(X ~Z-1)))#-2*coef(summary(lm(txt~Z-1)))[,2]
gamma2[2,1:K]<-0
gamma<-rbind(strength,strength2)
if(true.gam=="known"){gamma2<-gamma}

u<-y-X%*%delta-Z%*% (gamma2[1,])*sdelt
u2<-y-X%*%delta-Z%*% (gamma[1,])*sdelt
###################POWER CALCULATION###########
##################Kleibergen 2002#######################
sev<-as.vector((t(u)%*%X)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
svv<-((t(X)%*%X)/(n_obs-K))-((t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)/(n_obs-K))
see<-as.vector((t(u)%*%u)/(n_obs-K))-as.vector((t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u)/(n_obs-K))



pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-u%*%sev/see)
Ytilde<- Z %*% pi.tilde
P_Ytilda<-Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)


LM2<-t(u)%*%Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)%*%u/(see + svv[1,1]* sdelt^2+ svv[2,2]* (0.0*sdelt)^2 -2*(sdelt)*sev[1]-2*(0.0*sdelt)*sev[2]+2*sdelt*sdelt*0.0*svv[1,2])

LM.naive<-t(u)%*%Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)%*%u/see

if(LM2> qchisq(1-0.05,2,ncp=0,lower.tail=TRUE)){type1.er<-type1.er+1}
if(LM.naive> qchisq(1-0.05,2,ncp=0,lower.tail=TRUE)){type1.er.naive<-type1.er.naive+1}


#ncpalt.distK<-(as.vector(lambda%*% gamma3-(gamma2[1,])*sdelt)%*%t(Z)%*% Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)%*%Z%*% as.vector(lambda%*% gamma3-(gamma2[1,])*sdelt))/(see + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*sev[1]+2*(lambda[2]-0.0*sdelt) *sev[2]+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))

ncpalt.distK<-(as.vector(lambda%*% gamma-(gamma2[1,])*sdelt)%*%t(Z)%*% Ytilde%*%solve(t(Ytilde)%*% Ytilde)%*%t(Ytilde)%*%Z%*% as.vector(lambda%*% gamma-(gamma2[1,])*sdelt))/(seet + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*sevt[1]+2*(lambda[2]-0.0*sdelt) *sevt[2]+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))

#FF<-qf(1-0.05,2,n_obs-K-2,ncp= 0,lower.tail=TRUE)
FF<-qchisq(1-0.05,2,ncp= 0,lower.tail=TRUE)
#th.power.K[irep]<-1-pf(FF,2,n_obs-K-2,ncp= ncpalt.distK,lower.tail=TRUE)
th.power.K[i,irep]<-1-pchisq(FF,2,ncp= ncpalt.distK,lower.tail=TRUE)

###################Kleibergen 2006####################
#see2<-as.vector(t(u)%*%u/(n_obs-K))
pi.tilde<-solve(t(Z)%*%Z)%*%t(Z)%*%(X-u%*%sev/see)
D.comp<-solve(t(Z)%*%Z)%*%(qr.Q(qr(pi.tilde),complete=TRUE)[,-1])
middle.part<-t(D.comp)%*%t(Z)%*%Z%*% D.comp*(see + svv[1,1]* sdelt^2+ svv[2,2]* (0.0*sdelt)^2 -2*(sdelt)*sev[1]-2*(0.0*sdelt)*sev[2]+2*sdelt*sdelt*0.0*svv[1,2])

J_KLM<-t(u)%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%(u)
ncpJ<-0

#ncpalt.distKJ<-(((lambda%*% gamma3-c(sdelt, 0.0)%*%gamma3)%*%t(Z))%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%t((lambda%*% gamma3-c(sdelt, 0.0)%*%gamma3)%*%t(Z)))/(1 + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*.38+2*(lambda[2]-0.0*sdelt) *.38+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))*(see + svv[1,1]* sdelt^2+ svv[2,2]* (0.0*sdelt)^2 -2*(sdelt)*sev[1]-2*(0.0*sdelt)*sev[2]+2*sdelt*sdelt*0.0*svv[1,2])

ncpalt.distKJ<-(((lambda%*% gamma-c(sdelt, 0.0)%*%gamma3)%*%t(Z))%*%Z%*% D.comp%*%solve(middle.part)%*%t(D.comp)%*%t(Z)%*%t((lambda%*% gamma-c(sdelt, 0.0)%*%gamma3)%*%t(Z)))/(seet + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*sevt[1]+2*(lambda[2]-0.0*sdelt) *sevt[2]+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))*(see + svv[1,1]* sdelt^2+ svv[2,2]* (0.0*sdelt)^2 -2*(sdelt)*sev[1]-2*(0.0*sdelt)*sev[2]+2*sdelt*sdelt*0.0*svv[1,2])


FF1<-qchisq(1-0.03,K-2,ncp= 0,lower.tail=TRUE)
FF2<-qchisq(1-0.02,2,ncp= 0,lower.tail=TRUE)
th.power.KJ32[i,irep]<-1-pchisq(FF1,K-2,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,2,ncp= ncpalt.distK,lower.tail=TRUE)

if(J_KLM<qchisq(1-0.03,K-2,ncp= 0,lower.tail=TRUE)){ if( LM2>qchisq(1-0.02,2,ncp=0,lower.tail=TRUE)){k0632type1.er<-k0632type1.er+1} }else{k0632type1.er<-k0632type1.er+1}



FF1<-qchisq(1-0.01,K-2,ncp= 0,lower.tail=TRUE)
FF2<-qchisq(1-0.04,2,ncp= 0,lower.tail=TRUE)
th.power.KJ14[i,irep]<-1-pchisq(FF1,K-2,ncp= ncpalt.distKJ,lower.tail=TRUE)*pchisq(FF2,2,ncp= ncpalt.distK,lower.tail=TRUE)




if(J_KLM<qchisq(1-0.01,K-2,ncp= 0,lower.tail=TRUE)){ if( LM2>qchisq(1-0.04,2,ncp=0,lower.tail=TRUE)){k0614type1.er<-k0614type1.er+1} }else{k0614type1.er<-k0614type1.er+1}

##################AR#######################
AR<-t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u/(see + svv[1,1]* sdelt^2+ svv[2,2]* (0.0*sdelt)^2 -2*(sdelt)*sev[1]-2*(0.0*sdelt)*sev[2]+2*sdelt*sdelt*0.0*svv[1,2])/K
AR.unad<-t(u)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%u/(see)/K

if(AR> qf(1-0.05,K,n_obs-K-1,ncp=0,lower.tail=TRUE)){AR.type1.er<-AR.type1.er+1}
if(AR.unad> qf(1-0.05,K,n_obs-K-1,ncp=0,lower.tail=TRUE)){AR.unad.type1.er<-AR.unad.type1.er+1}

#ncpalt.dist<- (t(as.vector(lambda%*%gamma3-c(sdelt, 0.0)%*%gamma3))%*%t(Z)%*%Z%*%(as.vector(lambda%*%gamma3-c(sdelt, 0.0)%*%gamma3)))/(1 + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*.38+2*(lambda[2]-0.0*sdelt) *.38+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))

#ncpalt.dist<- (t(as.vector(lambda%*%gamma3-c(sdelt, 0.0)%*%gamma3))%*%t(Z)%*%Z%*%(as.vector(lambda%*%gamma3-c(sdelt, 0.0)%*%gamma3)))/(see + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*sev[2]+2*(lambda[2]-0.0*sdelt) *sev[2]+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))

ncpalt.dist<- (t(as.vector(lambda%*%gamma-c(sdelt, 0.0)%*%gamma2))%*%t(Z)%*%Z%*%(as.vector(lambda%*%gamma-c(sdelt, 0.0)%*%gamma2)))/(seet + 1 *(lambda[1]-sdelt)^2+1 *(lambda[2]-0.0*sdelt)^2+2*(lambda[1]-sdelt)*sevt[2]+2*(lambda[2]-0.0*sdelt) *sevt[2]+2*0.5*(lambda[2]-0.0*sdelt)* (lambda[1]-sdelt))

FF<-qf(1-0.05,K,n_obs-K-1,ncp= 0,lower.tail=TRUE)
th.power[i,irep]<-1-pf(FF,K,n_obs-K-1,ncp= ncpalt.dist,lower.tail=TRUE)


if(irep%%100==0){print(irep)}
  }
  em.power.K[i]<-type1.er/nrep
  em.power.KJ14[i]<-k0614type1.er/nrep
  em.power.KJ32[i]<-k0632type1.er/nrep
  em.power.AR[i]<-AR.type1.er/nrep

if(i%%100==0){print(i)}
}
th.mK32<-matrix(apply(th.power.KJ32,1,mean), nrow=length(temp1),ncol=length(temp2),byrow=TRUE)
th.mK14<-matrix(apply(th.power.KJ14,1,mean), nrow=length(temp1),ncol=length(temp2),byrow=TRUE)
th.m<-matrix((apply(th.power,1,mean)), nrow=length(temp1),ncol=length(temp2),byrow=TRUE)

return(list(beta= true.delta.vec,KJ14= th.mK14,KJ32= th.mK32, AR= th.m,beta1=temp1, beta2=temp2))

}

PP<-PowerPlot2D(K=20,n_obs=500,gamma= 1, delta =0.0)

plot(PP$beta1, PP$beta2,type="n",ylim=c(-.3,.3),xlim=c(-.25,.25),main=expression(paste("K=20,", delta, "=0.0,", gamma, "=1.0")), xlab= expression(beta[1]), ylab= expression(beta[2]))
contour(PP$beta1, PP$beta2, PP$KJ14,add=TRUE)
contour(PP$beta1, PP$beta2, PP$KJ32,add=TRUE,lty=2)
contour(PP$beta1, PP$beta2, PP$AR,add=TRUE,lty=3)



