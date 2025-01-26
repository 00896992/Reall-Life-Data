#R-Code for Robust Ridge Regression in presence of outliers (For real data)
#################Muhammad Suhail########################
#############Roll no P2-14########################
rm(list=ls())
set.seed(1990)
library(MASS)
p=4                              #No of explanatory variables
I=diag(p)
setwd("D:\\mphilwork\\Real life Data")      #To change working directory
data=read.csv("Tobacco.csv")
#data=data[1:30,]
#data=data[1:30,]
x1=data$x1          # use $ to add variables
x2=data$x2
x3=data$x3
x4=data$x4
#x5=data$x5            # use $ to add variables
#x6=data$x6
#x7=data$x7
#x8=data$x8
#x9=data$x9

n=length(x1)
x=cbind(x1,x2,x3,x4)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
y=data$y
y=(y-mean(y))/sd(y)

#y=scale(x,center = TRUE,scale = TRUE)
eigen(cor(x))
ev=eigen(cor(x))$values
vec=eigen(cor(x))$vectors
which.max(ev)
CN=max(ev)/min(ev)        #Computing condition number
CI=sqrt(max(ev)/min(ev))
beta=vec[,which.max(ev)]
C=t(x)%*%(x)
D=eigen(C)$vectors
Z=x%*%D
lam=t(Z)%*%Z
lam=round(lam,4)
lamda=diag(lam)
alpha=t(D)%*%beta

betahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=x%*%betahat
sigmahat=(sum((y-yhat)^2))/(n-p)
alphahat=solve(lam)%*%t(Z)%*%y
#y[5]=y[5]+10*sigmahat
#y[15]=y[15]-5*sigmahat
#y[20]=y[20]+20*sigmahat
#y[35]=y[35]-750*sigmahat
#y[40]=y[40]+1000*sigmahat

#Estimation of ridge parameter
#OLS estimator
alphahat.ols=alphahat

#Ridge Regression estimator of HK method (1970)

HK=sigmahat/(max(alphahat.ols^2))
K1=HK
alphahatk1=solve(lam+I*K1)%*%lam%*%alphahat.ols
alphahatk.HK = alphahatk1

# HKB (1976)
HKB=(p*sigmahat)/(sum(alphahat.ols^2))
K2=HKB
alphahatk2=solve(lam+I*K2)%*%lam%*%alphahat.ols
alphahatk.HKB = alphahatk2

# Kibria (2003) A.M,GM & Median

KAM=(sum(sigmahat/(alphahat.ols^2)))/p
K3=KAM
alphahatk3=solve(lam+I*K3)%*%lam%*%alphahat.ols
alphahatk.KAM = alphahatk3


KGM=sigmahat/(prod(alphahat.ols^2))^(1/p)
K4=KGM
alphahatk4=solve(lam+I*K4)%*%lam%*%alphahat.ols
alphahatk.KGM = alphahatk4


KMED=median(sigmahat/(alphahat.ols^2))
K5=KMED
alphahatk5=solve(lam+I*K5)%*%lam%*%alphahat.ols
alphahatk.KMED= alphahatk5


# Khalaf,Mansson & Shukur (2013)

kms=max(ev)/(sum((alphahat.ols^2)))
KMS=kms*HK
K6=KMS
alphahatk6=solve(lam+I*K6)%*%lam%*%alphahat.ols
alphahatk.KMS = alphahatk6

# LC estimator
LC=HK
K7=LC
q.num7=(t(y)%*%Z)%*%solve(lam+I*K7)%*%lam%*%alphahat.ols
q.den7=(t(y)%*%Z)%*%solve(lam+I*K7)%*%lam%*%solve(lam+I*K7)%*%lam%*%alphahat.ols
q7=q.num7/q.den7
q7=c(q7)
alphahatk7=q7*(solve(lam+I*K7)%*%lam%*%alphahat.ols)
alphahatk.LC = alphahatk7

# Selma Thoker 
ST=HK
K8=ST
q.opt.N8=sum(((alphahat.ols^2)*ev)/(ev+K8))
q.opt.D8=sum(((sigmahat*ev)+((alphahat.ols^2)*(ev^2)))/((ev+K8)^2))
q8=q.opt.N8/q.opt.D8
q8=c(q8)
k.opt.N8=((q8*(sum(sigmahat*ev)))+((q8-1)*sum((alphahat.ols^2)*(ev^2))))
k.opt.D8=sum((alphahat.ols^2)*ev)
k.opt8=k.opt.N8/k.opt.D8
alphahatk8=q8*(solve(lam+I*k.opt8)%*%lam%*%alphahat.ols)
alphahatk.ST = alphahatk8

#Shakir (2023)
KL1=sigmahat/(2*alphahat+(sigmahat/ev))
W=max(ev)/sum(abs(alphahat))
KL=W*KL1

NQD1=quantile(KL,0.05)
NQD2=quantile(KL,0.25)
NQD3=quantile(KL,0.50)
NQD4=quantile(KL,0.75)
NQD5=quantile(KL,0.95)
NQD6=quantile(KL,0.99)



K9=NQD1
q.num9=(t(y)%*%Z)%*%solve(lam+I*K9)%*%lam%*%alphahat.ols
q.den9=(t(y)%*%Z)%*%solve(lam+I*K9)%*%lam%*%solve(lam+I*K9)%*%lam%*%alphahat.ols
q9=q.num9/q.den9
q9=c(q9)
alphahatk9=q9*(solve(lam+I*K9)%*%lam%*%alphahat.ols)
alphahatk.NQD1= alphahatk9


K10=NQD2
q.num10=(t(y)%*%Z)%*%solve(lam+I*K10)%*%lam%*%alphahat.ols
q.den10=(t(y)%*%Z)%*%solve(lam+I*K10)%*%lam%*%solve(lam+I*K10)%*%lam%*%alphahat.ols
q10=q.num10/q.den10
q10=c(q10)
alphahatk10=q10*(solve(lam+I*K10)%*%lam%*%alphahat.ols)
alphahatk.NQD2= alphahatk10


K11=NQD3
q.num11=(t(y)%*%Z)%*%solve(lam+I*K11)%*%lam%*%alphahat.ols
q.den11=(t(y)%*%Z)%*%solve(lam+I*K11)%*%lam%*%solve(lam+I*K11)%*%lam%*%alphahat.ols
q11=q.num11/q.den11
q11=c(q11)
alphahatk11=q11*(solve(lam+I*K11)%*%lam%*%alphahat.ols)
alphahatk.NQD3= alphahatk11


K12=NQD4
q.num12=(t(y)%*%Z)%*%solve(lam+I*K12)%*%lam%*%alphahat.ols
q.den12=(t(y)%*%Z)%*%solve(lam+I*K12)%*%lam%*%solve(lam+I*K12)%*%lam%*%alphahat.ols
q12=q.num12/q.den12
q12=c(q12)
alphahatk12=q12*(solve(lam+I*K12)%*%lam%*%alphahat.ols)
alphahatk.NQD4= alphahatk12


K13=NQD5
q.num13=(t(y)%*%Z)%*%solve(lam+I*K13)%*%lam%*%alphahat.ols
q.den13=(t(y)%*%Z)%*%solve(lam+I*K13)%*%lam%*%solve(lam+I*K13)%*%lam%*%alphahat.ols
q13=q.num13/q.den13
q13=c(q13)
alphahatk13=q13*(solve(lam+I*K13)%*%lam%*%alphahat.ols)
alphahatk.NQD5= alphahatk13


K14=NQD6
q.num14=(t(y)%*%Z)%*%solve(lam+I*K14)%*%lam%*%alphahat.ols
q.den14=(t(y)%*%Z)%*%solve(lam+I*K14)%*%lam%*%solve(lam+I*K14)%*%lam%*%alphahat.ols
q14=q.num14/q.den14
q14=c(q14)
alphahatk14=q14*(solve(lam+I*K14)%*%lam%*%alphahat.ols)
alphahatk.NQD6= alphahatk14


comb.coeff=rbind(alphahat.ols,alphahatk.HK,alphahatk.HKB,alphahatk.KAM,alphahatk.KGM,alphahatk.KMED,
                 alphahatk.KMS,alphahatk.LC,alphahatk.ST,alphahatk.NQD1,alphahatk.NQD2,alphahatk.NQD3
                 ,alphahatk.NQD4,alphahatk.NQD5,alphahatk.NQD6)
comb.coeff=round(comb.coeff,4)

###########################################################################################################
#MSE of all estimators
#MSE of OLS and RR

KR1<-c(0,K1,K2,K3,K4,K5,K6) 
a1=rep(0,length(KR1)); b1=rep(0,length(KR1))                            #to obtain MSE of ridge estimators
for(i in 1:length(KR1)){
  
  a1[i]=sum((ev)/(ev+KR1[i])^2)
  b1[i]=sum(((KR1[i]^2)*(alphahat.ols^2))/(ev+KR1[i])^2)
  
}
MSE1=sigmahat*a1+b1

KTR1<-c(K7,K8,K9,K10,K11,K12,K13,K14) 
qTR<-c(q7,q8,q9,q10,q11,q12,q13,q14)
a2=rep(1,length(KTR1)); b2=rep(1,length(KTR1)); q1=rep(1,length(qTR))                             #to obtain MSE of ridge estimators
for(i in 1:length(KTR1)){
  
  
  q1[i]=qTR[i]
  a2[i]=((q1[i]^2)*(sigmahat)*sum((ev)/(ev+KTR1[i])^2))
  b2[i]=sum(((((q1[i]*ev)/((ev+KTR1[i]))-1)^2))*(alphahat.ols)^2)
  
} 

MSE2=a2+b2


MSE=c(MSE1,MSE2)
MSE=as.matrix(MSE)
MSE=round(MSE,8)
K.hat=c(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14)
write.csv(K.hat,file = "D:\\mphilwork\\Real life Data\\K.hat.Tobacco.csv")
row.names(MSE)=c("OLS","HK","HKB","KAM","KGM","KMED","KMS","LC","ST","NQD1","NQD2","NQD3","NQD4","NQD5","NQD6")
colnames(MSE)=c("MSE")
write.csv(MSE,file = "D:\\mphilwork\\Real life Data\\MSE.Tobacco.csv")
write.csv(comb.coeff,file = "D:\\mphilwork\\Real life Data\\Coefficient.Tobacco.csv")

cor.x=round(cor(x),4)
cor.y=round(cor(x,y),4)
cor.xy=cbind(cor.x,cor.y)
write.csv(cor.xy,file = "D:\\mphilwork\\Real life Data\\Correlation.Tobacco.csv")
CN
CI
ev
