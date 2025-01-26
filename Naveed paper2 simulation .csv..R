#========R-Code for Quantile based estimation of biasing parameters in  ridge regression model(2020)========
#==========Authors: Danish Wasim PhD scholar===========
rm(list = (ls()))
set.seed(1234)
sigma1 =c(0.1,1,5,10)   #standard deviation
ssize=c(20,60,120)        #sample size
pred=c(4,10)         #Independent variables
corr=c(0.85,0.95,0.99,0.999)  #Multicollinearity
N=10000  #Simulation runs or number of runs

MSE.r=matrix(0,15,length(corr))
MSE.s=array(0,dim = c(15,length(corr),length(sigma1)))
MSE.n=array(0,dim = c(15,length(corr),length(sigma1),length(ssize)))
MSE.p=array(0,dim = c(15,length(corr),length(sigma1),length(ssize),length(pred)))

#-----------------------
#Predictors loop
for(b in 1:length(pred)){
  p=pred[b]
  I=diag(p)
  
  #------------------------
  #Sample size loop
  for(a in 1:length(ssize)){
    n=ssize[a]
    
    z=matrix(0,n,p)
    for(i in 1:p){            #Generating stanndard normal random variable Z
      z[,i]=rnorm(n,0,1)
    }
    x=matrix(0,n,p)
    #--------------------------
    #Error variance loop
    for(s in 1:length(sigma1)){
      sigma=sigma1[s]
      for(q in 1:length(corr)){     #Correlation loop
        r=corr[q]
        
        for(i in 1:p){            #Generating X matrix
          x[,i]=sqrt(1-(r^2))*z[,i]+r*z[,p]
        }
        #x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
        
        C=t(x)%*%(x)
        eigen(C)
        ev=eigen(C)$values
        D=eigen(C)$vectors
        beta=D[,which.max(ev)]
        lamb=diag(t(D)%*%C%*%D)
        Z=x%*%D
        lam=t(Z)%*%Z
        alpha=t(D)%*%beta
        
        betahat_ols=matrix(ncol=N, nrow=p)
        alphahatk_M1 = matrix(ncol=N, nrow=p)
        alphahatk_M2 = matrix(ncol=N, nrow=p)
        alphahatk_M3 = matrix(ncol=N, nrow=p)
        alphahatk_M4 = matrix(ncol=N, nrow=p)
        alphahatk_M5 = matrix(ncol=N, nrow=p)
        alphahatk_M6 = matrix(ncol=N, nrow=p)
        alphahatk_M7 = matrix(ncol=N, nrow=p)
        alphahatk_M8 = matrix(ncol=N, nrow=p)
        alphahatk_M9 = matrix(ncol=N, nrow=p)
        alphahatk_M10 = matrix(ncol=N, nrow=p)
        alphahatk_M11 = matrix(ncol=N, nrow=p)
        alphahatk_M12 = matrix(ncol=N, nrow=p)
        alphahatk_M13 = matrix(ncol=N, nrow=p)
        alphahatk_M14 = matrix(ncol=N, nrow=p)
        
        
        for(i in 1:N){          #Simulation loop starts here
          e=rnorm(n,0,sigma)
          y=Z%*%alpha+e
          betahat=solve(t(x)%*%x)%*%t(x)%*%y
          yhat=x%*%betahat
          sigmahat=(sum((y-yhat)^2))/(n-p)
          sigmastarhat=(sum((y-yhat)^2))/(n-p+2)
          alphahat=solve(diag(lamb))%*%t(Z)%*%y
          
          yhat1=Z%*%alphahat
          sigmahat1=(sum((y-yhat1)^2))/(n-p)
          
          #OLS method
          
          betahat_ols[,i]=c(betahat)
          
          #HK method (1970)
          
          HK=sigmahat/(max(alphahat^2))
          K1=HK
          alphahatk<-solve(t(Z)%*%(Z)+I*K1)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M1[,i] = alphahatk
          
          #HKB method (1975)
          
          HKB=(p*sigmahat)/(sum(alphahat^2))
          K2=HKB
          alphahatk<-solve(t(Z)%*%(Z)+I*K2)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M2[,i] = alphahatk
          
          # Kibria (2003) A.M,GM & Median
          
          KAM=(sum(sigmahat/(alphahat^2)))/p
          K3=KAM
          alphahatk<-solve(t(Z)%*%(Z)+I*K3)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M3[,i] = alphahatk
          
          KGM=sigmahat/(prod(alphahat^2))^(1/p)
          K4=KGM
          alphahatk<-solve(t(Z)%*%(Z)+I*K4)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M4[,i] = alphahatk
          
          KMED=median(sigmahat/(alphahat^2))
          K5=KMED
          alphahatk<-solve(t(Z)%*%(Z)+I*K5)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M5[,i] = alphahatk
          
          
          # Khalaf,Mansson & Shukur (2013)
          
          KMSI=max(ev)/(sum((alphahat^2)))
          KMS1=KMSI*HK
          K6=KMS1
          alphahatk<-solve(t(Z)%*%(Z)+I*K6)%*%t(Z)%*%y
          alphahatk<-c(alphahatk)
          alphahatk_M6[,i] = alphahatk
          
          # LC estimator
          LC=HK
          K7=LC
          q.num7=(t(y)%*%Z)%*%solve(lam+I*K7)%*%(t(Z)%*%y)
          q.den7=(t(y)%*%Z)%*%solve(lam+I*K7)%*%lam%*%solve(lam+I*K7)%*%(t(Z)%*%y)
          q7=q.num7/q.den7
          q7=c(q7)
          alphahatk=q7*(solve(t(Z)%*%(Z)+I*K7)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M7[,i] = alphahatk
          
          # Selma Thoker 
          ST=HK
          K8=ST
          q.opt.N8=sum((alphahat^2)*ev)/(ev+K8)
          q.opt.D8=sum((sigmahat*ev+(alphahat^2)*ev^2)/(ev+K8)^2)
          q.opt8=q.opt.N8/q.opt.D8
          q.opt8=c(q.opt8)
          k.opt.N8=q.opt8*(sum(sigmahat*ev))+(q.opt8-1)*sum((alphahat^2)*ev)
          k.opt.D8=sum((alphahat^2)*ev)
          k.opt8=k.opt.N8/k.opt.D8
          k.opt8=c(k.opt8)
          alphahatk=q.opt8*(solve(t(Z)%*%(Z)+I*k.opt8)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M8[,i] = alphahatk
          
          # Quantile KL estimator (2020)
          KL1=sigmahat/(2*alphahat+(sigmahat/ev))
          W=max(ev)/sum(abs(alphahat))
          KL=W*KL1
          
          QTPKL1=quantile(KL,0.05)
          QTPKL2=quantile(KL,0.25)
          QTPKL3=quantile(KL,0.50)
          QTPKL4=quantile(KL,0.75)
          QTPKL5=quantile(KL,0.95)
          QTPKL6=quantile(KL,0.99)
          #NQW=HK1
          
        
          K9=QTPKL1
          q.num9=(t(y)%*%Z)%*%solve(lam+I*K9)%*%(t(Z)%*%y)
          q.den9=(t(y)%*%Z)%*%solve(lam+I*K9)%*%lam%*%solve(lam+I*K9)%*%(t(Z)%*%y)
          q9=q.num9/q.den9
          q9=c(q9)
          alphahatk=q9*(solve(t(Z)%*%(Z)+I*K9)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M9[,i] = alphahatk
          
          K10=QTPKL2
          q.num10=(t(y)%*%Z)%*%solve(lam+I*K10)%*%(t(Z)%*%y)
          q.den10=(t(y)%*%Z)%*%solve(lam+I*K10)%*%lam%*%solve(lam+I*K10)%*%(t(Z)%*%y)
          q10=q.num10/q.den10
          q10=c(q10)
          alphahatk=q10*(solve(t(Z)%*%(Z)+I*K10)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M10[,i] = alphahatk
          
          K11=QTPKL3
          q.num11=(t(y)%*%Z)%*%solve(lam+I*K11)%*%(t(Z)%*%y)
          q.den11=(t(y)%*%Z)%*%solve(lam+I*K11)%*%lam%*%solve(lam+I*K11)%*%(t(Z)%*%y)
          q11=q.num11/q.den11
          q11=c(q11)
          alphahatk=q11*(solve(t(Z)%*%(Z)+I*K11)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M11[,i] = alphahatk
          
          K12=QTPKL4
          q.num12=(t(y)%*%Z)%*%solve(lam+I*K12)%*%(t(Z)%*%y)
          q.den12=(t(y)%*%Z)%*%solve(lam+I*K12)%*%lam%*%solve(lam+I*K12)%*%(t(Z)%*%y)
          q12=q.num12/q.den12
          q12=c(q12)
          alphahatk=q12*(solve(t(Z)%*%(Z)+I*K12)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M12[,i] = alphahatk
          
          K13=QTPKL5
          q.num13=(t(y)%*%Z)%*%solve(lam+I*K13)%*%(t(Z)%*%y)
          q.den13=(t(y)%*%Z)%*%solve(lam+I*K13)%*%lam%*%solve(lam+I*K13)%*%(t(Z)%*%y)
          q13=q.num13/q.den13
          q13=c(q13)
          alphahatk=q13*(solve(t(Z)%*%(Z)+I*K13)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M13[,i] = alphahatk
          
          K14=QTPKL6
          q.num14=(t(y)%*%Z)%*%solve(lam+I*K14)%*%(t(Z)%*%y)
          q.den14=(t(y)%*%Z)%*%solve(lam+I*K14)%*%lam%*%solve(lam+I*K14)%*%(t(Z)%*%y)
          q14=q.num14/q.den14
          q14=c(q14)
          alphahatk=q14*(solve(t(Z)%*%(Z)+I*K14)%*%t(Z)%*%y)
          alphahatk=c(alphahatk)
          alphahatk_M14[,i] = alphahatk
          
          
         }     #Loop for N=2000 Runs ends
        
        MSE_ols=sum((betahat_ols-c(beta))^2)/N
        MSE1<-sum((alphahatk_M1-c(alpha))^2)/N
        MSE2<-sum((alphahatk_M2-c(alpha))^2)/N
        MSE3<-sum((alphahatk_M3-c(alpha))^2)/N
        MSE4<-sum((alphahatk_M4-c(alpha))^2)/N
        MSE5<-sum((alphahatk_M5-c(alpha))^2)/N
        MSE6<-sum((alphahatk_M6-c(alpha))^2)/N
        MSE7<-sum((alphahatk_M7-c(alpha))^2)/N
        MSE8<-sum((alphahatk_M8-c(alpha))^2)/N
        MSE9<-sum((alphahatk_M9-c(alpha))^2)/N
        MSE10<-sum((alphahatk_M10-c(alpha))^2)/N
        MSE11<-sum((alphahatk_M11-c(alpha))^2)/N
        MSE12<-sum((alphahatk_M12-c(alpha))^2)/N
        MSE13<-sum((alphahatk_M13-c(alpha))^2)/N
        MSE14<-sum((alphahatk_M14-c(alpha))^2)/N
       
        
        MSE=c(MSE_ols,MSE1,MSE2,MSE3,MSE4,MSE5,MSE6, MSE7,MSE8,MSE9,MSE10,MSE11,
              MSE12,MSE13,MSE14)
        
        MSE=round(MSE,10)
        as.matrix(MSE)
        MSE.r[,q]=MSE
      } #End of correlation loop
      MSE.s[,,s]=MSE.r
    } #End of error variance loop
    MSE.n[,,,a]=MSE.s
  }  #End of sample size loop
  MSE.p[,,,,b]=MSE.n
}   #End of predictors loop
col.names=c("0.85","0.95","0.99","0.999")
row.names=c("OLS","HK","HKB","KAM","KGM","KMED","KMS1","LC","ST","QTPKL1","QTPKL2","QTPKL3",
            "QTPKL4","QTPKL5","QTPKL6")
matrix.names=c("0.1","1","5","10")
matrix.names2=c("20","60","120")
matrix.names3=c("4","10")
dimnames(MSE.p)=list(row.names,col.names,matrix.names,matrix.names2,matrix.names3)
write.csv(MSE.p,file = "D:/mphilwork/R-results/Naveed.csv")

