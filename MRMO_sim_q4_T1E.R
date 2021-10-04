library(GGally)
library(beepr)

library(MRMO)

setwd("/Users/dyq/downloads")


n=559
maf=readRDS("maf.RDS")

p=19
q=4

############## 2+2
if(q==4){
  q1=2; q2=2
  sss=c(1,1,2,3)
  a=c(-1,-1.5,0,0)
  b=c(0,0,0,0)    #T1E
  #b=c(0,0.4,0,0)       #power scenario 1
  #b=c(0.3,0.3,0,0)     #2
  #b=c(0.3,0,0,0.3)     #3
  #b=c(0.2,0.3,0.3,0)      #4
  #b=c(0.25,0.25,0.25,0.25)    #5
}

#cor
rho1=0.5
rho2=0.3 #or -0.3      (additional settings in the supplementary materials: rho1=rho2=0; rho1=0.6, rho2=-0.6)
if(q==4){
  corZ=matrix(rho1,nrow=q,ncol=q)
  corZ[1:3,4]=corZ[4,1:3]=rho2
  diag(corZ)=1
}

which_bin=1:q1
if(q2==0){
  which_cont=c()
}else{
  which_cont=(q1+1):(q1+q2)
}

iter=5000


size_UV=size_MV=sd_UV=sd_MV=p_UV=p_MV=matrix(nrow=iter,ncol=q)
p_minP=p_Wald=p_aSPU=p_aSPU2=c()


now=Sys.time()

for(it in 1:iter){
  set.seed(it+4000)

  ################ Simulate X
  trunk=0.08
  size_GX=rnorm(p*50,sd=0.15)
  size_GX=size_GX[size_GX>trunk | size_GX< -trunk]
  size_GX=size_GX[1:p]
  size_GU=size_GX*0
  size_UX=runif(1)

  G=t(matrix(rbinom(n*p,2,prob=maf),nrow=p))

  U=G%*%size_GU+rnorm(nrow(G))
  X=G%*%size_GX + size_UX*U + rnorm(nrow(G))

  #Generate X (effect sizes are modified to control proportion of X explained by G)
  explained=0.2  #proportion of X explained by G
  if(explained>0){
    sdga=sd(G%*%size_GX)
    sii=sqrt(2*explained/(1-explained))/sdga
    size_GX=size_GX*sii
  }
  X=G%*%size_GX + size_UX*U + rnorm(nrow(G))
  var(G%*%size_GX)/var(X)


  ################ Simulate Y

  covZ=diag(sss)%*%corZ%*%diag(sss)                       #cov(Z)

  meanZ=matrix(rnorm(n*q),ncol=q)*0
  for(i in 1:q){
    meanZ[,i]=a[i]+b[i]*X
  }                                                        #mean(Z)

  size_UY=runif(q)

  Z=matrix(rnorm(n*q),ncol=q)
  Z=t(t(chol(covZ))%*%t(Z))+meanZ+U%*%t(size_UY)

  #probit link
  po=pnorm(Z[,which_bin])
  Y=rbinom(length(po),1,c(po))
  Y=matrix(Y,ncol=q1)

  Y=cbind(Y,Z[,which_cont])
  colMeans(Y)

  #analysis

  #get Xhat
  X1=X
  G1=G
  ll1=lm(X~G)
  Xhat0=ll1$fitted.values

  #2nd sample
  G2=t(matrix(rbinom(n*p,2,prob=maf),nrow=p))
  U2=G2%*%size_GU+rnorm(nrow(G2))
  X2=G2%*%size_GX + size_UX*U2 + rnorm(nrow(G2))
  meanZ2=matrix(rnorm(n*q),ncol=q)*0
  for(i in 1:q){
    meanZ2[,i]=a[i]+b[i]*X2
  }                                                        #mean(Z)
  Z2=matrix(rnorm(n*q),ncol=q)
  Z2=t(t(chol(covZ))%*%t(Z2))+meanZ2+U2%*%t(size_UY)
  #probit link
  po2=pnorm(Z2[,which_bin])
  Y2=rbinom(length(po2),1,c(po2))
  Y2=matrix(Y2,ncol=q1)
  Y2=cbind(Y2,Z2[,which_cont])

  GG2=cbind(rep(1,nrow(G2)),G2)
  Xhat=c(GG2%*%ll1$coefficients)

  #UVA, MCLE
  Ys=Y2
  N=n
  hoho=MRMO(G1,X1,G2,Y2,q1,q2,newton=1) #use newton=0 for gradient descent

  #save
  size_UV[it,]=hoho[1:q,1]
  sd_UV[it,]=hoho[1:q,2]
  p_UV[it,]=hoho[1:q,3]
  size_MV[it,]=hoho[1:q,4]
  sd_MV[it,]=hoho[1:q,5]
  p_MV[it,]=hoho[1:q,6]

  p_minP[it]=hoho[q+1,1]
  p_Wald[it]=hoho[q+1,2]
  p_aSPU[it]=hoho[q+1,3]   #permutation based aSPU
  p_aSPU2[it]=hoho[q+1,4]   #distribution based aSPU

  if(it%%10==0){
    now2=Sys.time()
    print(paste(it,"------------",round(now2-now,digits=3)))
  }
}

est=cbind(b,colMeans(size_UV,na.rm=TRUE),colMeans(sd_UV,na.rm=TRUE),colMeans(p_UV<0.05,na.rm=TRUE),colMeans(size_MV,na.rm=TRUE),colMeans(sd_MV,na.rm=TRUE),colMeans(p_MV<0.05,na.rm=TRUE))
colnames(est)=c("Truth","UVA mean","meanSD","rej","MVA mean","meanSD","rej")

est=rbind(est,est[1,])
est[q+1,]=NA
est[q+1,1:3]=ov
est[q+1,4]=mean(p_aSPU2<0.05,na.rm=TRUE)
row.names(est)=c(1:q,"minP/Wald/aSPU rej")

est2=round(est,digits=3)
est2

est2[nrow(est2),ncol(est2)]=sum(is.na(size_MV[,1]))

beep();beep();beep();beep();beep();Sys.sleep(1);beep();beep();beep();beep();beep();Sys.sleep(1);beep();beep();beep();beep();beep()

mincor=round(min(corZ-diag(diag(corZ))),digits=2)
maxcor=round(max(corZ-diag(diag(corZ))),digits=2)
nono=paste("MRMO_sim ",q1,"_",q2," iter",iter,"_n",n,"_cor",mincor,"_",maxcor,"_size",paste(b, collapse = "_"),".csv",sep="")

write.csv(est2,file=nono)
