#Compare q (number of source components)

#####Simulated Data Generation####
######sim1#####
q=5  
p=30
n=512
sd=0.1
hz=c(5,10,20,40,70)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=-sin(pi/hz[5]*c(1:n))+wn

S1=cbind(s1,s2,s3,s4,s5)
plot.ts(S1,nc=1,main="sim1")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X1<-S1 %*% A

######sim2#####
q=10  
p=30
n=512
sd=0.1
hz=c(5,10,20,40,50,70,80,100,120,140)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=sin(pi/hz[5]*c(1:n))+wn
s6=-cos(pi/hz[6]*c(1:n))+wn
s7=-sin(pi/hz[7]*c(1:n))+wn
s8=1.5*cos(pi/hz[8]*c(1:n))+wn
s9=-0.5*sin(pi/hz[9]*c(1:n))+wn
s10=-2*sin(pi/hz[10]*c(1:n))+wn

S2=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
plot.ts(S2,nc=1,main="sim2")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X2<-S2 %*% A




#Compare p (number of observed signals)

######sim3#####
q=5  
p=100
n=512
sd=0.1
hz=c(5,10,20,40,70)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=-sin(pi/hz[5]*c(1:n))+wn

S3=cbind(s1,s2,s3,s4,s5)
plot.ts(S3,nc=1,main="sim3")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X3<-S3 %*% A




#Compare n (signal length)

######sim4#####
q=5  
p=100
n=1024
sd=0.1
hz=c(5,10,20,40,70)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=-sin(pi/hz[5]*c(1:n))+wn

S4=cbind(s1,s2,s3,s4,s5)
plot.ts(S4,nc=1,main="sim4")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X4<-S4 %*% A




#Compare sd (standard deviation of white noise)

######sim5#####
q=5  
p=100
n=1024
sd=1
hz=c(5,10,20,40,70)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=-sin(pi/hz[5]*c(1:n))+wn

S5=cbind(s1,s2,s3,s4,s5)
plot.ts(S5,nc=1,main="sim5")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X5<-S5 %*% A



#Compare rng (similarity of source components)

######sim6#####
q=10  
p=30
n=512
sd=0.1
hz=c(5,6,10,12,20,21,40,42,50,51)
#source 5Hz, 10Hz,20Hz,40Hz,70Hz,WN(sd=0.1)
set.seed(1)
wn=rnorm(n,0,sd)
s1=2*sin(pi/hz[1]*c(1:n))+wn
s2=2*cos(pi/hz[2]*c(1:n))+wn
s3=cos(pi/hz[3]*c(1:n))+wn
s4=3*sin(pi/hz[4]*c(1:n))+wn
s5=sin(pi/hz[5]*c(1:n))+wn
s6=-cos(pi/hz[6]*c(1:n))+wn
s7=-sin(pi/hz[7]*c(1:n))+wn
s8=1.5*cos(pi/hz[8]*c(1:n))+wn
s9=-0.5*sin(pi/hz[9]*c(1:n))+wn
s10=-2*sin(pi/hz[10]*c(1:n))+wn

S6=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
plot.ts(S6,nc=1,main="sim6")
#Mixing matrix A
set.seed(1)
A<-matrix(rnorm(p*q,0,1),q,p)
#Mix two signals with the Mixing Matrix(A):
X6<-S6 %*% A

####proposed method####

#####improved ICAbyblock#####
library(ggpubr) #coef fuc
ICAbyBlock_fastica<-function(name,X,maxic){
  library(ggpubr) 
  library(fastICA)
  icabyblock1<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock2<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock3<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  
  for (rep in 1:10) {
    set.seed(rep)
    smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
    grp1<-X[,smp]
    grp2<-X[,-smp]
    for (nic in 1:maxic) {
      grp1.s <- fastICA(grp1, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,
                        method = "C", row.norm = FALSE, maxit = 500,
                        tol = 1e-06, verbose = TRUE)$S
      grp2.s <- fastICA(grp2, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,
                        method = "C", row.norm = FALSE, maxit = 500,
                        tol = 1e-06, verbose = TRUE)$S
      icabyblock1[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="pearson")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock2[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="kendall")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock3[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="spearman")[1:nic,(nic+1):(2*nic)])),2,max))
      
    }
  }
  
  pdf(file = paste0(name,"_fasticacoef_new.pdf"),height = 4,width = 12)
  par(mfrow=c(1,3))
  
  plot(icabyblock1[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="fastICA_pearson")
  for (i in 2:10) {
    lines(icabyblock1[i,],type="l",col=i)
  }
  
  plot(icabyblock2[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="fastICA_kendall")
  for (i in 2:10) {
    lines(icabyblock2[i,],type="l",col=i)
  }
  
  plot(icabyblock3[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="fastICA_spearman")
  for (i in 2:10) {
    lines(icabyblock3[i,],type="l",col=i)
  }
  
  dev.off()

}


ICAbyBlock_fastica("X1",X1,15)
ICAbyBlock_fastica("X2",X2,15)
ICAbyBlock_fastica("X3",X3,15)
ICAbyBlock_fastica("X4",X4,15)
ICAbyBlock_fastica("X5",X5,15)
ICAbyBlock_fastica("X6",X6,15)




ICAbyBlock_imax<-function(name,X,maxic){
  library(steadyICA)
  icabyblock1<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock2<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock3<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  for (rep in 1:10) {
    set.seed(rep)
    smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
    grp1<-X[,smp]
    grp2<-X[,-smp]
    for (nic in 2:maxic) {
      grp1.s <- infomaxICA(grp1, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
      grp2.s <- infomaxICA(grp2, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
      icabyblock1[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="pearson")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock2[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="kendall")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock3[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="spearman")[1:nic,(nic+1):(2*nic)])),2,max))
      
    }
  }
  
  pdf(file = paste0(name,"_infomaxcoef_new.pdf"),height = 4,width = 12)
  par(mfrow=c(1,3))
  
  plot(x=c(2:ncol(icabyblock1)),y=icabyblock1[1,2:ncol(icabyblock1)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="infomax_pearson")
  for (i in 2:10) {
    lines(x=c(2:ncol(icabyblock1)),y=icabyblock1[i,2:ncol(icabyblock1)],type="l",col=i)
  }
  
  plot(x=c(2:ncol(icabyblock2)),y=icabyblock2[1,2:ncol(icabyblock2)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="infomax_kendall")
  for (i in 2:10) {
    lines(x=c(2:ncol(icabyblock2)),y=icabyblock1[i,2:ncol(icabyblock2)],type="l",col=i)
  }
  
  plot(x=c(2:ncol(icabyblock3)),y=icabyblock2[1,2:ncol(icabyblock3)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="infomax_spearman")
  for (i in 2:10) {
    lines(x=c(2:ncol(icabyblock3)),y=icabyblock2[i,2:ncol(icabyblock3)],type="l",col=i)
  }
  
  dev.off()

}

ICAbyBlock_imax("X1",X1,15)
ICAbyBlock_imax("X2",X2,15)
ICAbyBlock_imax("X3",X3,15)
ICAbyBlock_imax("X4",X4,15)
ICAbyBlock_imax("X5",X5,15)
ICAbyBlock_imax("X6",X6,15)




ICAbyBlock_JADE<-function(X,maxic){
  library(JADE)
  icabyblock1<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock2<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  icabyblock3<-matrix(0,nrow=10,ncol=min(maxic,ncol(X)/2))
  for (rep in 1:10) {
    set.seed(rep*2)
    smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
    grp1<-X[,smp]
    grp2<-X[,-smp]
    for (nic in 1:(maxic)) {
      grp1.s <- JADE(grp1,nic,eps=1e-06,maxiter = 500)$S
      grp2.s <- JADE(grp2,nic,eps=1e-06,maxiter = 500)$S
      icabyblock1[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="pearson")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock2[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="kendall")[1:nic,(nic+1):(2*nic)])),2,max))
      icabyblock3[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="spearman")[1:nic,(nic+1):(2*nic)])),2,max))
    }
  }
  
  par(mfrow=c(1,1))
  plot(icabyblock1[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="JADE_pearson")
  for (i in 2:10) {
    lines(icabyblock1[i,],type="l",col=i)
  }
  
  plot(icabyblock2[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="JADE_kendall")
  for (i in 2:10) {
    lines(icabyblock2[i,],type="l",col=i)
  }
  
  plot(icabyblock3[1,],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main="JADE_spearman")
  for (i in 2:10) {
    lines(icabyblock3[i,],type="l",col=i)
  }
  
  return(list(pearson=icabyblock1,kendall=icabyblock2,spearman=icabyblock3))
}

ICAbyBlock_JADE(X6,10)


