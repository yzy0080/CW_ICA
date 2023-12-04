#After generating simulated component signals S and mixed EEG signals X.

######Coef Comparion test######
maxic=15
library(ggpubr) 
library(fastICA)
library(steadyICA)
###ICAbyBlock+Pearson###
icabyblock<-matrix(NA,nrow=maxic,ncol=min(maxic,ncol(X)/2))
set.seed(1)
smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
grp1<-X[,smp]
grp2<-X[,-smp]
for (nic in 2:maxic) {
  grp1.s <- fastICA(grp1, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
  grp2.s <- fastICA(grp2, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
  grp1.s <- infomaxICA(grp1, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
  grp2.s <- infomaxICA(grp2, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
  corr = sort(abs(as.vector(cor(cbind(grp1.s,grp2.s),method="pearson"))),decreasing = T)
  icabyblock[nic,(1:nic)] <- corr[seq((nic*2+2),(nic*2+nic*(2^2-2)),by = 2)]
}
par(mfrow=c(1,1))
plot(x=c(2:ncol(icabyblock)),y=icabyblock[1,2:ncol(icabyblock)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main=paste0("ICAbyBlocks_Pearson_Infomax"))
for (i in 2:maxic) {
  lines(x=c(2:ncol(icabyblock)),y=icabyblock[i,2:ncol(icabyblock)],type="l",col=i)
}

###ICAbyBlock+Spearman###
icabyblock<-matrix(NA,nrow=maxic,ncol=min(maxic,ncol(X)/2))
set.seed(1)
smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
grp1<-X[,smp]
grp2<-X[,-smp]
for (nic in 2:maxic) {
  grp1.s <- fastICA(grp1, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
  grp2.s <- fastICA(grp2, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
  grp1.s <- infomaxICA(grp1, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
  grp2.s <- infomaxICA(grp2, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
  corr = sort(abs(as.vector(cor(cbind(grp1.s,grp2.s),method="spearman"))),decreasing = T)
  icabyblock[nic,(1:nic)] <- corr[seq((nic*2+2),(nic*2+nic*(2^2-2)),by = 2)]
}
par(mfrow=c(1,1))
plot(x=c(2:ncol(icabyblock)),y=icabyblock[1,2:ncol(icabyblock)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main=paste0("ICAbyBlocks_Spearman_fastICA"))
for (i in 2:maxic) {
  lines(x=c(2:ncol(icabyblock)),y=icabyblock[i,2:ncol(icabyblock)],type="l",col=i)
}

###CW_ICA+Pearson###
replicate=10
icabyblock<-matrix(0,nrow=replicate,ncol=min(maxic,ncol(X)/2))
for (rep in 1:replicate) {
  set.seed(rep)
  smp<-sample(c(1:ncol(X)),size = ncol(X)/2)
  grp1<-X[,smp]
  grp2<-X[,-smp]
  for (nic in 2:maxic) {
    grp1.s <- fastICA(grp1, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
    grp2.s <- fastICA(grp2, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
    grp1.s <- infomaxICA(grp1, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
    grp2.s <- infomaxICA(grp2, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
    icabyblock[rep,nic] <- min(apply(as.matrix(abs(cor(cbind(grp1.s,grp2.s),method="pearson")[1:nic,(nic+1):(2*nic)])),2,max))
  }
}

#icabyblock<-icabyblock[,-1]
par(mfrow=c(1,1))
plot(x=c(2:ncol(icabyblock)),y=icabyblock[1,2:ncol(icabyblock)],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main=paste0("CWICA_Pearson_Infomax"))
for (i in 2:replicate) {
  lines(x=c(2:ncol(icabyblock)),y=icabyblock[i,2:ncol(icabyblock)],type="l",col=i)
}

###CW_ICA+Spearman###
replicate=10
CWICA_result = CWICA(replicate,X,maxic,"fastICA",plt=T)


######Accuracy test######
#Start method run
replicate<-c(5,10,25)
method = c("fastICA","Infomax")
result<-NULL

for (i in 1:length(replicate)) {
  accu<-rep(0,replicate[i])
  for (j in 1:length(method)) {
    for (rep in 1:replicate[i]) {
      set.seed(rep)
      A<-matrix(rnorm(p*q,0,1),q,p)
      X<-S %*% A
      ICAbyBlock(X,15,method[j],plt=T)
      DW(X,15,method[j],plt=T)
      accu[rep]<-CWICA(10,X,15,method[j],plt=T)
    }
   }
}

CWICA(10,X,15,method[1],plt=T)
ICAbyBlock(X,15,method[1],plt=T)
DW(X,15,method[1],plt=T)

#ICA_corr_y accuracy
for (r in 1:25) {
  set.seed(r)
  A<-matrix(rnorm(p*q,0,1),q,p)
  X<-S %*% A
  repetation = 30
  result = matrix(0,nrow = maxic-1,ncol = repetation)
  for (rep in 1:repetation){
    for (nic in 2:maxic) {
      set.seed(rep)
      #est.source <- fastICA(X, nic, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "C", row.norm = FALSE, maxit = 500,tol = 1e-06, verbose = TRUE)$S
      est.source <- infomaxICA(X, nic, whiten =TRUE, maxit=500, eps=1e-06)$S
      cor.result = rep(0,nic)
      for (s in 1:nic) {
        cor.result[s]= cor(x=s8,y=est.source[,s],method="pearson")
      }
      result[(nic-1),rep] = max(cor.result)
    }
  }
  par(mfrow=c(1,1))
  plot(x=c(2:maxic),y=result[,1],type="l",col=1,ylim=c(0,1),xlab="number of ICs",ylab="correlation",main=paste0("ICAcorry_Pearson_infomax_rep",r))
  for (i in 2:repetation) {
    lines(x=c(2:maxic),y=result[,i],type="l",col=i)
  }
}



library(readxl)
accuracy <- read_excel("Summary.xlsx",sheet = "accuracySimEEG")
accuracy$`Number of Dataset`<-as.factor(accuracy$`Number of Dataset`)
fast<-accuracy[which(accuracy$`ICA Method`=="fastICA"),]
info<-accuracy[which(accuracy$`ICA Method`!="fastICA"),]
colnames(fast)<-c("Numdata","ICA","Determine","Accuracy")
colnames(info)<-c("Numdata","ICA","Determine","Accuracy")

# library
library(ggplot2)

# Grouped
#FastICA
ggplot(fast, aes(fill=Determine, y=Accuracy, x=Numdata)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
        axis.title = element_text(size = (15)), axis.text= element_text(size = (15)))

#Infomax
ggplot(info, aes(fill=Determine, y=Accuracy, x=Numdata)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
        axis.title = element_text(size = (15)), axis.text= element_text(size = (15)))



#change mix signal parameters


#### p vs q  5 source with 20~100 mixtures ####
#### n vs q  5 source, 30 mix with 64~2048 length ####
#### freq vs q  5 source, 30 mix with 5~500hz ####
#### snr vs q  5 source, 30 mix with 0.5~50 ####
#Robustness plot
  
  #Insert table file
  library(re
  #pvsq
  pvsq <- read_excel("Summary.xlsx", sheet = "p")
  #fastICA
  data <- pvsq[which(pvsq$Algorithm=="FastICA"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:11]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(pvsq)]), 
                            as.numeric(data[2,3:ncol(pvsq)]), 
                            as.numeric(data[3,3:ncol(pvsq)]),
                            as.numeric(data[4,3:ncol(pvsq)])),
                      method = c(as.character(rep(pvsq[1,2], 9)),as.character(rep(pvsq[2,2], 9)),as.character(rep(pvsq[3,2], 9)),as.character(rep(pvsq[4,2], 9))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+
    ggtitle("p vs q s FastICA")+
    xlab("Number of mixed signals")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  #Infomax
  data <- pvsq[which(pvsq$Algorithm=="Infomax"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:11]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(pvsq)]), 
                            as.numeric(data[2,3:ncol(pvsq)]), 
                            as.numeric(data[3,3:ncol(pvsq)]),
                            as.numeric(data[4,3:ncol(pvsq)])),
                      method = c(as.character(rep(pvsq[1,2], 9)),as.character(rep(pvsq[2,2], 9)),as.character(rep(pvsq[3,2], 9)),as.character(rep(pvsq[4,2], 9))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+
    ggtitle("p vs q Infomax")+
    xlab("Number of mixed signals")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  
  #snrvsq
  snrvsq <- read_excel("Summary.xlsx", sheet = "snr")
  #fastICA
  data <- snrvsq[which(snrvsq$Algorithm=="FastICA"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 8)),
                                 as.character(rep(data[2,2], 8)),
                                 as.character(rep(data[3,2], 8)),
                                 as.character(rep(data[4,2], 8))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+
    ggtitle("sd vs q FastICA")+
    xlab("Signal ot Noise Ratio")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  #Infomax
  data <- snrvsq[which(snrvsq$Algorithm=="Infomax"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 8)),
                                 as.character(rep(data[2,2], 8)),
                                 as.character(rep(data[3,2], 8)),
                                 as.character(rep(data[4,2], 8))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+
    ggtitle("snr vs q Infomax")+
    xlab("Signal ot Noise Ratio")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  #nvsq
  nvsq <- read_excel("Summary.xlsx", sheet = "n")
  #fastICA
  data <- nvsq[which(nvsq$Algorithm=="FastICA"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 6)),
                                 as.character(rep(data[2,2], 6)),
                                 as.character(rep(data[3,2], 6)),
                                 as.character(rep(data[4,2], 6))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+
    ggtitle("n vs q FastICA")+
    xlab("Length of signal")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  #Infomax
  data <- nvsq[which(nvsq$Algorithm=="Infomax"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 6)),
                                 as.character(rep(data[2,2], 6)),
                                 as.character(rep(data[3,2], 6)),
                                 as.character(rep(data[4,2], 6))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+  
    ggtitle("n vs q Infomax")+
    xlab("Length of signal")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  #range
  rangevsq <- read_excel("Summary.xlsx", sheet = "range")
  #fastICA
  data <- rangevsq[which(rangevsq$Algorithm=="FastICA"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 7)),
                                 as.character(rep(data[2,2], 7)),
                                 as.character(rep(data[3,2], 7)),
                                 as.character(rep(data[4,2], 7))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+  
    ggtitle("range vs q s FastICA")+
    xlab("Range of Frequency")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  #Infomax
  data <- rangevsq[which(rangevsq$Algorithm=="Infomax"),]
  dataf <- data.frame(x = rep(as.numeric(colnames(data)[3:ncol(data)]),4),                    # Create data frame 
                      y = c(as.numeric(data[1,3:ncol(data)]), 
                            as.numeric(data[2,3:ncol(data)]),
                            as.numeric(data[3,3:ncol(data)]),
                            as.numeric(data[4,3:ncol(data)])),
                      method = c(as.character(rep(data[1,2], 7)),
                                 as.character(rep(data[2,2], 7)),
                                 as.character(rep(data[3,2], 7)),
                                 as.character(rep(data[4,2], 7))))
  ggplot(dataf, aes(x = x, y = y, col = method)) +           # Draw line plot with ggplot2
    geom_line(linewidth=2,aes(linetype=method, color=method))+
    geom_point(aes(shape=method),size=5)+
    scale_shape_manual(values=c(15,16,17,18))+  
    ggtitle("range vs q s FastICA")+
    xlab("Range of Frequency")+
    ylab("Estimated number of source signals")+
    ylim(0,15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black",arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
          axis.title = element_text(size = (10)), axis.text= element_text(size = (20)))
  
  
  
