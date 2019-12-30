rm(list=ls())
source('./STV_linear_functions.R')

TP_all = NULL
TN_all = NULL
##################
n = 200
TPN_all = NULL
for(n in c(200,500,1000,2000,5000,10000)){
  setwd(paste("./n",n,sep=""))
  prefix =  "STV"
  STV_zero1_all = read.table(paste(prefix,"_zero1.txt",sep="")) 
  STV_zero2_all = read.table(paste(prefix,"_zero2.txt",sep="")) 
  STV_zero1_mat1 = matrix(STV_zero1_all[,1], nr = 4)
  STV_zero1_mat2 = matrix(STV_zero1_all[,2], nr = 4)
  STV_zero1_mat3 = matrix(STV_zero1_all[,3], nr = 4)
  STV_zero2_mat1 = matrix(STV_zero2_all[,1], nr = 4)
  STV_zero2_mat2 = matrix(STV_zero2_all[,2], nr = 4)
  STV_zero2_mat3 = matrix(STV_zero2_all[,3], nr = 4)
  
  TP_mat1 = rbind(cbind(STV_zero1_mat1[1,],rep(1,ncol(STV_zero1_mat1))),
                  cbind(STV_zero1_mat2[1,],rep(2,ncol(STV_zero1_mat2))),
                  cbind(STV_zero1_mat3[1,],rep(3,ncol(STV_zero1_mat3))))
  TP_mat1 = as.data.frame(TP_mat1)
  colnames(TP_mat1) = c("TPR","beta")
  TP_mat1$n = rep(n,nrow(TP_mat1))
  TP_mat1$Definition = rep("ETPR",nrow(TP_mat1))
  
  TP_mat2 = rbind(cbind(STV_zero2_mat1[1,],rep(1,ncol(STV_zero2_mat1))),
                  cbind(STV_zero2_mat2[1,],rep(2,ncol(STV_zero2_mat2))),
                  cbind(STV_zero2_mat3[1,],rep(3,ncol(STV_zero2_mat3))))
  TP_mat2 = as.data.frame(TP_mat2)
  colnames(TP_mat2) = c("TPR","beta")
  TP_mat2$n = rep(n,nrow(TP_mat2))
  TP_mat2$Definition = rep("ITPR",nrow(TP_mat2))
  TP_all = rbind(TP_all,TP_mat1,TP_mat2)
  
  TN_mat1 = rbind(cbind(STV_zero1_mat1[2,],rep(1,ncol(STV_zero1_mat1))),
                  cbind(STV_zero1_mat2[2,],rep(2,ncol(STV_zero1_mat2))),
                  cbind(STV_zero1_mat3[2,],rep(3,ncol(STV_zero1_mat3))))
  TN_mat1 = as.data.frame(TN_mat1)
  colnames(TN_mat1) = c("TNR","beta")
  TN_mat1$n = rep(n,nrow(TN_mat1))
  TN_mat1$Definition = rep("ETNR",nrow(TN_mat1))
  
  TN_mat2 = rbind(cbind(STV_zero2_mat1[2,],rep(1,ncol(STV_zero2_mat1))),
                  cbind(STV_zero2_mat2[2,],rep(2,ncol(STV_zero2_mat2))),
                  cbind(STV_zero2_mat3[2,],rep(3,ncol(STV_zero2_mat3))))
  TN_mat2 = as.data.frame(TN_mat2)
  colnames(TN_mat2) = c("TNR","beta")
  TN_mat2$n = rep(n,nrow(TN_mat2))
  TN_mat2$Definition = rep("ITNR",nrow(TN_mat2))
  
  TN_all = rbind(TN_all,TN_mat1,TN_mat2)
  
  
  #####
  TPN_mat1 = rbind(cbind(t(STV_zero1_mat1[1:2,]),rep(1,ncol(STV_zero1_mat1))),
                  cbind(t(STV_zero1_mat2[1:2,]),rep(2,ncol(STV_zero1_mat2))),
                  cbind(t(STV_zero1_mat3[1:2,]),rep(3,ncol(STV_zero1_mat3))))
  TPN_mat1 = as.data.frame(TPN_mat1)
  colnames(TPN_mat1) = c("TPR","TNR","beta")
  TPN_mat1$n = rep(n,nrow(TPN_mat1))
  TPN_mat1$Definition = rep("Eestimation",nrow(TPN_mat1))
  
  TPN_mat2 = rbind(cbind(t(STV_zero2_mat1[1:2,]),rep(1,ncol(STV_zero2_mat1))),
                  cbind(t(STV_zero2_mat2[1:2,]),rep(2,ncol(STV_zero2_mat2))),
                  cbind(t(STV_zero2_mat3[1:2,]),rep(3,ncol(STV_zero2_mat3))))
  TPN_mat2 = as.data.frame(TPN_mat2)
  colnames(TPN_mat2) = c("TPR","TNR","beta")
  TPN_mat2$n = rep(n,nrow(TPN_mat2))
  TPN_mat2$Definition = rep("Inference",nrow(TPN_mat2))
  
  TPN_all = rbind(TPN_all,TPN_mat1,TPN_mat2)
}


##################
setwd(".")
save(TP_all,TN_all,file="summary_data.Rdata")
######################## end of the code ##############

library(ggplot2)

TP_all$n = factor(TP_all$n)
ggplot(TP_all, aes(x=n, y=TPR, fill=Definition)) + 
  geom_boxplot()+facet_wrap(~beta, scale="free")+ scale_fill_grey() + theme_classic()

TN_all$n = factor(TN_all$n)
ggplot(TN_all, aes(x=n, y=TNR, fill=Definition)) + 
  geom_boxplot()+facet_wrap(~beta, scale="free")+scale_fill_grey() + theme_classic()
for(j in 1:3){
  pdf(file = paste("TPR_beta",j,".pdf",sep=""),width = 6,height = 6)
  print(ggplot(TP_all[TP_all$beta==j,], aes(x=n, y=TPR, fill=Definition)) + 
          geom_boxplot()+ scale_fill_grey() + theme_classic()+ylim(0,1)
  )
  dev.off()
  pdf(file = paste("TNR_beta",j,".pdf",sep=""),width = 6,height = 6)
  print(ggplot(TN_all[TN_all$beta==j,], aes(x=n, y=TNR, fill=Definition)) + 
          geom_boxplot()+ scale_fill_grey() + theme_classic()+ylim(0,1))
  dev.off()
}


library(plyr)    
TPR.sum<-ddply(TP_all,.(n,Definition),
                    summarize, value = mean(TPR))

ggplot(TPR.sum, aes(x=n, y=value, col=Definition)) + geom_point()+scale_fill_grey() + theme_classic()+ylim(0,1)
for(j in 1:3){
  TPR.sum<-ddply(TP_all[TP_all$beta==j,],.(n,Definition),
                 summarize, value = mean(TPR))
  TNR.sum<-ddply(TN_all[TN_all$beta==j,],.(n,Definition),
                 summarize, value = mean(TNR))
  pdf(file = paste("newTPR_beta",j,".pdf",sep=""),width = 6,height = 6)
  print(ggplot(TPR.sum, aes(x=n, y=value, col=Definition)) + geom_point()+scale_fill_grey() + theme_classic()+ylim(0,1)
  )
  dev.off()
  pdf(file = paste("newTNR_beta",j,".pdf",sep=""),width = 6,height = 6)
  print(ggplot(TNR.sum, aes(x=n, y=value, col=Definition)) + geom_point()+scale_fill_grey() + theme_classic()+ylim(0,1))
  dev.off()
}


### summary table

library(plyr)    
TPR.mean<-ddply(TP_all,.(n,Definition,beta),
               summarize, TPRmean = mean(TPR), TPRsd=sd(TPR))

TNR.mean<-ddply(TN_all,.(n,Definition,beta),
                summarize, TNRmean = mean(TNR),TNRsd= sd(TNR))

TPN.mean <- ddply(TPN_all,.(n,Definition,beta),
                  summarize, TPRmean = mean(TPR), TPRsd=sd(TPR),TNRmean = mean(TNR),TNRsd= sd(TNR))


fun.aggregate = function(x) {
  return(sprintf('%0.0f (%0.0f)', mean(x)*1000, sd(x)*1000))
}
TPN.mean <- ddply(TPN_all,.(n,Definition,beta),
                  summarize, TPR = fun.aggregate(TPR), TNR=fun.aggregate(TNR))

j = 1

table_final = NULL
for(j in 1:3){
  dat1 = TPN.mean[TPN.mean$beta==j&TPN.mean$Definition=="ETPNR",]
  dat2 = TPN.mean[TPN.mean$beta==j&TPN.mean$Definition=="ITPNR",]
  table_final = rbind(table_final,dat1$TPR,dat2$TPR,dat1$TNR,dat2$TNR)
}
table_final = cbind(c(" ","$b_1$"," "," "," ","$b_2$"," "," "," ","$b_3$"," "," "),
                    rep(c("ETPR","ITPR","ETNR","ITNR"),3),
                    table_final)

colnames(table_final) = c("beta","measure","200","500","1000","2000","5000","10000")
rownames(table_final) = NULL
sink(file = 'latex_tbl.txt')
print(xtable::xtable(table_final), include.rownames=FALSE)
sink()