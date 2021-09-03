setwd("G:/R/1124/KIRC/methy")
#library(TCGAbiolinks)
library(SummarizedExperiment)
library("minfi")
library("impute")
library("wateRmelon")
library(survival)
library(Hmisc)
library(dplyr)
library(plyr)
library(caret)
library(glmnet)
library(ROCR)
library(survminer)
library(survivalROC)
library(rms)
library(foreign)
library(data.table)
TCGA <- as.data.frame(fread("KIRCtop20percentVarCgsiteMatrix.txt.gz"))
rownames(TCGA) <- TCGA[,1]
TCGA[,1] <- NULL

MTAB <- as.data.frame(fread("EMTAB_3247.raw.matrix107/E-MTAB-3274.raw.matrix107.txt.gz"))
rownames(MTAB) <- MTAB[,1]
MTAB[,1] <- NULL

x<-substr(colnames(TCGA),regexpr("TCGA",colnames(TCGA)),regexpr("T",colnames(TCGA))+14)
methy<-TCGA[,unlist(lapply(x,function(x) strsplit(as.character(x),"-")[[1]][4])) == "01"]
colnames(methy) <- substr(colnames(methy),regexpr("TCGA",colnames(methy)),regexpr("T",colnames(methy))+11)

clinSub1 <- read.table("Clinical BCR XML.merge.txt",header = T,stringsAsFactors = F,sep="\t")
clinSub2 <- clinSub1[!is.na(clinSub1$A8_New_Event_Time) & !is.na(clinSub1$A8_New_Event),]

trans <- as.data.frame(t(methy))
rownames(trans) <- unlist(lapply(rownames(trans), function(x) gsub("\\.","-",x)))
merg1<-merge(trans,clinSub2,by.x="row.names",by.y="A0_Samples",sort=F)

info<-as.data.frame(merg1$Row.names)
colnames(info)[1]<-"Sample"  #修改列名
info$Group<-merg1$A8_New_Event

mat=as.matrix(t(merg1[,2:(ncol(trans)+1)]))  #转换成表达矩阵
colnames(mat) <- merg1[,1]
mat=impute.knn(mat)  #这个填充NA数据
matData=mat$data    # 获取填充之后的表达矩阵
matData1=matData+0.00001

matData2=matData1[rowMeans(matData1)>0.005,]  #去除表达值很低的甲基化位点
matData3 = betaqn(matData2)   # 标准化
merg_Beta <- merge(as.data.frame(t(matData3)),clinSub2,by.x="row.names",by.y="A0_Samples",sort=F)
#write.table(matData3,file="norm.xls",sep="\t",quote=F)
grset=makeGenomicRatioSetFromMatrix(matData3,what="Beta")  #中间对象

M = getM(grset)
dmp <- dmpFinder(M, pheno=info$Group, type="categorical")
dmpDiff=dmp[(dmp$pval<0.05) & (is.na(dmp$pval)==F),]
dim(dmpDiff)


############# 外部验证集部分
GSEclin1<-read.table("E-MTAB-3274clinical.txt",header = T,sep = "\t",stringsAsFactors = F)
GSEclin <- GSEclin1[!is.na(GSEclin1$progression) & !is.na(GSEclin1$progression_free_survival),]
GSE <- t(MTAB[rownames(MTAB) %in% rownames(dmpDiff),])[!is.na(GSEclin1$progression) & !is.na(GSEclin1$progression_free_survival),]
GSEmerg <- cbind(GSE,GSEclin)
info<-as.data.frame(GSEmerg$Extract_Name)
colnames(info)[1]<-"Sample"  #修改列名
info$Group<-GSEmerg$progression
info$Sample <- as.character(info$Sample)


mat=as.matrix(t(GSEmerg[,1:ncol(GSE)]))  #转换成表达矩阵
colnames(mat) <- GSEmerg[,"Extract_Name"]
mat=impute.knn(mat)  #这个填充NA数据
matData=mat$data    # 获取填充之后的表达矩阵
matData1=matData+0.00001
matData2=matData1[rowMeans(matData1)>0.005,]  #去除表达值很低的甲基化位点
matData3GSE = betaqn(matData2)   # 标准化
merg_GSE <- cbind(t(matData3GSE),GSEclin)
#write.table(matData3,file="norm.xls",sep="\t",quote=F)
grset=makeGenomicRatioSetFromMatrix(matData3GSE,what="Beta")  #中间对象

M = getM(grset)
dmp <- dmpFinder(M, pheno=info$Group, type="categorical")
dmpDiff1=dmp[(dmp$pval<0.05) & (is.na(dmp$pval)==F),]
dim(dmpDiff1)
###################

cg <- intersect(rownames(dmpDiff),rownames(dmpDiff1))
merg <- merg_Beta
merg$A8_New_Event_Time <- as.numeric(merg$A8_New_Event_Time)
merg$status <- as.numeric(merg$A8_New_Event)
Lusurv <- Surv(time = merg$A8_New_Event_Time,event = merg$status)
merg$Lusurv <- with(merg,Lusurv)

UniCox <- function(x){
    FML <- as.formula(paste0('Lusurv~',x))
    Gcox<-coxph(FML,data = merg)
    Gsum<-summary(Gcox)
    HR <- round(Gsum$coefficients[,2],5)
    Pvalue <- round(Gsum$coefficients[,5],5)
    CI <- paste0(round(Gsum$conf.int[,3:4],5),collapse = "-")
    Unicox <- data.frame('Characteristics' = x,
                         'Hazard Ratio' = HR,
                         'CI95' = CI,
                         'P value' = Pvalue)
    return(Unicox)
}
VarNames <- cg
UniVar <- lapply(VarNames, UniCox)
UniVar_OS <- ldply(UniVar,data.frame)
dim(UniVar_OS)
table(UniVar_OS$P.value<0.05)


p<-0.00025
merg <- merg_Beta[,c("Row.names",cg,colnames(merg_Beta)[grep("A18",colnames(merg_Beta)):ncol(merg_Beta)])]  # 缺
merg$A8_New_Event_Time <- merg$A8_New_Event_Time+0.01
merg_GSE39279 <- merg_GSE[,c(cg,colnames(merg_GSE)[grep("^Extract_Name",colnames(merg_GSE)):ncol(merg_GSE)])]
merg_GSE39279$day <- merg_GSE39279$progression_free_survival*30
AUCmerg <- c()
AUCTrain <- c()
merg_GSE39279AUC <- c()
AUCTest <- c()
for (i in 1:1000) {
    set.seed(i)
    inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
    Train <- merg[ inTrain, ]
    Test <- merg[ -inTrain, ]
    fml <- as.formula(paste0(' ~ ',paste0(UniVar_OS$Characteristics[UniVar_OS$P.value<p],collapse = "+")))
    x <- as.matrix(Train[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])])
    y <- Surv(time = Train$A8_New_Event_Time,event = Train$A8_New_Event)
    cv.fit <- cv.glmnet(x, y, family="cox", alpha=1,type.measure = 'auc')
    Test$lasso.prob <- predict(cv.fit,type="response",  newx = as.matrix(Test[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])]), s = 'lambda.min')
    predTest <- prediction(Test$lasso.prob, Test$A8_New_Event)
    perfTest <- performance(predTest,"tpr","fpr")
    aucTest <- performance(predTest,"auc") # shows calculated AUC for model
    AUCTest <- c(AUCTest,aucTest@y.values[[1]])
	
    Train$lasso.prob <- predict(cv.fit,type="response",  newx = as.matrix(Train[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])]), s = 'lambda.min')
    predTrain <- prediction(Train$lasso.prob, Train$A8_New_Event)
    perfTrain <- performance(predTrain,"tpr","fpr")
    aucTrain <- performance(predTrain,"auc") # shows calculated AUC for model
    AUCTrain <- c(AUCTrain,aucTrain@y.values[[1]])	
    
    merg_GSE39279$lasso.prob <- predict(cv.fit,type="response",  newx = as.matrix(merg_GSE39279[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])]), s = 'lambda.min')
    merg_GSE39279predTrain <- prediction(merg_GSE39279$lasso.prob, merg_GSE39279$progression)
    merg_GSE39279perfTrain <- performance(merg_GSE39279predTrain,"tpr","fpr")
    merg_GSE39279aucTrain <- performance(merg_GSE39279predTrain,"auc") # shows calculated AUC for model
    merg_GSE39279AUC <- c(merg_GSE39279AUC,merg_GSE39279aucTrain@y.values[[1]])
    
    merg$lasso.prob <- predict(cv.fit,type="response",  newx = as.matrix(merg[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])]), s = 'lambda.min')
    predmerg <- prediction(merg$lasso.prob, merg$A8_New_Event)
    perfmerg <- performance(predmerg,"tpr","fpr")
    aucmerg <- performance(predmerg,"auc") # shows calculated AUC for model
    AUCmerg <- c(AUCmerg,aucmerg@y.values[[1]])
	print(i)
}
cb<-do.call("data.frame",list(AUCmerg,AUCTrain,merg_GSE39279AUC,AUCTest))
colnames(cb) <- c("merg","Train","GSE39279","Test")
subset(cb,cb[,1]>0.7 & cb[,2]>0.7 & cb[,3]>0.7 & cb[,4]>0.7)


         # merg     Train  GSE39279      Test
# 45  0.8249303 0.8821510 0.7055985 0.7091946
# 58  0.7933592 0.8411038 0.7215251 0.7083333
# 153 0.8038375 0.8659951 0.7065637 0.7004548
# 317 0.8062190 0.8483660 0.7215251 0.7284722
# 469 0.7851262 0.8033008 0.7075290 0.7546419
# 641 0.8240457 0.8535634 0.7031853 0.7698939
# 872 0.8194870 0.8363844 0.7022201 0.7911618
# 881 0.7879159 0.8010584 0.7022201 0.7576622

p<-0.00025
seed=641
set.seed(seed)
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
fml <- as.formula(paste0(' ~ ',paste0(UniVar_OS$Characteristics[UniVar_OS$P.value<p],collapse = "+")))
x <- as.matrix(Train[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])])
y <- Surv(time = Train$A8_New_Event_Time,event = Train$A8_New_Event)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1,type.measure = 'auc')
Test$lasso.prob <- predict(cv.fit,type="response",  newx = as.matrix(Test[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])]), s = 'lambda.min')
pred <- prediction(Test$lasso.prob, Test$A8_New_Event)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc") # shows calculated AUC for model
coe <- coef(cv.fit, s = "lambda.min")

rt <- merg[,c("Row.names",names(coe[match(coe@x,coe[,1]),]),tail(colnames(merg),20))]
rt[,"A8_New_Event_Time"]=rt[,"A8_New_Event_Time"]/365
rt$Lusurv <- Surv(rt$A8_New_Event_Time, rt$A8_New_Event)
fml <- as.formula(paste0('Lusurv~',paste(names(coe[match(coe@x,coe[,1]),]),collapse = " + ")))
multiCox <- coxph(fml,data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

sco<-c()
for (i in 1:nrow(merg)) {
    risk <- 0
    for (j in 1:length(outTab[,2])) {
        risk<-risk+merg[i,match(rownames(outTab),colnames(merg))[j]] * as.numeric(outTab[,2][j])
    }
    sco<-c(sco,risk)
}
cg_lasso_names <- rownames(outTab)
risk_score_formula <- paste(round(as.numeric(outTab[,2]),5),rownames(outTab),sep = "*",collapse = " + ")
# "2.65802*cg04092800 + 4.53723*cg09217923 + 2.13544*cg14427009 + 1.8933*cg11132272 + 1.68719*cg14202477 + 1.27751*cg12436377 + 2.68546*cg03010887 + 1.96867*cg04094846"

merg$Score <- sco
merg$OS_Time <- merg$A8_New_Event_Time/365
rocCol=c("red","green","blue")
aucText=c()

tiff("OverallROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=merg$OS_Time, status=merg$A8_New_Event, marker = merg$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)

roc=survivalROC(Stime=merg$OS_Time, status=merg$A8_New_Event, marker = merg$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=merg$OS_Time, status=merg$A8_New_Event, marker = merg$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()


set.seed(seed)
p<-0.00025
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]

rocCol=c("red","green","blue")
aucText=c()

tiff("TrainROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=Train$OS_Time, status=Train$A8_New_Event, marker = Train$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)

roc=survivalROC(Stime=Train$OS_Time, status=Train$A8_New_Event, marker = Train$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=Train$OS_Time, status=Train$A8_New_Event, marker = Train$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()

#######################  TestROC  ################
rocCol=c("red","green","blue")
aucText=c()

tiff("TestROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=Test$OS_Time, status=Test$A8_New_Event, marker = Test$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)

roc=survivalROC(Stime=Test$OS_Time, status=Test$A8_New_Event, marker = Test$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=Test$OS_Time, status=Test$A8_New_Event, marker = Test$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()

#######################  GSE39279 ROC  ################
score<-c()
for (i in 1:nrow(merg_GSE39279)) {
    risk <- 0
    for (j in 1:length(coe@x)) {
        risk<-risk+merg_GSE39279[i,match(names(coe[match(coe@x,coe[,1]),]),colnames(merg_GSE39279))[j]] * coe@x[j]
    }
    score<-c(score,risk)
}
merg_GSE39279$Score <- score
merg_GSE39279$time_rec <- merg_GSE39279$progression_free_survival/12

#score<-c()
#for (i in 1:nrow(merg_GSE39279)) {
#    risk <- 0
#    for (j in 1:length(outTab[,2])) {
#        risk<-risk+merg_GSE39279[i,match(rownames(outTab),colnames(merg_GSE39279))[j]] * as.numeric(outTab[,2][j])
#    }
#    score<-c(score,risk)
#}
#merg_GSE39279$Score <- score

rocCol=c("red","green","blue")
aucText=c()

tiff("GSE39279ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=merg_GSE39279$time_rec, status=merg_GSE39279$progression, marker = merg_GSE39279$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)

roc=survivalROC(Stime=merg_GSE39279$time_rec, status=merg_GSE39279$progression, marker = merg_GSE39279$Score, predict.time = 4, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=merg_GSE39279$time_rec, status=merg_GSE39279$progression, marker = merg_GSE39279$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()


##################  overallSurvival  ######################
merg$Risk[merg$Score>=median(sco)]<-"High"
merg$Risk[merg$Score<median(sco)]<-"Low"
merg$Lusurv <- Surv(time = merg$A8_New_Event_Time,event = merg$A8_New_Event)
fit <- survfit(Lusurv~Risk,data = merg)
sdf <- survdiff(Lusurv~Risk,data = merg)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("OverallSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.7,0.2,pp,cex=1.5)
dev.off()

set.seed(seed)   # Train survival
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
fit <- survfit(Lusurv~Risk,data = Train)
sdf <- survdiff(Lusurv~Risk,data = Train)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("TrainSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.5,0.7,pp,cex=1.5)
dev.off()

set.seed(seed)  # Test survival
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
fit <- survfit(Lusurv~Risk,data = Test)
sdf <- survdiff(Lusurv~Risk,data = Test)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("TestSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white",inset=c(-0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.6,0.6,pp,cex=1.5)
dev.off()

################ GSE39279 survival ######################
merg_GSE39279$Risk[merg_GSE39279$Score>=median(score)]<-"High"
merg_GSE39279$Risk[merg_GSE39279$Score<median(score)]<-"Low"
#merg_GSE39279$time_rec <- merg_GSE39279$time_rec*365
merg_GSE39279$Lusurv <- Surv(time = merg_GSE39279$day,event = merg_GSE39279$progression)
fit <- survfit(Lusurv~Risk,data = merg_GSE39279)
sdf <- survdiff(Lusurv~Risk,data = merg_GSE39279)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("GSE39279Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg_GSE39279$day)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("topright",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(-0.02,-0.03))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg_GSE39279$day)*0.8,0.5,pp,cex=1.5)
dev.off()

#####################  S线     #############################
merg$Status[merg$A8_New_Event==0] <- "No recurrence"
merg$Status[merg$A8_New_Event==1] <- "Recurrence"
df<-data.frame(Status=merg$Status,Group=merg$Risk,Score=merg$Score)
dat<-df[order(df$Score),]
dat$Rank<-1:nrow(dat)
tiff("Sline.tiff",units = "in",width = 5,height = 3,res = 600,compression = "lzw")
ggplot(data=dat, aes(Rank, Score, colour = Group))+ 
    geom_point(aes(shape=Group))+ scale_shape_manual(values=c(2, 16)) + scale_color_manual(values=c("red", "green"))+theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black",size = 1.5),
          plot.margin = unit(c(1,1,1,1), "cm"),axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
dev.off()

df<-data.frame(Status=merg$Status,Group=merg$Risk,Score=merg$Score,Time=merg$A8_New_Event_Time)
dat<-df[order(df$Score),]
dat$Rank<-1:nrow(dat)
dat$Status[merg$A8_New_Event==0] <- "No recurrence"
dat$Status[merg$A8_New_Event==1] <- "Recurrence"
tiff("Status_Time.tiff",units = "in",width = 5,height = 3,res = 600,compression = "lzw")
ggplot(data=dat, aes(Rank, Time, colour = Status))+ 
    geom_point(size=1)+ scale_color_manual(values=c("green", "red"))+theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA,size=1.5),panel.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
dev.off()
# 热图
library(pheatmap)
dat<-merg[order(merg$Score),]
heat<-dat[,colnames(dat) %in% cg_lasso_names]
xy <- t(apply(heat, 2, order))
colnames(xy) <- rownames(heat)
annotation_col = data.frame(Group=factor(c(rep("Low_risk",floor(nrow(dat)/2)),rep("High_risk",nrow(dat)-floor(nrow(dat)/2)))))
rownames(annotation_col) <- rownames(heat)
tiff("heatmap.tiff",units = "in",width = 6,height = 2,res = 600,compression = "lzw")
pheatmap(xy,color = colorRampPalette(c("green","black","red"))(100),show_colnames = F,cluster_cols = F,fontsize_row =5,annotation_col=annotation_col,annotation_colors = list(Group = c(High_risk = rgb(244/255,211/255,121/255), Low_risk = rgb(90/255,176/255,249/255))))
dev.off()



merg_GSE39279$Status[merg_GSE39279$progression==0] <- "No recurrence"
merg_GSE39279$Status[merg_GSE39279$progression==1] <- "Recurrence"
df<-data.frame(Status=merg_GSE39279$Status,Group=merg_GSE39279$Risk,Score=merg_GSE39279$Score)
dat<-df[order(df$Score),]
dat$Rank<-1:nrow(dat)
tiff("SlineGSE39279.tiff",units = "in",width = 5,height = 3,res = 600,compression = "lzw")
ggplot(data=dat, aes(Rank, Score, colour = Group))+ 
    geom_point(aes(shape=Group))+ scale_shape_manual(values=c(2, 16)) + scale_color_manual(values=c("red", "green"))+theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black",size = 1.5),
          plot.margin = unit(c(1,1,1,1), "cm"),axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
dev.off()

df<-data.frame(Status=merg_GSE39279$Status,Group=merg_GSE39279$Risk,Score=merg_GSE39279$Score,Time=merg_GSE39279$day)
dat<-df[order(df$Score),]
dat$Rank<-1:nrow(dat)
dat$Status[merg_GSE39279$progression==0] <- "No recurrence"
dat$Status[merg_GSE39279$progression==1] <- "Recurrence"
tiff("Status_TimeGSE39279.tiff",units = "in",width = 5,height = 3,res = 600,compression = "lzw")
ggplot(data=dat, aes(Rank, Time, colour = Status))+ 
    geom_point(size=1)+ scale_color_manual(values=c("green", "red"))+theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          panel.border = element_rect(colour = "black", fill=NA,size=1.5),panel.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
dev.off()

# 热图
library(pheatmap)
dat<-merg_GSE39279[order(merg_GSE39279$Score),]
heat<-dat[,colnames(dat) %in% cg_lasso_names]
xy <- t(apply(heat, 2, order))
colnames(xy) <- rownames(heat)
annotation_col = data.frame(Group=factor(c(rep("Low_risk",floor(nrow(dat)/2)),rep("High_risk",nrow(dat)-floor(nrow(dat)/2)))))
rownames(annotation_col) <- rownames(heat)
tiff("heatmapGSE39279.tiff",units = "in",width = 6,height = 2,res = 600,compression = "lzw")
pheatmap(xy,color = colorRampPalette(c("green","black","red"))(100),show_colnames = F,cluster_cols = F,fontsize_row =5,annotation_col=annotation_col,annotation_colors = list(Group = c(High_risk = rgb(244/255,211/255,121/255), Low_risk = rgb(90/255,176/255,249/255))))
dev.off()


#######################
set.seed(seed)
p <- 0.00025
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
fml <- as.formula(paste0(' ~ ',paste0(UniVar_OS$Characteristics[UniVar_OS$P.value<p],collapse = "+")))
x <- as.matrix(Train[,as.character(subset(UniVar_OS,UniVar_OS[,4]<p)[,1])])
y <- Surv(time = Train$A8_New_Event_Time,event = Train$A8_New_Event)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1,type.measure = 'auc')

tiff("lambda.tiff",units = "in",width = 6,height = 3.8,res = 600,compression = "lzw")
plot(cv.fit$glmnet.fit, "lambda")
box(lwd=3)
dev.off()

tiff("lambdaMin.tiff",units = "in",width = 6,height = 3.8,res = 600,compression = "lzw")
plot(cv.fit)
box(lwd=3)
dev.off()


######################### 临床基线表 #################################
merg$age_at_initial_pathologic_diagnosis <- as.numeric(merg$age_at_initial_pathologic_diagnosis)
merg$age[merg$age_at_initial_pathologic_diagnosis >60 & merg$age_at_initial_pathologic_diagnosis < 150] <- ">60"
merg$age[merg$age_at_initial_pathologic_diagnosis <=60] <- "<=60"
set.seed(seed)
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
choose<- c("A18_Sex","age","A6_Stage","A5_M","A3_T","A4_N","A7_Grade","A9_Cancer_Status","hemoglobin_result","laterality","platelet_qualitative_result","race_list.race","serum_calcium_result","white_cell_count_result")
Base <- subset(merg,select = choose)
Freq <- lapply(Base,table)
Prop <- lapply(Freq, prop.table)
Char <- NULL
for (i in 1:length(Freq)){
    Character <- c(names(Freq[i]),names(Freq[[i]]))
    Noc <- c(NA,paste0(Freq[[i]],"(",round(Prop[[i]],4)*100,")"))
    Characteristics <- data.frame("Characteristics" = Character,"Number of cases(%)" = Noc)
    Char <- rbind(Char,Characteristics)
}

Base <- subset(Train,select = choose)
Freq <- lapply(Base,table)
Prop <- lapply(Freq, prop.table)
CharTrain <- NULL
for (i in 1:length(Freq)){
    Character <- c(names(Freq[i]),names(Freq[[i]]))
    Noc <- c(NA,paste0(Freq[[i]],"(",round(Prop[[i]],4)*100,")"))
    Characteristics <- data.frame("Characteristics" = Character,"Number of cases(%)" = Noc)
    CharTrain <- rbind(CharTrain,Characteristics)
}

Base <- subset(Test,select = choose)
Freq <- lapply(Base,table)
Prop <- lapply(Freq, prop.table)
CharTest <- NULL
for (i in 1:length(Freq)){
    Character <- c(names(Freq[i]),names(Freq[[i]]))
    Noc <- c(NA,paste0(Freq[[i]],"(",round(Prop[[i]],4)*100,")"))
    Characteristics <- data.frame("Characteristics" = Character,"Number of cases(%)" = Noc)
    CharTest <- rbind(CharTest,Characteristics)
}

merg.Train<-left_join(Char,CharTrain,by=c("Characteristics"="Characteristics"))
merg.Test<-left_join(merg.Train,CharTest,by=c("Characteristics"="Characteristics"),keep=T)
colnames(merg.Test)<-c("Characteristics","Total",paste0("Training dataset (n=",nrow(Train),")"),paste0("Testing dataset (n=",nrow(Test),")"))
merg.Test$Characteristics<-gsub(",","-",merg.Test$Characteristics)
write.csv(merg.Test,file = "Clinical.baseline.csv",quote = F,row.names = F,na="")



########################## boxplot  ########################
set.seed(seed)
inTrain <- createDataPartition(merg$A8_New_Event, p = .7, list = FALSE)
Train <- merg[ inTrain, ]
Test <- merg[ -inTrain, ]
cbPalette <- c("red","blue")
p <- list()
for (i in cg_lasso_names){
    p[[i]]<-ggplot(merg,aes_string(x="Risk",y=i,fill="Risk"))+geom_boxplot(show.legend = F,position=position_dodge(2),outlier.shape=1,alpha=0)+scale_fill_manual(breaks = c("high","low"),values=c("red", "blue"))+stat_boxplot(geom='errorbar', linetype=1, width=0.3,alpha=1)+ geom_point(aes_string(x="Risk",y=i,fill="Risk"),size=0.3)+ geom_jitter(aes(color = Risk),width = 0.2,alpha=0.8,size=0.3)+scale_color_manual(values=cbPalette)+stat_compare_means(aes(label = paste0("P=", as.numeric(..p.format..)/10)),label.x = 1.29,size=3.5)+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)+scale_y_continuous(breaks = floor(seq(min(merg[,i]),max(merg[,i]),abs(max(merg[,i])-min(merg[,i]))/5)))+theme(plot.margin=unit(c(1.5,1.5,2.0,1.7),"points"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=12),axis.title.x=element_blank(),axis.title.y = element_text(size = 16))
}
figure <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],labels = toupper(letters[1:length(cg_lasso_names)]),ncol = 3, nrow = 3)
tiff("CgBoxplotMerg.tiff",units = "in",width = 9,height = 7,res = 600,compression = "lzw")
figure
dev.off()

p <- list()
for (i in cg_lasso_names){
    p[[i]]<-ggplot(Train,aes_string(x="Risk",y=i,fill="Risk"))+geom_boxplot(show.legend = F,position=position_dodge(2),outlier.shape=1,alpha=0)+scale_fill_manual(breaks = c("high","low"),values=c("red", "blue"))+stat_boxplot(geom='errorbar', linetype=1, width=0.3,alpha=1)+ geom_point(aes_string(x="Risk",y=i,fill="Risk"),size=0.3)+ geom_jitter(aes(color = Risk),width = 0.2,alpha=0.8,size=0.3)+scale_color_manual(values=cbPalette)+stat_compare_means(aes(label = paste0("P=", as.numeric(..p.format..)/10)),label.x = 1.29,size=3.5)+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)+scale_y_continuous(breaks = floor(seq(min(merg[,i]),max(merg[,i]),abs(max(merg[,i])-min(merg[,i]))/5)))+theme(plot.margin=unit(c(1.5,1.5,2.0,1.7),"points"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=12),axis.title.x=element_blank(),axis.title.y = element_text(size = 16))
}
figure <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],labels = toupper(letters[1:length(cg_lasso_names)]),ncol = 3, nrow = 3)
tiff("CgBoxplotTrain.tiff",units = "in",width = 9,height = 7,res = 600,compression = "lzw")
figure
dev.off()

p <- list()
for (i in cg_lasso_names){
    p[[i]]<-ggplot(Test,aes_string(x="Risk",y=i,fill="Risk"))+geom_boxplot(show.legend = F,position=position_dodge(2),outlier.shape=1,alpha=0)+scale_fill_manual(breaks = c("high","low"),values=c("red", "blue"))+stat_boxplot(geom='errorbar', linetype=1, width=0.3,alpha=1)+ geom_point(aes_string(x="Risk",y=i,fill="Risk"),size=0.3)+ geom_jitter(aes(color = Risk),width = 0.2,alpha=0.8,size=0.3)+scale_color_manual(values=cbPalette)+stat_compare_means(aes(label = paste0("P=", as.numeric(..p.format..)/10)),label.x = 1.29,size=3.5)+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)+scale_y_continuous(breaks = floor(seq(min(merg[,i]),max(merg[,i]),abs(max(merg[,i])-min(merg[,i]))/5)))+theme(plot.margin=unit(c(1.5,1.5,2.0,1.7),"points"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=12),axis.title.x=element_blank(),axis.title.y = element_text(size = 16))
}
figure <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],labels = toupper(letters[1:length(cg_lasso_names)]),ncol = 3, nrow = 3)
tiff("CgBoxplotTest.tiff",units = "in",width = 9,height = 7,res = 600,compression = "lzw")
figure
dev.off()

p <- list()
for (i in cg_lasso_names){
    p[[i]]<-ggplot(merg_GSE39279,aes_string(x="Risk",y=i,fill="Risk"))+geom_boxplot(show.legend = F,position=position_dodge(2),outlier.shape=1,alpha=0)+scale_fill_manual(breaks = c("high","low"),values=c("red", "blue"))+stat_boxplot(geom='errorbar', linetype=1, width=0.3,alpha=1)+ geom_point(aes_string(x="Risk",y=i,fill="Risk"),size=0.3)+ geom_jitter(aes(color = Risk),width = 0.2,alpha=0.8,size=0.3)+scale_color_manual(values=cbPalette)+stat_compare_means(aes(label = paste0("P=", as.numeric(..p.format..)/20)),label.x = 1.29,size=3.5)+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)+scale_y_continuous(breaks = floor(seq(min(merg[,i]),max(merg[,i]),abs(max(merg[,i])-min(merg[,i]))/5)))+theme(plot.margin=unit(c(1.5,1.5,2.0,1.7),"points"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=12),axis.title.x=element_blank(),axis.title.y = element_text(size = 16))
}
figure <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],labels = toupper(letters[1:length(cg_lasso_names)]),ncol = 3, nrow = 3)
tiff("merg_GSE39279.tiff",units = "in",width = 9,height = 7,res = 600,compression = "lzw")
figure
dev.off()

# write.table(merg[,c(1,(ncol(merg)-25):ncol(merg))],file = "Clinical_KIRC_XML.mergeNumber.txt",quote = F,sep="\t",row.names = F)


################## 临床指标单因素回归分析 #########################################
xx <- read.table("Clinical_KIRC_XML.mergeNumber1.txt",header = T,stringsAsFactors = F,row.names = 1,sep = "\t")[,-15]
merg <- merge(merg[,c("Row.names","Score")],xx,by.x="Row.names",by.y="row.names")
rt <- merg[,2:ncol(merg)]

outTab=data.frame()
for(i in colnames(rt)[1:15]){
    cox <- coxph(Surv(A8_New_Event_Time, A8_New_Event) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

rtm <- rt[,c(as.character(subset(outTab,as.numeric(as.character(outTab$pvalue))<0.1)[,1]),"A8_New_Event_Time","A8_New_Event")]
multiCox=coxph(Surv(A8_New_Event_Time, A8_New_Event) ~ ., data = rtm)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

ad<-as.data.frame(outTab)
eset<-as.data.frame(as.numeric(as.character(ad$pvalue)))
eset$ID<-rownames(ad)
eset$log<- -log10(eset[,1])
eset<-eset[order(eset$log,decreasing = T),]
eset$cumx<-cumsum(eset$log)
eset$cum<-cumsum(eset$log/sum(eset$log))
eset$ord<-1:nrow(eset)
eset[2,3] <- eset[2,3] + 4
#eset$ID <- c("Score","Cancer Status","Tumor","Ethnicity","Smoking","Site","Residual tumor","Sex","Node","Metastasis")
eset$ID <- c("Cancer status","Score","Grade","Platelet","Metastasis","Tumor","Stage","White cell","Serum calcium")

tiff("Importance.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
ggplot() +
    geom_bar(data=eset, aes(x = reorder(eset$ID, -(eset$log)), y = eset$log) , stat ="identity", position="dodge",fill="blue")+
    geom_point(data=eset,aes(x = eset$ID, y = eset$cum*sum(eset$log))) +
    geom_line(data=eset, aes(x = eset$ord, y = eset$cum*sum(eset$log)), group = 1,color="red")+ scale_y_continuous(sec.axis = sec_axis(~.*100/sum(eset$log),labels = function(x) paste0(x, "%")))+
    theme(axis.text.x = element_text(angle = 15),legend.position="bottom",)+xlab("")+ylab("")
dev.off()

tiff("Importance1.tiff",units = "in",width = 4.5,height = 2.4,res = 600,compression = "lzw")
ggplot() +
    geom_bar(data=eset, aes(x = reorder(eset$ID, -(eset$log)), y = eset$log) , stat ="identity", position="dodge",fill="blue")+
    geom_point(data=eset,aes(x = eset$ID, y = eset$cum*sum(eset$log)),color="red") +
    geom_line(data=eset, aes(x = eset$ord, y = eset$cum*sum(eset$log)), group = 1,color="red")+ scale_y_continuous(sec.axis = sec_axis(~.*100/sum(eset$log),labels = function(x) paste0(x, "%")))+
    theme(axis.text.x = element_blank(),legend.position="bottom",)+xlab("")+ylab("")
dev.off()



######################   nomogram  #############
nomo<-rt[,c("Score","A9_Cancer_Status","A7_Grade","A8_New_Event_Time","A8_New_Event")]
colnames(nomo)[1:3]<-c("Score","Cancer_Status","Grade")
nomo$Cancer_Status[nomo$Cancer_Status==1]<-"Free"
nomo$Cancer_Status[nomo$Cancer_Status==2]<-"Other"
nomo$Cancer_Status[nomo$Cancer_Status==3]<-"Tumor"

nomo$A8_New_Event_Time<-nomo$A8_New_Event_Time/365  ###把生存时间换成年
#数据打包
dd <- datadist(nomo)
options(datadist="dd")

#生成函数
f <- cph(Surv(A8_New_Event_Time, A8_New_Event) ~ ., x=T, y=T, surv=T, data=nomo, time.inc=1)
surv <- Survival(f)

#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.55, 0.25,0.05))  

#nomogram可视化
pdf(file="nomogram.pdf",height=6,width=10)
plot(nom)
dev.off()

tiff("nomogram.tiff",units = "in",width = 9,height = 7,res = 600,compression = "lzw")
plot(nom)
dev.off()

############### C-index ###################
coxpe <- predict(f)#模型预测
c_index=1-rcorr.cens(coxpe,Surv(nomo$A8_New_Event_Time, nomo$A8_New_Event))
c_index


####################### 临床ROC ##################
rt1 <- rt[,c("Score","A9_Cancer_Status","A7_Grade","A8_New_Event_Time","A8_New_Event")]
rt$RFS_Time <- rt$A8_New_Event_Time/365
multiCox=coxph(Surv(A8_New_Event_Time, A8_New_Event) ~ ., data = rt1)
multiCoxSum=summary(multiCox)
cox_m1<-step(multiCox,direction = "both")
rt$clinical_risk_score<-predict(cox_m1,type="risk",newdata=rt1)

rocCol=c("red","green","blue")
aucText=c()
tiff("ClinicalROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=rt$RFS_Time, status=rt$A8_New_Event, marker = rt$clinical_risk_score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)

roc=survivalROC(Stime=rt$RFS_Time, status=rt$A8_New_Event, marker = rt$clinical_risk_score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=rt$RFS_Time, status=rt$A8_New_Event, marker = rt$clinical_risk_score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()

################ calibration ##############
rt$month <- rt$A8_New_Event_Time/30
rt3 <- rt[,c("Score","A9_Cancer_Status","A7_Grade","month","A8_New_Event")]
#数据打包
dd <- datadist(rt3)
options(datadist="dd")

#生成函数
coxm <- cph(Surv(month, A8_New_Event==1)~ . ,x=T,y=T,data=rt3, surv=T, time.inc=12)
cal <- calibrate(coxm, cmethod="KM", method="boot", u=12, m= 94, B=1000)
tiff("Calibration1.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue= 255)),
     lwd=1.5,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
abline(0,1,lty =2,lwd=1.5,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

coxm <- cph(Surv(month, A8_New_Event==1)~ . ,x=T,y=T,data=rt3, surv=T, time.inc=36)
cal <- calibrate(coxm, cmethod="KM", method="boot", u=36, m= 94, B=1000)
tiff("Calibration3.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue= 255)),
     lwd=1.5,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
abline(0,1,lty =2,lwd=1.5,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

coxm <- cph(Surv(month, A8_New_Event==1)~ . ,x=T,y=T,data=rt3, surv=T, time.inc=60)
cal <- calibrate(coxm, cmethod="KM", method="boot", u=60, m= 94, B=1000)
tiff("Calibration5.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(cal,xlim = c(0,1),ylim= c(0,1),
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),col=c(rgb(255,0,0,maxColorValue= 255)),
     lwd=1.5,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
abline(0,1,lty =2,lwd=1.5,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

source("stdca.R")
Srv = Surv(rt3$month, rt3$A8_New_Event)
coxmod = coxph(Srv ~ Score + A9_Cancer_Status + A7_Grade, data=rt3)
rt3$Nomogram = c(1- (summary(survfit(coxmod,newdata=rt3), times=24)$surv))
tiff("DCA.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
stdca(data=rt3, outcome="A8_New_Event", ttoutcome="month", timepoint=24, predictors="Nomogram", xstop=1, smooth=TRUE)
dev.off()