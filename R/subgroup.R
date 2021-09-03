###################### Stage #############################
stagesub1<-subset(merg,merg$A6_Stage %in% c("Stage I","Stage II"))
stagesub2<-subset(merg,merg$A6_Stage %in% c("Stage III","Stage IV"))

fit1 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = stagesub1)
fit2 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = stagesub2)


sdf <- survdiff(Lusurv~Risk,data = stagesub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("Stage12Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit1,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomright",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.2,0.4,pp,cex=1.5)
dev.off()

sdf <- survdiff(Lusurv~Risk,data = stagesub2)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("Stage34Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit2,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("center",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.6,0.3,pp,cex=1.5)
dev.off()


rocCol=c("red","green","blue")
aucText=c()
tiff("Stage12ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=stagesub1$OS_Time, status=stagesub1$A8_New_Event, marker = stagesub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=stagesub1$OS_Time, status=stagesub1$A8_New_Event, marker = stagesub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=stagesub1$OS_Time, status=stagesub1$A8_New_Event, marker = stagesub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和stagesub1即可



rocCol=c("red","green","blue")
aucText=c()
tiff("Stage34ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=stagesub2$OS_Time, status=stagesub2$A8_New_Event, marker = stagesub2$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=stagesub2$OS_Time, status=stagesub2$A8_New_Event, marker = stagesub2$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=stagesub2$OS_Time, status=stagesub2$A8_New_Event, marker = stagesub2$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  



################################  gender subgroup ############################
gendersub1<-subset(merg,merg$A18_Sex=="FEMALE")
gendersub2<-subset(merg,merg$A18_Sex=="MALE")
fit1 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = gendersub1)
fit2 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = gendersub2)

sdf <- survdiff(Lusurv~Risk,data = gendersub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("GenderFemaleSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit1,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomright",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.2,0.4,pp,cex=1.5)
dev.off()

sdf <- survdiff(Lusurv~Risk,data = gendersub2)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("GenderMaleSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit2,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "RFS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.6,0.9,pp,cex=1.5)
dev.off()

rocCol=c("red","green","blue")
aucText=c()
tiff("GenderFemaleROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=gendersub1$OS_Time, status=gendersub1$A8_New_Event, marker = gendersub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=gendersub1$OS_Time, status=gendersub1$A8_New_Event, marker = gendersub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=gendersub1$OS_Time, status=gendersub1$A8_New_Event, marker = gendersub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和stagesub1即可

rocCol=c("red","green","blue")
aucText=c()
tiff("GenderMaleROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=gendersub2$OS_Time, status=gendersub2$A8_New_Event, marker = gendersub2$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=gendersub2$OS_Time, status=gendersub2$A8_New_Event, marker = gendersub2$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=gendersub2$OS_Time, status=gendersub2$A8_New_Event, marker = gendersub2$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和gendersub2即可



########################## Age #####################
agesub1<-subset(merg,merg$age_at_initial_pathologic_diagnosis>60)
agesub2<-subset(merg,merg$age_at_initial_pathologic_diagnosis<=60)
fit1 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = agesub1)
fit2 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = agesub2)

sdf <- survdiff(Lusurv~Risk,data = agesub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("AgeOver60Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit1,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.7,0.2,pp,cex=1.5)
dev.off()

sdf <- survdiff(Lusurv~Risk,data = agesub2)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("AgeLess60Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit2,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(-0.01,-0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.5,0.7,pp,cex=1.5)
dev.off()

rocCol=c("red","green","blue")
aucText=c()
tiff("AgeOver60ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=agesub1$OS_Time, status=agesub1$A8_New_Event, marker = agesub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=agesub1$OS_Time, status=agesub1$A8_New_Event, marker = agesub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=agesub1$OS_Time, status=agesub1$A8_New_Event, marker = agesub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和agesub1即可

rocCol=c("red","green","blue")
aucText=c()
tiff("AgeLess60ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=agesub2$OS_Time, status=agesub2$A8_New_Event, marker = agesub2$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=agesub2$OS_Time, status=agesub2$A8_New_Event, marker = agesub2$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=agesub2$OS_Time, status=agesub2$A8_New_Event, marker = agesub2$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()   #修改tiff文件名和gendersub2即可




##########################  Sites #################################
leftsub1<-subset(merg,merg$laterality == "Left")
rightsub1<-subset(merg,merg$laterality == "Right")
fit1 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = leftsub1)
fit2 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = rightsub1)

sdf <- survdiff(Lusurv~Risk,data = leftsub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("LeftSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit1,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.7,0.2,pp,cex=1.5)
dev.off()

sdf <- survdiff(Lusurv~Risk,data = rightsub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("RightSurvival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit2,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomright",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(-0.01,-0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.5,0.75,pp,cex=1.5)
dev.off()

rocCol=c("red","green","blue")
aucText=c()
tiff("LeftROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=leftsub1$OS_Time, status=leftsub1$A8_New_Event, marker = leftsub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=leftsub1$OS_Time, status=leftsub1$A8_New_Event, marker = leftsub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=leftsub1$OS_Time, status=leftsub1$A8_New_Event, marker = leftsub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和agesub1即可

rocCol=c("red","green","blue")
aucText=c()
tiff("RightROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=rightsub1$OS_Time, status=rightsub1$A8_New_Event, marker = rightsub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=rightsub1$OS_Time, status=rightsub1$A8_New_Event, marker = rightsub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=rightsub1$OS_Time, status=rightsub1$A8_New_Event, marker = rightsub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()   #修改tiff文件名和gendersub2即可





##########################  Grade #################################
G12sub1<-subset(merg,merg$A7_Grade %in% c("G1","G2"))
G34sub1<-subset(merg,merg$A7_Grade %in% c("G3","G4"))
fit1 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = G12sub1)
fit2 <- survfit(Surv(A8_New_Event_Time, A8_New_Event) ~ Risk, data = G34sub1)

sdf <- survdiff(Lusurv~Risk,data = G12sub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("G12Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit1,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(0.01,0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.7,0.2,pp,cex=1.5)
dev.off()

sdf <- survdiff(Lusurv~Risk,data = G34sub1)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
tiff("G34Survival.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
plot(fit2,conf.int="none",col = c("red","blue"),lwd = 5,mark.time = T,xlab = "Days",ylab = "OS",las=1,cex.axis=1.5,cex.lab=1.45,yaxt="n",bty='l',ylim=c(0, 1), yaxs="i",xlim=c(0, max(merg$A8_New_Event_Time)), xaxs="i")   #bty去掉边框
pp <- ifelse(p.val < 0.001,paste0("P = ",scientific_format(1)(p.val)),paste0("P = ",round(p.val,3)))
legend("bottomleft",c("High","Low"),col = c("red","blue"),lwd = 5,text.font = 1,cex=1.3,box.col = "white",bg = "white", inset=c(-0.01,-0.015))
axis(2,at=c(1,0.8,0.6,0.4,0.2,0),c(100,80,60,40,20,0),las=1,cex.axis=1.5)
axis(1,lwd=3,labels = F)
axis(2,lwd=3,labels = F)
text(max(merg$A8_New_Event_Time)*0.5,0.5,pp,cex=1.5)
dev.off()

rocCol=c("red","green","blue")
aucText=c()
tiff("G12ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=G12sub1$OS_Time, status=G12sub1$A8_New_Event, marker = G12sub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=G12sub1$OS_Time, status=G12sub1$A8_New_Event, marker = G12sub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=G12sub1$OS_Time, status=G12sub1$A8_New_Event, marker = G12sub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()  #修改tiff文件名和agesub1即可

rocCol=c("red","green","blue")
aucText=c()
tiff("G34ROC.tiff",units = "in",width = 4.5,height = 3.8,res = 600,compression = "lzw")
roc=survivalROC(Stime=G34sub1$OS_Time, status=G34sub1$A8_New_Event, marker = G34sub1$Score, predict.time = 5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd=3,cex.lab=1, cex.axis=1.5, cex.main=0.6, cex.sub=0.6,  pch=19, col.lab="black",cex.main=1.5,family="serif")
box(lwd=3)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(c(0,1),c(0,1),col = "gray", lty = 8 ,lwd=2)
roc=survivalROC(Stime=G34sub1$OS_Time, status=G34sub1$A8_New_Event, marker = G34sub1$Score, predict.time = 3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=G34sub1$OS_Time, status=G34sub1$A8_New_Event, marker = G34sub1$Score, predict.time = 1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex=0.8)
dev.off()   #修改tiff文件名和gendersub2即可
