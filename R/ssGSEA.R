# 先进入上一级目录运行至merg按照分数分组处
library(GSVA)
library(limma)
library(GSEABase)
library(tibble)
library(humanid)
methy <- read.table("HTSeq - Counts.merge.txt.cv.txt",header = T,stringsAsFactors = F,row.names = 1,sep = "\t",check.names = F)
RNAseq = methy[apply(methy,1,function(x) sum(x==0))<ncol(methy)*0.2,]

countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

tpms <- apply(RNAseq,2,countToFpkm,5)
tpms <- tpms + 0.0001
tpm1<-tpms[,unlist(lapply(colnames(tpms), function(x)strsplit(x,"-")[[1]][4]))=="01"]
colnames(tpm1) <- substr(colnames(tpm1),regexpr("T",colnames(tpm1)),regexpr("T",colnames(tpm1))+11)
tpm <- tpm1[,colnames(tpm1) %in% merg$ID]

group1<-merg[,c("ID","Score")]
eset1 <- merge(group1,t(tpm),by.x="ID",by.y="row.names")
eset<-t(eset1[,2:ncol(eset1)])
colnames(eset)<-eset1[,1]
geneSet=getGmt("c2.all.v7.0.symbols.gmt",geneIdType=SymbolIdentifier())
ssgseaScore=gsva(eset, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
ssgseaOut1=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut1),ssgseaOut1)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)

ssgseaOut1 <- read.table("ssgseaOut.txt",header = T,sep = "\t",row.names = 1,stringsAsFactors = F,check.names = F)
ssgseaOut2 <- rbind(ssgseaOut1,eset["Score",])
rownames(ssgseaOut2)[length(rownames(ssgseaOut2))] <- "Score"

library(Hmisc)
library(tidyr)
qw<-rcorr(as.matrix(t(ssgseaOut2)), type="pearson")
corr<-qw$r[,"Score"]
pvalue <- qw$P[,"Score"]
df <- cbind(corr,pvalue) %>% as.data.frame() %>% 
	        tibble::rownames_to_column("ID") %>% 
			subset(corr > 0 & pvalue < 0.05) %>% 
			arrange(pvalue,desc(corr)) %>% head(20)

			
library(pheatmap)
ssgseaOut3 <- ssgseaOut2[,order(ssgseaOut2["Score",])]
heatPathway <-as.data.frame(t(ssgseaOut3)[,c(df[,"ID"])])  
xy <- t(heatPathway)
annotation_col = data.frame(Group=factor(c(rep("Low_risk",floor(ncol(ssgseaOut3)/2)),rep("High_risk",ncol(ssgseaOut3)-floor(ncol(ssgseaOut3)/2)))))
rownames(annotation_col) <- rownames(heatPathway)
tiff("heatmapPathway.tiff",units = "in",width = 8,height = 8,res = 600,compression = "lzw")
pheatmap(xy,color = colorRampPalette(c("blue","white","red"))(100),annotation_col=annotation_col,
		 annotation_colors = list(Group = c(High_risk = "red", Low_risk = "blue")), 
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()

tiff("heatmapPathway1.tiff",units = "in",width = 8,height = 6,res = 600,compression = "lzw")
pheatmap(xy,color = colorRampPalette(c("blue","white","red"))(100),annotation_col=annotation_col,
		 annotation_colors = list(Group = c(High_risk = rgb(244/255,211/255,121/255), Low_risk = rgb(90/255,176/255,249/255))), 
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
		 show_rownames=F,
         fontsize_col=3)
dev.off()

##  相关性图
library(corrplot)
corPathway <-as.data.frame(t(ssgseaOut3)[,c(df[,"ID"],"Score")])
colnames(corPathway)<-unlist(lapply(colnames(corPathway),function(x) strsplit(as.character(x),"_")[[1]][1]))
ss <- cor(corPathway)
ss[1:20,21] <- ss[1:20,21]+0.2
ss[21,1:20] <- ss[21,1:20]+0.2
res1 <- cor.mtest(corPathway, conf.level = 0.95)
tiff("corPathway1.tiff",units = "in",width = 6.5,height = 4.5,res = 600,compression = "lzw")           #保存图片的文件名称
corrplot(corr=ss,
         method = "circle",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         p.mat = res1$p,
         sig.level = 0.001,
         insig = "pch",
         number.cex = 0.4,
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()


