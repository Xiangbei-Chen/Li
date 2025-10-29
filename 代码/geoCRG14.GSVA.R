

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("SpatialExperiment")

BiocManager::install("GSVA", version = "3.17")  # 对应 Bioconductor 3.17 / R 4.4

#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(devtools)
library(IRanges)

expFile="normalize.txt"              #表达数据文件
clusterFile="cluster.txt"            #分型的结果文件
gmtFile="c2.cp.kegg.symbols.gmt"     #基因集文件
setwd("C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/18.GSVA")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA分析
ssgseaScore=gsva(data, geneSets, method='gsva')
#对GSVA的打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))

#通路差异分析
outTab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	pvalue=test$p.value
	t=test$statistic
	if(pvalue<0.05){
		Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
		outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
	}
}

#绘制柱状图
termNum=10      #展示通路的数目
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot_KEGG.pdf", width=14, height=14)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("blue3", "red3"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, legend="right", title="",
		xlab="Term", ylab="t value of GSVA score, C2 vs C1",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()


###https://zhuanlan.zhihu.com/p/506912398?utm_id=0

###以上为GSVA和GSEA的富集分析的比较，写的很详细。请注意



