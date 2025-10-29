

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)              #引用包
expFile="normalize.txt"     #表达数据文件
geneFile="gene.txt"         #基因列表文件
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//06.crgExp_Variable")     #设置工作目录

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因列表文件，获取铜死亡相关基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出铜死亡基因的表达数据
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="crgGeneExp.txt", sep="\t", quote=F, col.names=F)



