
#install.packages("VennDiagram")


library(VennDiagram)       #引用包
diseaseFile="hubGenes_MMturquoise.txt"                 #疾病共表达分析的结果文件
clusterFile="cluster.hubGenes_MMturquoise.txt"         #分型共表达分析的结果文件
setwd("F://WHW//Test_AZ_2022.12.31.6th//21.venn")     #设置工作目录
geneList=list()

#读取疾病共表达分析的结果文件
rt=read.table(diseaseFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["Disease WGCNA"]]=uniqGene     #将疾病WGCNA的核心基因保存到geneList里面

#读取分型共表达分析的结果文件
rt=read.table(clusterFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["Cluster WGCNA"]]=uniqGene     #把分型WGCNA的核心基因保存到geneList里面

#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex=1)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集核心基因的列表
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)

