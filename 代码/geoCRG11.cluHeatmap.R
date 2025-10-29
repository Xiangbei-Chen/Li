

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="cluster.txt"      #分型的结果文件
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//15.clusterHeatmap")      #设置工作目录

#读取输入文件
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$Cluster),]

#准备热图数据
data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

#定义热图注释的颜色
bioCol=c("#DFB6BC","#8EB5E0","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
crgCluCol=bioCol[1:length(levels(factor(Type$Cluster)))]
names(crgCluCol)=levels(factor(Type$Cluster))
ann_colors[["Cluster"]]=crgCluCol

#绘制热图
pdf("heatmap.pdf", width=6, height=3.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#407BD0",2), "white", rep("#A32A31",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

#把数据转换成ggplot2输入文件
data=melt(rt, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Cluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Cluster",
	     palette = crgCluCol,
	     width=0.8,
	     add="point")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=10, height=5)
print(p1)
dev.off()


