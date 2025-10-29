

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

clusterFile="cluster.txt"           #分型的结果文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
setwd("C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/17.clusterImm")      #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", row.names(immune))
data=immune[group=="Treat",,drop=F]

#读取分型的结果文件
Cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(Cluster))
rt=cbind(data[sameSample,,drop=F], Cluster[sameSample,"Cluster",drop=F])
rt=rt[order(rt$Cluster, decreasing=F),]
conNum=nrow(rt[rt$Cluster=="C1",])
treatNum=nrow(rt[rt$Cluster=="C2",])

##########绘制柱状图##########
data=t(rt[,-ncol(rt)])
pdf(file="barplot.pdf", width=14.5, height=8)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data, col=col, xaxt="n", yaxt="n", ylab="Relative Percent", cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"C1",cex=2)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5 , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"C2",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=rt
data=melt(data, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Cluster))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Cluster",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Cluster",
				  add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()



