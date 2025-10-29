

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)      #引用包
expFile="diffGeneExp.txt"          #表达数据文件
workDir="C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/14.cluster"     #工作目录
setwd(workDir)      #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#去除对照组样品, 只保留实验组样品
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]

#聚类
maxK=9     #设置最大的k值
results=ConsensusClusterPlus(data,
              maxK=maxK,  #最大的K值，形成一系列梯度
              reps=50,   #重复抽样的数目
              pItem=0.8,  #选择80%的样本进行重复抽样
              pFeature=1,  #选择100%的基因进行重复抽样
              title=workDir,  #输出结果的文件夹名字，包含了输出的图片
              clusterAlg="km",   #层次聚类的算法
              distance="euclidean",  #距离矩阵的算法
              seed=123456,  #seed, 随机种子，用于重复结果
              plot="png")

#一致性打分
calcICL(results, title="consensusScore", plot="png")

#输出分型结果
clusterNum=2        #分几类，根据前面的图形判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)

#教程：https://cloud.tencent.com/developer/article/2019378

#consensus 累计分布图 CDF (consensus010.png)：对于每个K对应的consensus matrix,  采用100个bin的柱状图来计算累计分布，CDF图可以用来帮助决定最佳的K值

##选取最下的

#delta area plot (consensus011.png)：对于每个K, 计算K和K-1相比，CDF 曲线下面积的相对变化，对于K=2, 因为没有K=1, 所以是totla CDF curve area

##选取增加不明显的点作为最佳的K值,类似单细胞的pc.number



