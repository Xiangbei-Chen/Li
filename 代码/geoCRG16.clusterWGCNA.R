#引用包
library(limma)
library(gplots)
library(WGCNA)

expFile="normalize.txt"       #表达数据文件
clusterFile="cluster.txt"     #分型的结果文件
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//20.clusterWGCNA_Variable")     #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.25)))
data=data[selectGenes,]

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=data[,nameC1,drop=F]
dataC2=data[,nameC2,drop=F]
conCount=ncol(dataC1)
treatCount=ncol(dataC2)
data=cbind(dataC1, dataC2)
datExpr0=t(data)

###检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###样品聚类
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "cluster01.sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 20000, col = "red")
dev.off()

###删除剪切线以下的样品
clust=cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples=(clust==1)
datExpr0=datExpr0[keepSamples, ]



# ===== 抽取临床性状 =====
traitData=data.frame(C1=c(rep(1,conCount),rep(0,treatCount)),
                     C2=c(rep(0,conCount),rep(1,treatCount)))
rownames(traitData)=colnames(data)
sameSample=intersect(rownames(datExpr0), rownames(traitData))
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]




###样品聚类,得到样品聚类的热图
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="cluster02.sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="cluster03.scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9

###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# ===== 构建 TOM =====
softPower <- 17  # ← 根据 sft 分析结果选择更优的 power
adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM


# 层次聚类
geneTree = hclust(as.dist(dissTOM), method = "average")

# 输出聚类树
pdf("04.gene_clustering.pdf", width = 12, height = 9)
plot(geneTree,
     main = "Gene clustering on TOM-based dissimilarity",
     xlab = "", sub = "", labels = FALSE)
dev.off()


###动态剪切模块识别
minModuleSize = 100      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="cluster05.Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###对模块进行聚类,找出模块之间的相似性
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="cluster06.Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

#绘制模块基因的热图
moduleColors=dynamicColors
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
select = sample(nGenes, size=1000)
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method="average")
selectColors = moduleColors[select]
#sizeGrWindow(9,9)
plotDiss=selectTOM^softPower
diag(plotDiss)=NA
myheatcol = colorpanel(250, "red", "orange", "lemonchiffon")    #设置热图颜色(白色背景）
pdf(file="cluster07.TOMplot.pdf", width=7, height=7)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()


###模块与性状数据热图
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="cluster08.Module_trait.pdf", width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3.5, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


###计算MM和GS值
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###批量输出性状和模块散点图
trait="C2"
traitColumn=match(trait,traitNames)  
for (module in modNames){
	column = match(module, modNames)
	moduleGenes = moduleColors==module
	if (nrow(geneModuleMembership[moduleGenes,]) > 1){
	    outPdf=paste("cluster09.", trait, "_", module,".pdf",sep="")
	    pdf(file=outPdf,width=7,height=7)
	    par(mfrow = c(1,1))
	    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
	                         xlab = paste("Module Membership in", module, "module"),
	                         ylab = paste("Gene significance for ",trait),
	                         main = paste("Module membership vs. gene significance\n"),
	                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
	    abline(v=0.8,h=0.5,col="red")
	    dev.off()
	}
}


###输出GS_MM数据
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "cluster.GS_MM.xls",sep="\t",row.names=F)


###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("cluster.module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###输出每个模块的核心基因
geneSigFilter=0.5         #基因重要性的过滤条件
moduleSigFilter=0.8      #基因与模块相关性的过滤条件
cli="GS.C2"
for(mol in unique(geneInfo$moduleColor)){
	geneInfoMol=geneInfo[geneInfo$moduleColor==mol,]
	mmi=paste0("MM", mol)
	geneInfoMol2=geneInfoMol[((abs(geneInfoMol[,mmi])>moduleSigFilter) & (abs(geneInfoMol[,cli])>geneSigFilter)),]
	write.table(geneInfoMol2[,1], file =paste0("cluster.hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}




