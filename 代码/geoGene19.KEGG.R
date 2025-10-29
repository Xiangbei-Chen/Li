
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("digest")
#install.packages("GOplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")

pvalueFilter=0.05       #p值过滤条件
p.adjustFilter=1     #矫正后的p值过滤条件

#定义颜色
colorSel="p.adjust"
if(p.adjustFilter>0.05){
	colorSel="pvalue"
}
	
setwd("C:\\Users\\DELL\\Desktop\\22.KEGG")      #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)      #读取输入文件

#提取差异基因的名称, 将基因名称转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=9, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#获取通路信息
kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
#读取基因的差异信息
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图的参数
circ <- circle_dat(kegg, genelist)
termNum =8       #显示通路的数目
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=200      #显示基因的数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)

#绘制通路的圈图
pdf(file="KEGGcircos.pdf", width=11, height=6)
GOCircle(circ, table.legend=T, label.size=5, nsub=termNum)
dev.off()

#绘制和弦图
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGGchord.pdf", width=12, height=12.6)
GOChord(chord, 
        space = 0.001,           #基因之间的距离
        gene.order = 'logFC',    #基因的排序方式(按照logFC值对基因进行排序)
        gene.space = 0.25,       #基因名称跟圆圈的相对距离
        gene.size = 5,           #基因名称字体大小 
        border.size = 0.1,       #线条的粗细
        process.label = 6)       #通路字体的大小
dev.off()

#通路的聚类图
pdf(file="KEGGcluster.pdf", width=15, height=10)
GOCluster(circ, 
          kegg$Term[1:termNum], 
          lfc.space = 0.2,        #logFC与树之间的空隙大小
          lfc.width = 1,          #logFC的圆圈宽度
          term.space = 0.2,       #logFC与通路间的空隙大小
          term.width = 1)         #通路圆圈的宽度
dev.off()          




