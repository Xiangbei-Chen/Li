

install.packages("org.Hs.eg.db")
install.packages("stringi")
install.packages("ggplot2")
install.packages("digest")
install.packages("GOplot")

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")

pvalueFilter=0.05        #p值过滤条件
p.adjustFilter=1      #矫正后的p值过滤条件

#定义颜色
colorSel="p.adjust"
if(p.adjustFilter>0.05){
	colorSel="pvalue"
}

setwd("C:\\Users\\DELL\\Desktop\\21.GO")       #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取差异基因的名称, 将基因名称转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<p.adjustFilter),]
#保存显著富集的结果
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel, drop=TRUE) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#气泡图
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio", label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#获取GO的信息
go=data.frame(Category=GO$ONTOLOGY, ID=GO$ID, Term=GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)
#读取基因的logFC
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图参数
circ <- circle_dat(go, genelist)
termNum =8       #设置展示GO数目
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum=200      #限定基因数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)

#绘制GO的圈图
pdf(file="GOcircos.pdf", width=13, height=6)
GOCircle(circ, table.legend=T, label.size=5, nsub=termNum)
dev.off()

#绘制和弦图
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOChord.pdf", width=11, height=11.5)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名称与圆圈之间的距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 6)       #GO字体大小
dev.off()

#GO聚类图
pdf(file="GOcluster.pdf", width=14, height=11)
GOCluster(circ, 
          go$Term[1:termNum], 
          lfc.space = 0.2,        #logFC与树之间的空隙大小
          lfc.width = 1,          #logFC的圆圈宽度
          term.space = 0.2,       #logFC与GO之间空隙的大小
          term.width = 1)         #GO圆圈的宽度
dev.off()



