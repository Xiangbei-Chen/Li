

#install.packages("RCircos")


library("RCircos")       #引用包
setwd("C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/09.Rcircos")      #设置工作目录

#初始化圈图
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#设置圈图的参数
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#输出文件
pdf(file="RCircos.pdf", width=10, height=10)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#读取基因注释文件，标注基因的名称
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


 

