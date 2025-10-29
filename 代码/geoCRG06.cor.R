# 安装包（如未安装）
# install.packages("corrplot")
# install.packages("circlize")

# 加载库
library(corrplot)
library(circlize)

# 设置工作目录和输入文件路径
inputFile <- "diffGeneExp.txt"
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//10.cor")

# 读取数据
data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 只保留Treat组
group <- gsub("(.*)\\_(.*)", "\\2", colnames(data))
data <- data[, group == "Treat", drop = FALSE]
rt <- t(data)

# 计算相关性矩阵
cor1 <- cor(rt)
cor1[cor1 == 1] <- 0  # 避免自相关显示为1

# 设置统一配色（与 corrplot 相同）
col_vector <- colorRampPalette(c("#8BACD1", "white", "#C17F9E"))(100)

# 将相关性值 [-1, 1] 映射到颜色索引 [1, 100]
scaled_cor <- round((cor1 + 1) / 2 * 99) + 1
scaled_cor[cor1 == 0] <- NA  # 对于0的值不分配颜色

# 生成颜色矩阵
col1 <- matrix(NA, nrow = nrow(cor1), ncol = ncol(cor1))
for (i in 1:nrow(cor1)) {
  for (j in 1:ncol(cor1)) {
    if (!is.na(scaled_cor[i, j])) {
      col1[i, j] <- col_vector[scaled_cor[i, j]]
    }
  }
}

# 绘制圈图
pdf(file = "circos.pdf", width = 7, height = 7)
par(mar = c(2, 2, 2, 4))
circos.par(gap.degree = c(3, rep(2, nrow(cor1) - 1)), start.degree = 180)
chordDiagram(cor1, grid.col = rainbow(ncol(rt)), col = col1, transparency = 0.5, symmetric = TRUE)
par(xpd = TRUE)
colorlegend(col_vector, vertical = TRUE, labels = c(1, 0, -1), xlim = c(1.1, 1.3), ylim = c(-0.4, 0.4))
dev.off()
circos.clear()

# 绘制相关性矩阵图（与圈图配色一致）
pdf(file = "corrplot.pdf", width = 7, height = 7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col = colorRampPalette(c("#8BACD1", "white", "#C17F9E"))(50)
)
dev.off()
