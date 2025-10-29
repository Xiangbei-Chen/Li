# install.packages("BiocManager")
# BiocManager::install("limma")

library(limma)

# ==== 设置路径与文件 ====
expFile <- "geneMatrix.txt"
conFile <- "s1.txt"
treatFile <- "s2.txt"
setwd("C://Users//41939//Desktop//大创项目文件夹//子宫内膜异位症//Reanalysis_GSE51981//05.normalize-GSE51981")

# ==== 读取表达数据 ====
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# ==== 合并重复基因名 ====
data <- avereps(data)

# ==== 读取样本组信息 ====
s1 <- read.table(conFile, header = FALSE, sep = "\t", check.names = FALSE)
sampleName1 <- as.vector(s1[, 1])
conData <- data[, sampleName1]

s2 <- read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
sampleName2 <- as.vector(s2[, 1])
treatData <- data[, sampleName2]

# ==== 合并表达矩阵 ====
rt <- cbind(conData, treatData)

# ==== 构建带分组标签的列名 ====
group_labels <- c(rep("Control", length(sampleName1)), rep("Treat", length(sampleName2)))
tagged_colnames <- paste0(colnames(rt), "_", group_labels)

# ==== 判断是否为 log2 转换 ====
stats <- summary(as.numeric(as.matrix(rt)))
print(stats)

# ==== 写出文件：包含标签行（作为第一行） ====
if (stats["Max."] < 25) {
  message("🟡 表达数据为 log2 转换，将执行 2^x 还原...")
  
  expr_data <- 2^rt
  
  # 添加分组标签行
  expr_out <- rbind(id = tagged_colnames, expr_data)
  norm_out <- rbind(id = tagged_colnames, rt)
  
  write.table(expr_out, file = "raw_expr.txt", sep = "\t", quote = FALSE, col.names = FALSE)
  write.table(norm_out, file = "normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)
  
} else {
  message("🟢 表达数据为线性表达（非 log 转换），直接保存为 normalize.txt")
  
  norm_out <- rbind(id = tagged_colnames, rt)
  write.table(norm_out, file = "normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)
}
