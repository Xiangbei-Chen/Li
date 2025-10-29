# ===== 引用所需包 =====
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

# ===== 设置路径与文件 =====
expFile = "diffGeneExp.txt"           # 表达数据
immFile = "CIBERSORT-Results.txt"     # 免疫细胞浸润结果
setwd("C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/13.immuneCor")

# ===== 读取并整理表达数据 =====
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, -1]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)

# ===== 提取 Treat 组样本 =====
group = gsub("(.*)_(.*)", "\\2", colnames(data))
data = data[, group == "Treat", drop = FALSE]
data = t(data)

# ===== 读取免疫细胞浸润数据 =====
immune = read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# ===== 取交集样本，确保配对 =====
sameSample = intersect(rownames(data), rownames(immune))
data = data[sameSample, , drop = FALSE]
immune = immune[sameSample, , drop = FALSE]

# ===== 计算 Spearman 相关性 =====
outTab = data.frame()
for(cell in colnames(immune)) {
  if(sd(immune[, cell]) == 0) next
  for(gene in colnames(data)) {
    x = as.numeric(immune[, cell])
    y = as.numeric(data[, gene])
    corT = cor.test(x, y, method = "spearman")
    cor = corT$estimate
    pvalue = corT$p.value
    text = ifelse(pvalue < 0.001, "***",
                  ifelse(pvalue < 0.01, "**",
                         ifelse(pvalue < 0.05, "*", "")))
    outTab = rbind(outTab, data.frame(Gene = gene, Immune = cell, cor = cor, text = text, pvalue = pvalue))
  }
}

# ===== 绘制相关性热图（默认字体，无加粗）=====
outTab$cor = as.numeric(outTab$cor)

pdf(file = "cor.pdf", width = 8, height = 6)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  geom_text(aes(label = text), col = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", "*  p<0.05", "\n\nCorrelation")) +
  scale_x_discrete(position = "bottom")
dev.off()
