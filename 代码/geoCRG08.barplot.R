# ===== 引用包 =====
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

# ===== 设置路径和文件名 =====
inputFile = "CIBERSORT-Results.txt"
setwd("C:/Users/41939/Desktop/大创项目文件夹/子宫内膜异位症/Reanalysis_GSE51981/12.barplot")

# ===== 读取免疫浸润数据 =====
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# ===== 样本分组信息提取 =====
con = grepl("_Control", rownames(rt), ignore.case = TRUE)
treat = grepl("_Treat", rownames(rt), ignore.case = TRUE)
conData = rt[con, ]
treatData = rt[treat, ]
conNum = nrow(conData)
treatNum = nrow(treatData)
data = t(rbind(conData, treatData))

# ===== 设置柱状图颜色 =====
mypalette <- colorRampPalette(brewer.pal(8, "Set2"))
col <- mypalette(nrow(data))
group_col_control <- "#66C2A5"
group_col_treat <- "#FC8D62"

# ===== 绘制 barplot =====
pdf(file = "barplot.pdf", width = 14.5, height = 8.5)
par(las = 1, mar = c(8, 5, 4, 16), mgp = c(3, 0.1, 0), cex.axis = 1.5)
a1 = barplot(data, col = col, xaxt = "n", yaxt = "n", ylab = "Relative Percent", cex.lab = 1.8)
a2 = axis(2, tick = FALSE, labels = FALSE)
axis(2, a2, paste0(a2 * 100, "%"))
par(srt = 0, xpd = TRUE)
rect(a1[1] - 0.5, -0.01, a1[conNum] + 0.5, -0.06, col = group_col_control)
text(a1[conNum] / 2, -0.035, "Control", cex = 1.8)
rect(a1[conNum] + 0.5, -0.01, a1[length(a1)] + 0.5, -0.06, col = group_col_treat)
text((a1[length(a1)] + a1[conNum]) / 2, -0.035, "Treat", cex = 1.8)
legend(par("usr")[2] * 0.98, par("usr")[4], legend = rownames(data), col = col, pch = 15, bty = "n", cex = 1.2)
dev.off()

# ===== 小提琴图部分（immune.diff.pdf）=====
# 转换为 ggplot2 格式
Type = gsub("(.*)_(.*)", "\\2", rownames(rt))
data2 = cbind(as.data.frame(t(data)), Type = Type)
data2 = melt(data2, id.vars = c("Type"))
colnames(data2) = c("Group", "CellType", "Fraction")
data2$Fraction <- as.numeric(data2$Fraction)
data2 <- na.omit(data2)

# 计算显著性符号位置
marker <- data2 %>%
  group_by(CellType, Group) %>%
  summarise(Max = max(Fraction), .groups = "drop") %>%
  arrange(CellType, desc(Max)) %>%
  distinct(CellType, .keep_all = TRUE)
marker$X <- 1:nrow(marker)
data2$CellType <- factor(data2$CellType, levels = marker$CellType)

# 绘图
p <- ggplot(data2, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_violin(scale = "width", alpha = 0.6, size = 0.5, width = 1) +
  geom_boxplot(width = 0, position = position_dodge(1), outlier.shape = NA, alpha = 0.1, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", shape = 16, size = 3, color = "white", position = position_dodge(1)) +
  scale_fill_manual(values = c("Control" = "#377EB8", "Treat" = "#E41A1C")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.75), breaks = seq(0, 0.6, 0.2)) +
  stat_compare_means(aes(group = Group),
                     method = "wilcox",
                     label = "p.signif",
                     label.y = marker$Max + 0.03,
                     size = 6,
                     symnum.args = list(cutpoint = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "NS"))) +
  geom_segment(data = marker,
               aes(x = X - 0.2, y = Max + 0.02,
                   xend = X + 0.2, yend = Max + 0.02),
               size = 0.8) +
  labs(x = NULL, y = "Relative Proportion", fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.border = element_rect(size = 1.2),
        legend.position = c(0.08, 0.85),
        legend.text = element_text(size = 16),
        legend.title = element_blank())

# 输出小提琴图
pdf(file = "immune.diff.pdf", width = 18, height = 10)
print(p)
dev.off()
