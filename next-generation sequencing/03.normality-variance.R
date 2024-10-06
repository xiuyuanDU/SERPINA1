dir <- setwd('F:/TCGA/workpath/01serpina1expression')
library(reshape2)

# Data preprocessing and normality assumption testing
# 'SERPINA1.txt' contains sample names and metric names
# 'SERPINA1group.txt' contains sample names and grouping information
# Read files, merge group information, and reshape data
serpina1 <- read.table('SERPINA1.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
group <- read.table('SERPINA1group.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
serpina1 <- melt(merge(serpina1, group, by = 'sampleID'), id = c('sampleID', 'group'))

# We want to check if there is a significant difference in observed species index between group1 and group2
# Selecting groups to compare
group_12 <- subset(serpina1, variable == 'SERPINA1' & group %in% c('Lymphocytic Thyroiditis', 'No Thyroiditis'))
group_12$group <- factor(group_12$group)
head(group_12, 10)

# Normality test
# Method 1: QQ plot to assess normality
library(car)

pdf("01.QQ-plot.pdf") 
qqPlot(lm(value ~ group, data = group_12), simulate = TRUE, main = 'QQ Plot', labels = FALSE)
dev.off() 

# Method 2: Shapiro-Wilk test; data is normally distributed if both p-values are greater than 0.05
shapiro <- tapply(group_12$value, group_12$group, shapiro.test)
shapiro
shapiro[['Lymphocytic Thyroiditis']]$p.value
shapiro[['Normal']]$p.value

# Homogeneity of variance test
# Levene's test
library(car)
serpina1$group <- as.factor(serpina1$group)
leveneTest(value ~ group, data = serpina1)

# Wilcoxon test
wilcox.test(value ~ group, data = serpina1, var.equal = TRUE)

wilcox.test(value ~ group, data = serpina1, var.equal = TRUE, exact = TRUE)

wilcox.test(value ~ group, data = serpina1, var.equal = TRUE, conf.int = TRUE)

library(ggpubr)
# Save as PDF format
pdf("01.ggboxplot_output.pdf", width = 7.5, height = 8)  
P <- ggboxplot(serpina1, x = "group", y = "value", fill = "group",
               alpha = 1,
               width = 0.5,
               legend = "none",
               notch = TRUE,
               font.legend = c(20, "bold", "black"),
               ylab = "SERPINA1", xlab = "",
               font.y = 30,
               add = "jitter", add.params = list(color = "grey"),
               x.text.angle = 0, y.text.angle = 90,
               font.tickslab = c(22, "plain", "black")
) +
  scale_x_discrete(labels = c('Lymphocytic Thyroiditis' = 'Lymphocytic thyroiditis', 'No Thyroiditis' = 'No Thyroiditis')) +
  ggtitle("Thyroid disorder history") + theme(plot.title = element_text(hjust = 0.5, size = 30))

library(rstatix)
stat.test <- serpina1 %>% wilcox_test(value ~ group) %>% add_significance()
stat.test <- stat.test %>% add_y_position()
P + stat_pvalue_manual(stat.test, label = "Wilcoxon, p = {p}", tip.length = 0.01, size = 10) + ylim(0, 20)

dev.off()  