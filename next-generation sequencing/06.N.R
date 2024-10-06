setwd('F:/TCGA/workpath/03.N')
sample_data <- read.xlsx("572sample.xlsx", sheet = 1)
sample_data$pathologic_N <- substr(sample_data$pathologic_N, 1, 2)
sample_data <- sample_data[sample_data$pathologic_N != "NX", ]
str(sample_data)
write.xlsx(sample_data, file = "Nsample.xlsx", row.Names = FALSE)

library(reshape2)

# Data preprocessing and normality assumption testing
serpina1 <- read.table('pathologic_N.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
pathologic_N <- read.table('pathologic_Ngroup.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
serpina1 <- melt(merge(serpina1, pathologic_N, by = 'sampleID'), id = c('sampleID', 'pathologic_N'))

# Comparing groups
group_12 <- subset(serpina1, variable == 'SERPINA1' & pathologic_N %in% c('N0', 'N1'))
group_12$group <- factor(group_12$pathologic_N)
head(group_12, 10)

# Normality test using QQ plot
library(car)
pdf("03.QQ-plot.pdf") 
qqPlot(lm(value ~ pathologic_N, data = group_12), simulate = TRUE, main = 'QQ Plot(pathologic_N)', labels = FALSE)
dev.off() 

# Shapiro-Wilk test
shapiro <- tapply(group_12$value, group_12$pathologic_N, shapiro.test)
shapiro

# Levene's test for variance homogeneity
library(car)
serpina1$pathologic_N <- as.factor(serpina1$pathologic_N)
leveneTest(value ~ pathologic_N, data = serpina1)

# Wilcoxon test
wilcox.test(value ~ pathologic_N, data = serpina1, var.equal = FALSE)

library(ggpubr)
pdf("03.ggboxplot_output.pdf", width = 6, height = 8)  
P <- ggboxplot(serpina1, x = "pathologic_N", y = "value", fill = "pathologic_N", alpha = 1, width = 0.5, 
               legend = "none", notch = TRUE, 
               ylab = "SERPINA1", xlab = "pathologic_N", 
               font.y = 30, font.x = 30,
               add = "jitter", add.params = list(color = "grey"),
               x.text.angle = 0, y.text.angle = 90,
               font.tickslab = c(22, "plain", "black")) +
  scale_x_discrete(labels = c('N0' = 'N0', 'N1' = 'N1')) +
  ggtitle("Lymphocytic thyroiditis") + theme(plot.title = element_text(hjust = 0.5, size = 30)) + ylim(0, 20)

stat.test <- serpina1 %>% wilcox_test(value ~ pathologic_N) %>% add_significance()
stat.test <- stat.test %>% add_y_position()
P + stat_pvalue_manual(stat.test, label = "Wilcoxon, p = {p}", tip.length = 0.01, size = 10) 
dev.off()  

# Repeat the process for the normal data
serpina1 <- read.table('pathologic_N(NOMAL).txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
pathologic_N <- read.table('pathologic_Ngroup(NOMAL).txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
serpina1 <- melt(merge(serpina1, pathologic_N, by = 'sampleID'), id = c('sampleID', 'pathologic_N'))

group_12 <- subset(serpina1, variable == 'SERPINA1' & pathologic_N %in% c('N0', 'N1'))
group_12$group <- factor(group_12$pathologic_N)
head(group_12, 10)

# Normality test using QQ plot
pdf("03.QQ-plot(nomal).pdf") 
qqPlot(lm(value ~ pathologic_N, data = group_12), simulate = TRUE, main = 'QQ Plot(pathologic_N)', labels = FALSE)
dev.off() 

# Shapiro-Wilk test
shapiro <- tapply(group_12$value, group_12$pathologic_N, shapiro.test)
shapiro

# Levene's test for variance homogeneity
serpina1$pathologic_N <- as.factor(serpina1$pathologic_N)
leveneTest(value ~ pathologic_N, data = serpina1)

# Wilcoxon test
wilcox.test(value ~ pathologic_N, data = serpina1, var.equal = FALSE)

pdf("03.ggboxplot_output(NOMAL).txt.pdf", width = 6, height = 8)  
P <- ggboxplot(serpina1, x = "pathologic_N", y = "value", fill = "pathologic_N", alpha = 1, width = 0.5,
               legend = "none", notch = TRUE, 
               ylab = "SERPINA1", xlab = "pathologic_N", 
               font.y = 30, font.x = 30,
               add = "jitter", add.params = list(color = "grey"),
               x.text.angle = 0, y.text.angle = 90,
               font.tickslab = c(22, "plain", "black")) +
  scale_x_discrete(labels = c('N0' = 'N0', 'N1' = 'N1')) +
  ggtitle("No Thyroiditis") + theme(plot.title = element_text(hjust = 0.5, size = 30)) + ylim(0, 20)

stat.test <- serpina1 %>% wilcox_test(value ~ pathologic_N) %>% add_significance()
stat.test <- stat.test %>% add_y_position()
P + stat_pvalue_manual(stat.test, label = "Wilcoxon, p = {p}", tip.length = 0.01, size = 10)
dev.off()  

