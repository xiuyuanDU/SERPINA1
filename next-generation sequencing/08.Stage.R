setwd('F:/TCGA/workpath/06.Stage')
library(openxlsx)

# Read sample data
sample_data <- read.xlsx("572sample.xlsx", sheet = 1)

# Replace certain stages
sample_data$pathologic_stage <- ifelse(sample_data$pathologic_stage %in% c("Stage IVA", "Stage IVC"), "Stage IV", sample_data$pathologic_stage)

# Remove rows with NA in pathologic_stage
sample_data <- sample_data[!is.na(sample_data$pathologic_stage), ]

# Convert pathologic_stage to a factor with specified levels
sample_data$pathologic_stage <- factor(sample_data$pathologic_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

# Extract subsets
HT <- sample_data[sample_data$group == "Lymphocytic Thyroiditis", ]
normal <- sample_data[sample_data$group == "No Thyroiditis", ]

# Save subsets to Excel files
write.xlsx(HT, file = "HTsample_pathologic_stage.xlsx", rowNames = FALSE)
write.xlsx(normal, file = "NORMALsample_pathologic_stage.xlsx", rowNames = FALSE)

# Q-Q Plot for Lymphocytic Thyroiditis
sample_data <- read.xlsx("HTsample_pathologic_stage.xlsx", sheet = 1)
library(car)
pdf("HTQQ-plot(pathologic_stage).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$pathologic_stage), data = sample_data, simulate = TRUE, main = "Q-Q Plot (pathologic_stage)", labels = FALSE)
dev.off()

# Levene's test
library(car)
leveneTest(sample_data$SERPINA1 ~ sample_data$pathologic_stage, center = median)

# Kruskal-Wallis test
kruskal.test(sample_data$SERPINA1 ~ sample_data$pathologic_stage, sample_data)

# Choose color palette
library(RColorBrewer)
mycol <- brewer.pal(5, "Blues")
library(ggpubr)

P <- ggbarplot(sample_data, x = "pathologic_stage", y = "SERPINA1", fill = "pathologic_stage",
               add = c("mean_se", "dotplot"), legend = "none", palette = mycol) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))

# Pairwise group comparisons
library(rstatix)
stat.test <- sample_data %>% pairwise_wilcox_test(SERPINA1 ~ pathologic_stage, p.adjust.method = "bonf")
stat.test <- stat.test %>% add_y_position()

# Add p-values to the plot
P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.3)))

# Adding p-value comparison
P <- P + stat_compare_means(method = "kruskal.test", aes(label = "p.format"), label.x = 1, label.y = 28, size = 10)
P <- P + ggtitle("Lymphocytic thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30))

# Save as PDF
pdf("ggboxplot_output(pathologic_stage).pdf", width = 10, height = 8)  
P
dev.off() 

# Q-Q Plot for Normal
sample_data <- read.xlsx("NORMALsample_pathologic_stage.xlsx", sheet = 1)
library(car)
pdf("pathologic_stagesample_QQ-plot(ETE).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$pathologic_stage), data = sample_data, simulate = TRUE, main = "Q-Q Plot (pathologic_stage)", labels = FALSE)
dev.off()

# Levene's test for Normal
library(car)
leveneTest(sample_data$SERPINA1 ~ sample_data$pathologic_stage, center = median)

# Kruskal-Wallis test for Normal
kruskal.test(sample_data$SERPINA1 ~ sample_data$pathologic_stage, sample_data)

# Choose color palette
mycol <- brewer.pal(6, "Blues")
library(ggpubr)
sample_data$pathologic_stage <- factor(sample_data$pathologic_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

P <- ggbarplot(sample_data, x = "pathologic_stage", y = "SERPINA1", fill = "pathologic_stage",
               add = c("mean_se", "dotplot"), legend = "none", palette = mycol) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))

# Pairwise group comparisons for Normal
library(rstatix)
stat.test <- sample_data %>% pairwise_wilcox_test(SERPINA1 ~ pathologic_stage, p.adjust.method = "bonf")
stat.test <- stat.test %>% add_y_position()

# Add p-values to the Normal plot
P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.3)))

# Adding p-value comparison for Normal
P <- P + stat_compare_means(method = "kruskal.test", aes(label = "p.format"), label.x = 1, label.y = 29, size = 10)
P <- P + ggtitle("No Thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30))

# Save as PDF for Normal
pdf("ggboxplot_output(NORMAL).pdf", width = 10, height = 8)  
P
dev.off()


