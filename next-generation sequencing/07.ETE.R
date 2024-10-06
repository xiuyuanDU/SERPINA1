setwd('F:/TCGA/workpath/05.ETE')
library(openxlsx)
sample_data <- read.xlsx("572sample.xlsx", sheet = 1)
sample_data$ETE <- substr(sample_data$extrathyroid_carcinoma_present_extension_status, 1, 4)
sample_data <- sample_data[!(is.na(sample_data$ETE) | sample_data$ETE == "Very"), ]
str(sample_data)
write.xlsx(sample_data, file = "570.sample.xlsx", row.Names = FALSE)
sample_data <- read.xlsx("570.sample.xlsx", sheet = 1)
sample_data$ETE <- as.factor(sample_data$ETE)
ETE <- sample_data$ETE

# Extract subsets
HT <- sample_data[sample_data$group == "Lymphocytic Thyroiditis", ]
normal <- sample_data[sample_data$group == "No Thyroiditis", ]
write.xlsx(HT, file = "HTsmple_ETE.xlsx", rowNames = FALSE)
write.xlsx(normal, file = "NORMALsmple_ETE.xlsx", rowNames = FALSE)

# QQ plot for Lymphocytic Thyroiditis
sample_data <- read.xlsx("HTsmple_ETE.xlsx", sheet = 1)
library(car)
pdf("HTQQ-plot(ETE).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$ETE), data = sample_data, simulate = TRUE,
       main = "Q-Q Plot(ETE)", labels = FALSE)
dev.off()

# Levene test
library(car)
leveneTest(sample_data$SERPINA1 ~ sample_data$ETE, center = median)
# Kruskal test
kruskal.test(sample_data$SERPINA1 ~ sample_data$ETE, sample_data)

# Color palette
library(RColorBrewer)
mycol <- brewer.pal(5, "Blues")
library(ggpubr)
P <- ggbarplot(sample_data, x = "ETE", y = "SERPINA1", fill = "ETE",
               add = c("mean_se", "dotplot"),
               legend = "none", 
               palette = mycol) +
  theme(axis.text.x = element_text(size = 20),  
        axis.text.y = element_text(size = 20))

# Pairwise comparisons
library(rstatix)
stat.test <- sample_data %>%
  pairwise_wilcox_test(SERPINA1 ~ ETE, p.adjust.method = "bonf") %>%
  add_y_position()

P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15, hide.ns = TRUE)

P <- P + stat_compare_means(method = "kruskal.test",
                            aes(label = "p.format"),
                            label.x = 1, label.y = 28, size = 10) +
  ggtitle("Lymphocytic thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

pdf("ggboxplot_output(HT).pdf", width = 10, height = 8)  
P
dev.off()

# QQ plot for No Thyroiditis
sample_data <- read.xlsx("NORMALsmple_ETE.xlsx", sheet = 1)
pdf("NORMALsmple_QQ-plot(ETE).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$ETE), data = sample_data, simulate = TRUE,
       main = "Q-Q Plot(ETE)", labels = FALSE)
dev.off()

# Levene test for No Thyroiditis
leveneTest(sample_data$SERPINA1 ~ sample_data$ETE, center = median)
kruskal.test(sample_data$SERPINA1 ~ sample_data$ETE, sample_data)

# Color palette and bar plot for No Thyroiditis
P <- ggbarplot(sample_data, x = "ETE", y = "SERPINA1", fill = "ETE",
               add = c("mean_se", "dotplot"), 
               legend = "none", 
               palette = mycol) +
  theme(axis.text.x = element_text(size = 20),  
        axis.text.y = element_text(size = 20))

stat.test <- sample_data %>%
  pairwise_wilcox_test(SERPINA1 ~ ETE, p.adjust.method = "bonf") %>%
  add_y_position()

P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15, hide.ns = TRUE)

P <- P + stat_compare_means(method = "kruskal.test", 
                            aes(label = "p.format"), 
                            label.x = 1, label.y = 28, size = 10) +
  ggtitle("No Thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

pdf("ggboxplot_output(NORMAL).pdf", width = 10, height = 8)  
P
dev.off()


