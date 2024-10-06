dir <- setwd('F:/TCGA/workpath/02T')
sample_data <- read.xlsx("No Thyroiditis.xlsx", sheet = 1)
# Keep the first two characters of each cell value
sample_data$pathologic_T <- substr(sample_data$pathologic_T, 1, 2)
sample_data <- sample_data[sample_data$pathologic_T != "TX", ]
sample_data$pathologic_T <- as.factor(sample_data$pathologic_T) # Convert character to factor
# Print the processed result
print(sample_data)
pathologic_T <- sample_data$pathologic_T
# Normality test and variance homogeneity test
library(car)
pdf("02.QQ-plot(pathologic_T).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$pathologic_T), # Data linear fitting with lm()
       data = sample_data, simulate = T,
       main = "Q-Q Plot(pathologic_T)", labels = F)
dev.off()
# Levene's test
library(car)
leveneTest(sample_data$SERPINA1 ~ sample_data$pathologic_T, center = median)

# Kruskal's test
kruskal.test(sample_data$SERPINA1 ~ sample_data$pathologic_T, sample_data)
# Select color palette
# install.packages("RColorBrewer")
library(RColorBrewer)
mycol <- brewer.pal(5, "Blues")
library(ggpubr)
P <- ggbarplot(sample_data, x = "pathologic_T", y = "SERPINA1", fill = "pathologic_T",
               add = c("mean_se", "dotplot"), # Add standard error and dot plot
               legend = "none", # Remove legend
               palette = mycol) + # Set color palette
  theme(axis.text.x = element_text(size = 20),  # Set x-axis label font size to 20
        axis.text.y = element_text(size = 20)) # Set y-axis label font size to 20

# Pairwise comparison between groups
library(rstatix)
pairwise_wilcox_test(sample_data, SERPINA1 ~ pathologic_T, p.adjust.method = "bonf") # Use Bonferroni method for p-value adjustment
# Organize p-values
stat.test <- sample_data %>%
  pairwise_wilcox_test(SERPINA1 ~ pathologic_T, p.adjust.method = "bonf")
stat.test <- stat.test %>% add_y_position()
# Plot with p-values
P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10)
# Hide non-significant p-values
P <- P + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15,                      
                            hide.ns = T) + # Hide non-significant p-values  
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.3))) # Adjust y-axis
# Add test p-value
P <- P + stat_compare_means(method = "kruskal.test", # Statistical method
                            aes(label = "p.format"), # Display format
                            label.x = 1, label.y = 28, # Position
                            size = 10) # Size
P <- P + ggtitle("No Thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30), # Center title at the top with size 30
        axis.title.x = element_text(size = 30), # Set x-axis title size
        axis.title.y = element_text(size = 30)) # Set y-axis title size

# Save as PDF format
pdf("01.ggboxplot_output(normal).pdf", width = 10, height = 8)  
P
dev.off() 

#########################################################################################
# In WPS, split the data into two groups by group
sample_data <- read.xlsx("Lymphocytic Thyroiditis.xlsx", sheet = 1)
# Keep the first two characters of each cell value
sample_data$pathologic_T <- substr(sample_data$pathologic_T, 1, 2)
sample_data <- sample_data[sample_data$pathologic_T != "TX", ]
sample_data$pathologic_T <- as.factor(sample_data$pathologic_T) # Convert character to factor
# Print the processed result
print(sample_data)
pathologic_T <- sample_data$pathologic_T
# Normality test and variance homogeneity test
library(car)
pdf("02.2.QQ-plot(pathologic_T).pdf") 
qqPlot(lm(sample_data$SERPINA1 ~ sample_data$pathologic_T), # Data linear fitting with lm()
       data = sample_data, simulate = T,
       main = "Q-Q Plot(pathologic_T)", labels = F)
dev.off()
# Levene's test
library(car)
leveneTest(sample_data$SERPINA1 ~ sample_data$pathologic_T, center = median)

# Kruskal's test
kruskal.test(sample_data$SERPINA1 ~ sample_data$pathologic_T, sample_data)
# Select color palette
# install.packages("RColorBrewer")
library(RColorBrewer)
mycol <- brewer.pal(5, "Blues")
library(ggpubr)
P <- ggbarplot(sample_data, x = "pathologic_T", y = "SERPINA1", fill = "pathologic_T",
               add = c("mean_se", "dotplot"), # Add standard error and dot plot
               legend = "none", # Remove legend
               palette = mycol) + # Set color palette
  theme(axis.text.x = element_text(size = 20),  # Set x-axis label font size to 20
        axis.text.y = element_text(size = 20)) # Set y-axis label font size to 20

# Pairwise comparison between groups
library(rstatix)
pairwise_wilcox_test(sample_data, SERPINA1 ~ pathologic_T, p.adjust.method = "bonf") # Use Bonferroni method for p-value adjustment
# Organize p-values
stat.test <- sample_data %>%
  pairwise_wilcox_test(SERPINA1 ~ pathologic_T, p.adjust.method = "bonf")
stat.test <- stat.test %>% add_y_position()
# Plot with p-values
P <- P + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.03, step.increase = 0.15, size = 10)
# Hide non-significant p-values
P <- P + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, step.increase = 0.15,                      
                            hide.ns = T) + # Hide non-significant p-values  
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.3))) # Adjust y-axis

# Add test p-value
P <- P + stat_compare_means(method = "kruskal.test", # Statistical method
                            aes(label = "p.format"), # Display format
                            label.x = 1, label.y = 28, # Position
                            size = 10) # Size
P <- P + ggtitle("Lymphocytic Thyroiditis") +
  theme(plot.title = element_text(hjust = 0.5, size = 30), # Center title at the top with size 30
        axis.title.x = element_text(size = 30), # Set x-axis title size
        axis.title.y = element_text(size = 30)) # Set y-axis title size

# Save as PDF format
pdf("01.ggboxplot_output(HT).pdf", width = 10, height = 8)  
P
dev.off() 
