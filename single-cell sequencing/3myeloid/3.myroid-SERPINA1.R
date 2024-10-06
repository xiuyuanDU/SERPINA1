setwd('F:/scRNA/workpath/3.myeloid')
Subset_myeloid <- readRDS("Subset_myeloid.rds")

# Extract SERPINA1 gene expression
SERPINA1_expression <- FetchData(Subset_myeloid, vars = c("SERPINA1"), cells = Cells(Subset_myeloid))

# Check the extraction result and ensure it's a numeric vector
if (is.list(SERPINA1_expression)) {
  SERPINA1_expression <- unlist(SERPINA1_expression)
}

SERPINA1_expression[is.na(SERPINA1_expression)] <- 0
SERPINA1_expression <- as.numeric(SERPINA1_expression)

# Extract cluster information
cluster_info <- Subset_myeloid$seurat_clusters

# Create a data frame with SERPINA1 expression and cluster information
data_df <- data.frame(SERPINA1_Expression = SERPINA1_expression, Cluster = cluster_info)

# Load necessary packages
library(ggplot2)
library(multcompView)
seurat_clusters <- Subset_myeloid@meta.data[["seurat_clusters"]]

# Perform one-way ANOVA
anova_result <- aov(SERPINA1_Expression ~ seurat_clusters, data = data_df)
summary(anova_result)

# If ANOVA is significant, perform Tukey HSD test
if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result, "seurat_clusters")
  print(tukey_result)
  
  # Extract Tukey HSD results
  tukey_data <- as.data.frame(tukey_result$seurat_clusters)
  
  # Add significance markers
  tukey_letters <- multcompLetters4(anova_result, tukey_result)
  tukey_significance <- tukey_letters$seurat_clusters$Letters
}

# Calculate standard error for each group
SE <- with(data_df, tapply(SERPINA1_Expression, seurat_clusters, function(x) sd(x) / sqrt(length(x))))

# Calculate mean expression for each cluster
mean_expression <- with(data_df, tapply(SERPINA1_Expression, seurat_clusters, mean))

# Create a data frame with means and standard errors
mean_se_df <- data.frame(seurat_clusters = names(mean_expression), Mean = mean_expression, SE = SE)

# Add significance markers to the data frame
mean_se_df$Significance <- tukey_significance[match(mean_se_df$seurat_clusters, names(tukey_significance))]

# Plot the bar chart
mean_se_df$seurat_clusters <- factor(mean_se_df$seurat_clusters, levels = 0:15)

# Create ggplot
p <- ggplot(mean_se_df, aes(x = seurat_clusters, y = Mean, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  xlab("Cell Clusters") +
  ylab("Mean SERPINA1 Expression") +
  ggtitle("Comparison of Mean SERPINA1 Expression Across Cell Clusters") +
  theme_minimal() +
  scale_x_discrete(labels = 0:15)

# Display the plot
print(p)

pdf(file = "10.Subset_myeloid_SERPINA1_Expression.pdf", width = 8, height = 5)
p
dev.off()