setwd('F:/scRNA/workpath/2.thyrocytes')
Subset_thyrocytes <- readRDS("Subset_thyrocytes.rds")

# Extract SERPINA1 gene expression
SERPINA1_expression <- FetchData(Subset_thyrocytes, vars = c("SERPINA1"), cells = Cells(Subset_thyrocytes))

# Check the extraction result and ensure it is a numeric vector
if (is.list(SERPINA1_expression)) {
  SERPINA1_expression <- unlist(SERPINA1_expression)
}

SERPINA1_expression[is.na(SERPINA1_expression)] <- 0
SERPINA1_expression <- as.numeric(SERPINA1_expression)

# Extract cluster information
cluster_info <- Subset_thyrocytes@meta.data[["cellType"]]

# Create a data frame containing SERPINA1 expression and cluster information
data_df <- data.frame(SERPINA1_Expression = SERPINA1_expression, Cluster = cluster_info)

# Select specific cluster types
# selected_Clusters <- c("")
# data_df <- data_df[data_df$Tissue %in% selected_Clusters, ]

# Load necessary packages
library(ggplot2)
library(multcompView)

# Perform one-way ANOVA
anova_result <- aov(SERPINA1_Expression ~ Cluster, data = data_df)
summary(anova_result)

# If ANOVA is significant, perform Tukey HSD test
if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result, "Cluster")
  print(tukey_result)
  
  # Extract Tukey HSD results
  tukey_data <- as.data.frame(tukey_result$Cluster)
  
  # Only retain comparisons with "tumor_With_HT" and the other three groups
  # tukey_filtered <- tukey_data[grepl("tumor_With_HT", rownames(tukey_data)), ]
  
  # Add significance markers
  tukey_letters <- multcompLetters4(anova_result, tukey_result)
  tukey_significance <- tukey_letters$Cluster$Letters
}

# Calculate standard error for each group
SE <- with(data_df, tapply(SERPINA1_Expression, Cluster, function(x) sd(x) / sqrt(length(x))))

# Calculate mean expression for each cluster
mean_expression <- with(data_df, tapply(SERPINA1_Expression, Cluster, mean))

# Create a data frame with mean values and standard errors
mean_se_df <- data.frame(Cluster = names(mean_expression), Mean = mean_expression, SE = SE)

# Add significance markers to the data frame
mean_se_df$Significance <- tukey_significance[match(mean_se_df$Cluster, names(tukey_significance))]

# Plot bar chart
mean_se_df$Cluster <- factor(mean_se_df$Cluster)

# Plot ggplot figure
p3 <- ggplot(mean_se_df, aes(x = Cluster, y = Mean, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  xlab("Cluster") +
  ylab("Mean SERPINA1 Expression") +
  ggtitle("Mean SERPINA1 Expression Across Subtypes") +
  theme_minimal()

# Display the figure
print(p3)

pdf(file = "07.Cluster thyrocytes.combined.pdf", width = 5.5, height = 7)
p3
dev.off()

# Save the figure
# ggsave("08_SERPINA1_mean_expression_with_significance.pdf", width = 15, height = 10)