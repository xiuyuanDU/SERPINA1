LRP1_expression <- FetchData(Subset_thyrocytes, vars = c("LRP1"), cells = Cells(Subset_thyrocytes))
if (is.list(LRP1_expression)) {
  LRP1_expression <- unlist(LRP1_expression)
}

LRP1_expression[is.na(LRP1_expression)] <- 0
LRP1_expression <- as.numeric(LRP1_expression)

cluster_info <- Subset_thyrocytes$seurat_clusters

data_df <- data.frame(LRP1_expression = LRP1_expression, Cluster = cluster_info)

# selected_Clusters <- c("")
# data_df <- data_df[data_df$Tissue %in% selected_Clusters, ]

library(ggplot2)
library(multcompView)

anova_result <- aov(LRP1_expression ~ Cluster, data = data_df)
summary(anova_result)

if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
  tukey_result <- TukeyHSD(anova_result, "Cluster")
  print(tukey_result)
  
  tukey_data <- as.data.frame(tukey_result$Cluster)
  
  # tukey_filtered <- tukey_data[grepl("tumor_With_HT", rownames(tukey_data)), ]
  
  tukey_letters <- multcompLetters4(anova_result, tukey_result)
  tukey_significance <- tukey_letters$Cluster$Letters
}

SE <- with(data_df, tapply(LRP1_expression, Cluster, function(x) sd(x) / sqrt(length(x))))

mean_expression <- with(data_df, tapply(LRP1_expression, Cluster, mean))

mean_se_df <- data.frame(Cluster = names(mean_expression), Mean = mean_expression, SE = SE)

mean_se_df$Significance <- tukey_significance[match(mean_se_df$Cluster, names(tukey_significance))]

mean_se_df$Cluster <- factor(mean_se_df$Cluster, levels = 0:12)

p3 <- ggplot(mean_se_df, aes(x = Cluster, y = Mean, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  xlab("Cluster") +
  ylab("Mean LRP1 expression") +
  ggtitle("Comparison of Mean LRP1 Expression Across Cell Clusters") +
  theme_minimal() +
  scale_x_discrete(labels = 0:12)

print(p3)

pdf(file = "07.LRP1_expression.pdf", width = 7, height = 5)
p3
dev.off()
# ggsave("08_SERPINA1_mean_expression_with_significance.pdf", width = 15, height = 10)