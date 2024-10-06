# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(Seurat)

# Ensure that you have already executed previous relevant code
# Set cluster resolution
Subset_myeloid <- FindClusters(Subset_myeloid, resolution = 1)

# Define gene lists for dot plot
genes <- list( 
  "myeloid cells" = c("LYZ", "FCER1G", "TYROBP"),
  "CD16+/" = c("FCGR3A"),
  "CD14+ monocyte" = c("CD14", "S100A8", "S100A9"),
  "Macrophage" = c("CD68", "CSF1R", "MSR1"), # Corrected the typo here
  "M2" = c("MRC1", "CD163"),
  "M1" = c("IL1B", "CD86"),
  "mDCs" = c("LAMP3", "CD1C", "CD1E", "CCL17"),
  "SERPINA1" = c("SERPINA1"), # Removed leading space
  "pDCs" = c("IRF7", "IRF8", "SPIB"),
  "proliferating macrophages" = c("GTSE1", "UBE2C", "TOP2A")
)

# Create PDF for the plot
pdf(file = "10.ann_cluster_marker.pdf", width = 16, height = 10)

# Create DotPlot using Seurat's DotPlot function
dot_plot <- DotPlot(
  object = Subset_myeloid, 
  features = genes, 
  group.by = "RNA_snn_res.1", 
  dot.scale = 10,  # Adjust dot size
  scale.by = "count",  # or "radius" based on preference
  cols = brewer.pal(n = 9, name = "YlOrRd"),  # Using a color palette
  dot.min = 0.1 # minimum size of dots
)

# Add labels and format the plot
dot_plot + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    axis.title.x = element_text(face = "bold", size = 14),  
    axis.title.y = element_text(face = "bold", size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),  
    axis.text.y = element_text(face = "bold", size = 12),  
    legend.position = "bottom",  
    legend.title = element_text(face = "bold", size = 12),  
    legend.text = element_text(face = "bold", size = 12),  
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2),  
    panel.background = element_rect(fill = "lightgray"),  
    text = element_text(face = "bold")  
  )

# Close PDF device
dev.off()

