# Set working directory
dir <- setwd('F:/TCGA/workpath/01serpina1expression')

# Load openxlsx package
library(openxlsx)

# Read the first worksheet of the Excel file
sample_data <- read.xlsx("sample_with_SERPINA1.xlsx", sheet = 1)

# Create a new column "SERPINA1 expression"
sample_data$`SERPINA1 expression` <- ifelse(sample_data$SERPINA1 < 15, "low", "high")

# Save the updated data to the same Excel file
write.xlsx(sample_data, file = "572sample.xlsx", sheetName = "Sheet1", overwrite = TRUE)