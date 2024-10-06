#Remove potential doublets
setwd('F:/scRNA/workpath/DoulbletFinder')

NT=readRDS("NT.rds")
pT_WHT1=readRDS("pT_WHT1.rds")
pT_WOHT4=readRDS("pT_WOHT4.rds")
T_WHT1=readRDS("T_WHT1.rds")
T_WHT2=readRDS("T_WHT2.rds")
T_WHT3=readRDS("T_WHT3.rds")
T_WOHT4=readRDS("T_WOHT4.rds")
T_WOHT5=readRDS("T_WOHT5.rds")
T_WOHT6=readRDS("T_WOHT6.rds")


#View the current metadata columns
colnames(NT@meta.data)
#Retain key metadata columns
key_metadata <- c("orig.ident", "nCount_RNA", "nFeature_RNA","stim", "percent.mt","percent.rb")


NT@meta.data <- NT@meta.data[, key_metadata, drop = FALSE]
pT_WHT1@meta.data <- pT_WHT1@meta.data[, key_metadata, drop = FALSE]
pT_WOHT4@meta.data <- pT_WOHT4@meta.data[, key_metadata, drop = FALSE]
T_WHT1@meta.data <- T_WHT1@meta.data[, key_metadata, drop = FALSE]
T_WHT2@meta.data <- T_WHT2@meta.data[, key_metadata, drop = FALSE]
T_WHT3@meta.data <- T_WHT3@meta.data[, key_metadata, drop = FALSE]
T_WOHT4@meta.data <- T_WOHT4@meta.data[, key_metadata, drop = FALSE]
T_WOHT5@meta.data <- T_WOHT5@meta.data[, key_metadata, drop = FALSE]
T_WOHT6@meta.data <- T_WOHT6@meta.data[, key_metadata, drop = FALSE]


# View the default assay
DefaultAssay(NT)
# Modify the default assay, for example, set it to RNA
DefaultAssay(NT) <- "RNA"
DefaultAssay(pT_WHT1) <- "RNA"
DefaultAssay(pT_WOHT4) <- "RNA"
DefaultAssay(T_WHT1) <- "RNA"
DefaultAssay(T_WHT2) <- "RNA"
DefaultAssay(T_WHT3) <- "RNA"
DefaultAssay(T_WOHT4) <- "RNA"
DefaultAssay(T_WOHT5) <- "RNA"
DefaultAssay(T_WOHT6) <- "RNA"

#Remove unnecessary items
NT@assays[["SCT"]] <- NULL
pT_WHT1@assays[["SCT"]] <- NULL
pT_WOHT4@assays[["SCT"]] <- NULL
T_WHT1@assays[["SCT"]] <- NULL
T_WHT2@assays[["SCT"]] <- NULL
T_WHT3@assays[["SCT"]] <- NULL
T_WOHT4@assays[["SCT"]] <- NULL
T_WOHT5@assays[["SCT"]] <- NULL
T_WOHT6@assays[["SCT"]] <- NULL

# Merge data
thyroid.combined <- merge(NT, y = list(pT_WOHT4, pT_WHT1, T_WHT1, T_WHT2, T_WHT3, T_WOHT4, T_WOHT5, T_WOHT6), add.cell.ids = c("NT", "pT_WOHT4", "pT_WHT1", "T_WHT1", "T_WHT2", "T_WHT3", "T_WOHT4", "T_WOHT5", "T_WOHT6"), project = "thyroid")
thyroid.combined <- JoinLayers(thyroid.combined)

# Check for duplicate cells in the dataset and confirm there are none
duplicates <- duplicated(thyroid.combined@meta.data)
sum(duplicates)  # View how many duplicates there are

setwd('F:/scRNA/workpath')
saveRDS(thyroid.combined, file = "thyroid.combined.rds")





