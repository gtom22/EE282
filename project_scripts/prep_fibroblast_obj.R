library(Seurat)
library(dplyr)
library(harmony)
library(readxl)
library(ggplot2)
library(miloR)
library(tidyverse)
library(patchwork)
library(igraph)
library(RColorBrewer)
library(cowplot)
library(scater)
library(scran)
library(tibble)
library(ComplexHeatmap)
library(ggbeeswarm) 
library(scales)
library(tidyr)
library(circlize)
library(biomaRt)
library(pheatmap)
library(future)
library(org.Hs.eg.db)
library(EnhancedVolcano)

options(future.globals.maxSize = 50 * 1024^3) 

brca <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/combined_BRCA.rds")
brca_fibroblasts <- subset(brca, subset = cell_type == "fibroblast of mammary gland")
hbca_fibroblasts <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/HBCA_fibroblasts_fixed.rds")

set.seed(123)
cells_to_keep <- sample(colnames(hbca_fibroblasts), size = 50000)
hbca_fibroblasts <- hbca_fibroblasts[, cells_to_keep]



counts <- GetAssayData(hbca_fibroblasts, layer = "counts")

ensembl_ids <- rownames(counts)
mapping <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL")

new_names <- ifelse(is.na(mapping), ensembl_ids, mapping)
new_names <- make.unique(new_names)

rownames(counts) <- new_names

hbca_fibroblasts_renamed <- CreateSeuratObject(
  counts = counts,
  meta.data = hbca_fibroblasts@meta.data
)

sum(grepl("^ENSG", rownames(hbca_fibroblasts_renamed))) 


hbca_fibroblasts <- SCTransform(hbca_fibroblasts_renamed)


saveRDS(hbca_fibroblasts, "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/HBCA_fibroblasts_corrected_SCT.rds")
objs <- list(brca_fibroblasts, hbca_fibroblasts)

usable_features <- intersect(rownames(brca_fibroblasts), rownames(hbca_fibroblasts))
length(usable_features)

features <- usable_features
features <- usable_features[1:min(6000, length(usable_features))]

markers <- c("LUM", "MMP3", "MMP12")
features <- usable_features[1:min(5000, length(usable_features))]

features <- union(features, markers[markers %in% usable_features])

sct_genes_1 <- rownames(objs[[1]][["SCT"]])
sct_genes_2 <- rownames(objs[[2]][["SCT"]])
features <- intersect(features, sct_genes_1)
features <- intersect(features, sct_genes_2)
features <- as.character(features)
features <- features[!is.na(features)]
features <- unique(features)

# Verify markers made it through
print(markers %in% features)
fibroblast_markers <- c("LUM", "MMP3", "MMP12", "COL1A1", "COL1A2", "COL3A1", 
                        "DCN", "COL6A1", "COL6A2", "POSTN", "THBS1", "MMP1", 
                        "MMP10", "MMP12", "ACTA2", "TAGLN")
for (marker in fibroblast_markers) {
  if (marker %in% sct_genes_1 | marker %in% sct_genes_2) {
    features <- c(features, marker)
  }
}

features <- unique(features)
features <- as.character(features)
features <- features[!is.na(features)]

hbca_fibroblasts <- subset(hbca_fibroblasts, features = features)
brca_fibroblasts <- subset(brca_fibroblasts, features = features)

VariableFeatures(hbca_fibroblasts) <- features
VariableFeatures(brca_fibroblasts) <- features



merged <- merge(
  x = hbca_fibroblasts,
  y = brca_fibroblasts
)


DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- features

merged <- ScaleData(merged, verbose = FALSE)

merged <- RunPCA(merged, verbose = FALSE)

DefaultAssay(merged) <- "SCT"
integrated <- IntegrateLayers(
  object = merged,
  method = HarmonyIntegration
)

saveRDS(integrated, "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/integrated.rds")

fibroblasts <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/integrated.rds")


# filter the fibroblasts
cells_to_keep <- is.na(fibroblasts$risk_status) | fibroblasts$risk_status %in% c("AR", "HR-BR1", "HR-cBR1")
subsetted <- fibroblasts[, cells_to_keep]
subsetted$risk_status_updated <- subsetted$risk_status

subsetted$risk_status_updated[is.na(subsetted$risk_status) & 
                                !is.na(subsetted$custom_labels)] <- "Other"
other_cells <- subsetted[, subsetted$risk_status_updated == "Other"]
subsetted_for_matching <- subsetted[, subsetted$risk_status_updated %in% c("AR", "HR-BR1", "HR-cBR1")]


metadata <- subsetted_for_matching@meta.data
patients <- unique(metadata[, c("donor_id", "donor_age", "risk_status_updated")])


# Get age distribution of HR-BR1 patients (remove those with NA age)
hr_br1_patients <- patients[patients$risk_status_updated == "HR-BR1" & 
                              !is.na(patients$risk_status_updated) & 
                              !is.na(patients$donor_age), ]
hr_br1_ages <- hr_br1_patients$donor_age


age_min <- min(hr_br1_ages, na.rm = TRUE)
age_max <- max(hr_br1_ages, na.rm = TRUE)
print(summary(hr_br1_ages))

ar_patients <- patients[
  patients$risk_status_updated == "AR" &
    !is.na(patients$donor_age),
]

ar_matched <- ar_patients[ar_patients$donor_age >= age_min & 
                            ar_patients$donor_age <= age_max, ]


meta <- fibroblasts@meta.data

donor_level <- unique(meta[, c("donor_id", "donor_age", "risk_status")])



# Combine with your previously defined patients
final_patients <- unique(rbind(hr_br1_patients,
                               ar_matched,
                               nee_hrcbr1))


# Combine the final patient set
print(table(final_patients$risk_status, useNA = "ifany"))

meta_full <- fibroblasts@meta.data
donor_level <- unique(meta_full[, c("donor_id", "donor_age", "risk_status")])

nee_donors_all <- unique(donor_level$donor_id[grepl("^Nee", donor_level$donor_id)])

keep_donors <- unique(c(hr_br1_patients$donor_id,
                        ar_matched$donor_id,
                        nee_donors_all))

final_subset <- subsetted_for_matching[, subsetted_for_matching$donor_id %in% keep_donors]
final_subset_with_other <- merge(final_subset, other_cells)
fibroblasts <- final_subset_with_other

fibroblasts$group <- NA

# For cells with brca_status (HBCA dataset)
brca_status <- c("HR-BR1", "HR-cBR1")
fibroblasts$group[fibroblasts$risk_status %in% brca_status] <- "BRCA1"
fibroblasts$group[fibroblasts$risk_status == "AR"] <- "Control"

# For cells with sample (BRCA dataset)
fibroblasts$group[grepl("_BRCA", fibroblasts$sample)] <- "BRCA1"
fibroblasts$group[grepl("_Ctrl", fibroblasts$sample)] <- "Control"

# Drop unknown and NA
fibroblasts <- subset(fibroblasts, subset = !is.na(group))


# filter out Kumar BRCA
is_kumar <- grepl("^Kumar", fibroblasts$donor_id)

# Identify Kumar BRCA1
remove_cells <- is_kumar & fibroblasts$group == "BRCA1"

# Subset object
fibroblasts <- subset(
  fibroblasts,
  cells = colnames(fibroblasts)[!remove_cells]
)

# Verify
print("Unique groups:")
print(unique(fibroblasts$group))
print("\nGroup counts:")
print(table(fibroblasts$group))


fibroblasts <- FindVariableFeatures(fibroblasts, verbose = FALSE)
fibroblasts <- ScaleData(fibroblasts, verbose = FALSE)
fibroblasts <- RunPCA(fibroblasts, verbose = FALSE)


fibroblasts <- RunUMAP(fibroblasts, dims = 1:15)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.2)

DimPlot(fibroblasts, label = TRUE)
DimPlot(fibroblasts, group.by = "group")
DimPlot(fibroblasts, group.by = "donor_id")



# Cluster 3, 4, 7 from one patient each
predefined_markers <- unique(c(
  # Fibro-prematrix
  "CFD", "FOS", "GPX3", "PLA2G2A", "WISP2", "MGST1", "MFAP5",
  
  # Fibro-SFRP4
  "SFRP4", "CLU", "G0S2", "MGP", "OGN", "ADIRF", "IGFBP5",
  
  # Fibro-major
  "MMP3", "CXCL1", "GEM", "KDM6B", "CEBPB", "TNFAIP6", "CXCL2",
  
  # Fibro-matrix
  "COL3A1", "POSTN", "COL1A1", "IGF1", "IGFBP2", "ADAM12", "TNC"
)





desired_clusters <- c("0", "1", "2", "3", "4", "5")
meaningful_fibroblasts <- subset(fibroblasts, subset = seurat_clusters %in% desired_clusters)

meaningful_fibroblasts <- RunUMAP(meaningful_fibroblasts, dims = 1:15)
meaningful_fibroblasts <- FindNeighbors(meaningful_fibroblasts, dims = 1:15)
meaningful_fibroblasts <- FindClusters(meaningful_fibroblasts, resolution = 0.1)

DimPlot(meaningful_fibroblasts, group.by = "seurat_clusters", label = TRUE)

desired_clusters <- c("0", "1", "2")
meaningful_fibroblasts <- subset(meaningful_fibroblasts, subset = seurat_clusters %in% desired_clusters)


# # redo integration
DefaultAssay(meaningful_fibroblasts) <- "SCT"   

meaningful_fibroblasts <- NormalizeData(meaningful_fibroblasts)
meaningful_fibroblasts <- FindVariableFeatures(meaningful_fibroblasts)
meaningful_fibroblasts <- ScaleData(meaningful_fibroblasts)
meaningful_fibroblasts <- RunPCA(meaningful_fibroblasts, npcs = 30)
meaningful_fibroblasts <- harmony::RunHarmony(meaningful_fibroblasts, group.by.vars = "donor_id")

meaningful_fibroblasts <- RunUMAP(meaningful_fibroblasts, reduction = "harmony", dims = 1:15)
meaningful_fibroblasts <- FindNeighbors(meaningful_fibroblasts, reduction = "harmony", dims = 1:15)
meaningful_fibroblasts <- FindClusters(meaningful_fibroblasts, resolution = 0.2)

desired_clusters <- c("0", "1", "2", "3")
meaningful_fibroblasts <- subset(meaningful_fibroblasts, subset = seurat_clusters %in% desired_clusters)


# # redo integration
DefaultAssay(meaningful_fibroblasts) <- "SCT"  
meaningful_fibroblasts <- NormalizeData(meaningful_fibroblasts)
meaningful_fibroblasts <- FindVariableFeatures(meaningful_fibroblasts)
meaningful_fibroblasts <- ScaleData(meaningful_fibroblasts)
meaningful_fibroblasts <- RunPCA(meaningful_fibroblasts, npcs = 30)
meaningful_fibroblasts <- harmony::RunHarmony(meaningful_fibroblasts, group.by.vars = "donor_id")

meaningful_fibroblasts <- RunUMAP(meaningful_fibroblasts, reduction = "harmony", dims = 1:15)
meaningful_fibroblasts <- FindNeighbors(meaningful_fibroblasts, reduction = "harmony", dims = 1:15)
meaningful_fibroblasts <- FindClusters(meaningful_fibroblasts, resolution = 0.2)

# Check the results
DimPlot(meaningful_fibroblasts, group.by = "seurat_clusters", label = TRUE)
DimPlot(meaningful_fibroblasts, group.by = "donor_id")
DimPlot(meaningful_fibroblasts, group.by = "group")


### Perform label transfer from fibroblast cluster labels (kumar paper)
fibroblast_cellstates <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/fibroblast_cellstates.rds")

counts <- GetAssayData(fibroblast_cellstates, layer = "counts")

ensembl_ids <- rownames(counts)
mapping <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL")

new_names <- ifelse(is.na(mapping), ensembl_ids, mapping)
new_names <- make.unique(new_names)

rownames(counts) <- new_names

hbca_fibroblasts_renamed <- CreateSeuratObject(
  counts = counts,
  meta.data = fibroblast_cellstates@meta.data
)



fibroblast_cellstates <- hbca_fibroblasts_renamed




Idents(fibroblast_cellstates) <- "author_cell_type"
fibroblast_cellstates$fib_cluster <- Idents(fibroblast_cellstates)
DefaultAssay(fibroblast_cellstates) <- "RNA"
DefaultAssay(meaningful_fibroblasts) <- "RNA"

fibroblast_cellstates <- NormalizeData(fibroblast_cellstates, verbose = FALSE)
fibroblast_cellstates <- FindVariableFeatures(fibroblast_cellstates, nfeatures = 3000, verbose = FALSE)

meaningful_fibroblasts <- NormalizeData(meaningful_fibroblasts, verbose = FALSE)
meaningful_fibroblasts <- FindVariableFeatures(meaningful_fibroblasts, nfeatures = 3000, verbose = FALSE)
features <- intersect(rownames(fibroblast_cellstates[["RNA"]]), rownames(meaningful_fibroblasts[["RNA"]]))
features <- head(features, 3000)


anchors <- FindTransferAnchors(
  reference = fibroblast_cellstates,
  query = meaningful_fibroblasts,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  features = features,
  dims = 1:20,
  k.anchor = 20
)
ref.labels <- setNames(fibroblast_cellstates$fib_cluster, colnames(fibroblast_cellstates))

predictions <- TransferData(
  anchorset = anchors,
  refdata = ref.labels,
  dims = 1:20,
  k.weight = 20
)
meaningful_fibroblasts <- AddMetaData(meaningful_fibroblasts, predictions)

# incorrect cell state assignment transition to scvi approach
saveRDS(meaningful_fibroblasts, "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/ref_fibroblasts_cellstates.rds")
