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
options(future.globals.maxSize = 40 * 1024^3)  # 20 GiB



hbca_macrophages <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/HBCA_macrophages.rds")

counts <- GetAssayData(hbca_macrophages, layer = "counts")

ensembl_ids <- rownames(counts)
mapping <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL")

new_names <- ifelse(is.na(mapping), ensembl_ids, mapping)
new_names <- make.unique(new_names)

rownames(counts) <- new_names

hbca_macrophages_renamed <- CreateSeuratObject(
  counts = counts,
  meta.data = hbca_macrophages@meta.data
)

sum(grepl("^ENSG", rownames(hbca_macrophages_renamed)))  



hbca_macrophages <- SCTransform(hbca_macrophages_renamed)
DefaultAssay(hbca_macrophages) <- "SCT"


hbca_macrophages <- ScaleData(hbca_macrophages)
hbca_macrophages <- RunPCA(hbca_macrophages)
hbca_macrophages <- RunUMAP(hbca_macrophages, dims = 1:15)
hbca_macrophages <- FindNeighbors(hbca_macrophages, dims = 1:15)
hbca_macrophages <- FindClusters(hbca_macrophages, resolution = 0.2)

# label transfer
macrophage_ref <- readRDS("/share/crsp/lab/dalawson/share/HBCA_BRCA1_Immune_Analysis/myeloidcell.rds")
macrophage_ref <- SCTransform(macrophage_ref, verbose = TRUE)
DefaultAssay(macrophage_ref) <- "SCT"
common_features <- intersect(rownames(macrophage_ref[["SCT"]]), rownames(hbca_macrophages[["SCT"]]))


anchors <- FindTransferAnchors(
  reference = macrophage_ref, 
  query = hbca_macrophages, 
  normalization.method = "SCT",
  features = common_features,
  dims = 1:30,
  reference.reduction = "pca"
)
predictions <- TransferData(anchorset = anchors, refdata = macrophage_ref$cellstate, dims = 1:30)
hbca_macrophages <- AddMetaData(hbca_macrophages, metadata = predictions)

saveRDS(hbca_macrophages, "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/HBCA_macrophages_annotated.rds")