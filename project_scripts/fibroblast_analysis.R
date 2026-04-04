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


meaningful_fibroblasts <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/fibroblast_cellstates_corrected.rds")
meaningful_fibroblasts <- SCTransform(meaningful_fibroblasts, assay = "originalexp")



# DNA damage markers
ddr_genes <- c("ATM", "ATR", "CHEK1", "CHEK2", "PARP1", "TP53", "MDC1")

genes_to_plot <- intersect(ddr_genes, rownames(meaningful_fibroblasts))

DotPlot(meaningful_fibroblasts, features = genes_to_plot) + 
  coord_flip() + 
  scale_color_gradientn(colors = c("lightgrey", "blue", "darkblue")) +
  theme_minimal() +
  labs(
    title = "DNA Damage Response Markers",
    subtitle = "Fibroblast states",
    x = "Gene",
    y = "Fibroblast State"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# pre-cafs
precaf_genes <- read.csv("/share/crsp/lab/dalawson/share/BRCA_NatGen_Objects/Kai_grant_analysis_10312023/preCAF_signature_100genes.csv")

precaf_genes <- precaf_genes$Names  
valid_precaf_genes <- list(intersect(precaf_genes, rownames(meaningful_fibroblasts)))

meaningful_fibroblasts <- AddModuleScore(
  object = meaningful_fibroblasts,
  features = valid_precaf_genes,
  name = "preCAF_Score"
)

colnames(meaningful_fibroblasts@meta.data)[which(colnames(meaningful_fibroblasts@meta.data) == "preCAF_Score1")] <- "preCAF_Score"
meaningful_fibroblasts$updated_states <- as.character(meaningful_fibroblasts$predicted.id)
threshold <- quantile(meaningful_fibroblasts$preCAF_Score, 0.80)
meaningful_fibroblasts$updated_states[meaningful_fibroblasts$preCAF_Score > threshold] <- "pre-CAF"

# plot precafs
DefaultAssay(meaningful_fibroblasts) <- "SCT"

meaningful_fibroblasts <- ScaleData(meaningful_fibroblasts)
meaningful_fibroblasts <- RunPCA(meaningful_fibroblasts)
meaningful_fibroblasts <- RunUMAP(meaningful_fibroblasts, dims = 1:15)
meaninful_fibroblasts <- FindNeighbors(meaningful_fibroblasts, dims = 1:15)
DimPlot(meaningful_fibroblasts, 
        group.by = "updated_states", 
        reduction = "umap", 
        repel = TRUE) +        
  theme_classic() + 
  labs(title = "Fibroblast States including pre-CAFs")

Idents(meaningful_fibroblasts) <- "updated_states"


DotPlot(meaningful_fibroblasts, 
        features = genes_to_plot, 
        split.by = "group", 
        cols = c("lightgrey", "blue", "red")) + 
  coord_flip() +
  theme_bw() +
  labs(
    title = "DNA Damage Response: Carrier vs Non-Carrier",
    subtitle = "Fibroblast states including pre-CAFs",
    x = "Repair Genes",
    y = "Cell State"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# NF-kB markers
nfkb_ddr_genes <- c(
  # Core NF-kB Complex
  "RELA", "NFKB1", "NFKB2", "RELB", 
  
  # Upstream Activators (NEMO/IKK complex linked to DDR)
  "IKBKG", "CHUK", "IKBKB", 
  
  # Downstream NF-kB targets / SASP effectors
  "IL6", "CXCL8", "CCL2", "TNFAIP3", "PTGS2"
)

genes_to_plot <- intersect(nfkb_ddr_genes, rownames(meaningful_fibroblasts))

DotPlot(meaningful_fibroblasts, 
        features = genes_to_plot, 
        split.by = "group", 
        cols = c("lightgrey", "blue", "red")) + 
  coord_flip() +
  theme_bw() +
  labs(
    title = "NF-kb: Carrier vs Non-Carrier",
    subtitle = "Fibroblast states including pre-CAFs",
    x = "NF-kb Genes",
    y = "Cell State"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# correlation between dna damage marker and NF-kb
meaningful_fibroblasts <- AddModuleScore(
  meaningful_fibroblasts, 
  features = list(intersect(ddr_genes, rownames(meaningful_fibroblasts))), 
  name = "DDR_Stress_Score"
)
meta <- meaningful_fibroblasts@meta.data
cor_val <- cor(meta$preCAF_Score, meta$DDR_Stress_Score1, method = "pearson")
ggplot(meta, aes(x = preCAF_Score, y = DDR_Stress_Score1, color = updated_states)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    title = "pre-CAF DNA damage Correlation",
    subtitle = paste("Pearson Correlation:", round(cor_val, 3)),
    x = "pre-CAF Module Score",
    y = "DDR Marker Score"
  )


library(SingleCellExperiment)

# milo analysis
sce <- as.SingleCellExperiment(meaningful_fibroblasts, assay = "SCT")
data_milo <- Milo(sce)

colData(data_milo)$cellID <- rownames(colData(data_milo))

data_milo <- buildGraph(data_milo, d = 5, k = 30)
data_milo <- makeNhoods(data_milo, k=60, d=5, prop = 0.15, refined=FALSE)
plotNhoodSizeHist(data_milo, bins=150)

data_milo <- countCells(
  data_milo,
  meta.data = as.data.frame(colData(data_milo)),
  sample    = "cellID"
)

colData(data_milo)$cellID <- NULL

meta_df <- as_tibble(as.data.frame(colData(data_milo)), rownames="cellID") %>%
  dplyr::select(cellID, group) %>%
  distinct() %>%
  mutate(
    group = factor(group, levels = c("Control", "BRCA1"))
  ) %>%
  column_to_rownames("cellID") 



design_df <- as_tibble(as.data.frame(colData(data_milo)), rownames="cellID") %>%
  dplyr::select(cellID, group) %>%
  distinct() %>%
  mutate(
    group = factor(group, levels = c("Control", "BRCA1"))
  ) %>%
  column_to_rownames("cellID")

meta_cells <- as_tibble(design_df, rownames = "cellID")
data_milo <- calcNhoodDistance(data_milo, d = 5)
saveRDS(data_milo, "/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/fibroblasts_milo.rds")

data_milo <- readRDS("/share/crsp/lab/dalawson/prabhakg/BRCA_fibroblasts/fibroblasts_milo.rds")


colData(data_milo)$sample <- colData(data_milo)$donor_id
data_milo <- countCells(
  data_milo,
  meta.data = as.data.frame(colData(data_milo)),
  sample    = "sample"
)

col_df <- as.data.frame(colData(data_milo)) %>%
  dplyr::select(sample, group)

design_df_correct <- col_df %>%
  group_by(sample) %>%
  summarise(
    comp_group = mean(group == "BRCA1")
  ) %>%
  column_to_rownames("sample")


cd <- as.data.frame(colData(data_milo)) %>%
  rownames_to_column("cellID") %>%
  mutate(group = factor(group, levels=c("Control","BRCA1")))



data_milo <- countCells(data_milo, meta.data = cd, sample = "sample")

design_df <- as.data.frame(colData(data_milo)) %>%
  rownames_to_column("cellID") %>%
  distinct(sample, group) %>%
  mutate(
    group = factor(as.character(group),
                         levels = c("Control","BRCA1"))  
  ) %>%
  filter(!is.na(group)) %>%              
  column_to_rownames("sample")

if (anyDuplicated(rownames(design_df))) {
  design_df <- design_df[!duplicated(rownames(design_df)), ]
}
stopifnot(all(colnames(nhoodCounts(data_milo)) %in% rownames(design_df)))
qr(model.matrix(~group, design_df))$rank  

milo_res <- testNhoods(data_milo, design = ~ group, design.df = design_df)


data_milo <- buildNhoodGraph(data_milo)
plotNhoodGraphDA(data_milo, milo_res, alpha = 0.1, size_range=c(2,6))

milo_res_modified <- milo_res

# Now groupNhoods should work
milo_res_endogroups <- groupNhoods(data_milo, milo_res_modified, 
                                   max.lfc.delta = 2, 
                                   overlap = 1)


p1 <- plotNhoodGroups(data_milo, milo_res_endogroups, 
                      size_range=c(1,3)) 

milo_res_endogroups <- annotateNhoods(data_milo, milo_res_endogroups, 'updated_states')

p2 <- plotDAbeeswarm(milo_res_endogroups, group.by = 'NhoodGroup') +
  facet_grid(updated_states~., scales="free", space="free")

alpha <- 0.10  
tot_cells <- Matrix::rowSums(nhoodCounts(data_milo))
keep_nhoods <- names(tot_cells)[tot_cells > 0]

df <- as.data.frame(milo_res_endogroups)
if (!"Nhood" %in% names(df)) df <- rownames_to_column(df, "Nhood")
df <- df %>% filter(Nhood %in% keep_nhoods)

# choose the adj p-value column present
sig_col <- dplyr::case_when(
  "SpatialFDR" %in% names(df) ~ "SpatialFDR",
  "PValueAdj"  %in% names(df) ~ "PValueAdj",
  TRUE                         ~ "PValue"
)

df <- df %>%
  mutate(
    status = case_when(
      .data[[sig_col]] < alpha & logFC > 0 ~ "Enriched in BRCA",
      .data[[sig_col]] < alpha & logFC < 0 ~ "Enriched in Control",
      TRUE                                 ~ "Not significant"
    ),
    status = factor(status,
                    levels = c("Enriched in BRCA","Not significant","Enriched in Control"))
  )

df$updated_states <- reorder(df$updated_states, df$logFC, FUN = median, na.rm = TRUE)

df$size <- tot_cells[df$Nhood]



L <- max(2, quantile(abs(df$logFC), 0.99, na.rm = TRUE)) 
p_grad <- ggplot(df, aes(x = logFC, y = updated_states)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_vline(xintercept = c(-2, 2), linetype = 3, linewidth = 0.3) +
  ggbeeswarm::geom_quasirandom(
    aes(color = logFC, size = size), 
    width = 0.2, stroke = 0
  ) +
  scale_alpha_manual(values = c(ns = 0.35, sig = 1), guide = "none") + 
  scale_size(range = c(1.2, 3.0), guide = "none") +
  scale_color_gradient2(
    low = "#3b4cc0", mid = "grey85", high = "#b40426", midpoint = 0,
    limits = c(-L, L), oob = squish, name = "logFC (Control vs BRCA)"
  ) +
  labs(x = "logFC (Control vs BRCA)", y = "Cell state") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
