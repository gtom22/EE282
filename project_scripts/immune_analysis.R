library(Seurat)
library(future)
library(org.Hs.eg.db)
library(tidyr)
library(ggplot2)




hbca_immune <- readRDS("/share/crsp/lab/dalawson/prabhakg/immune_ihbca/macro_subcluster.rds")


hbca_immune$group <- ifelse(hbca_immune$risk_status == "HR-BR1", "BRCA1",
                                 ifelse(hbca_immune$risk_status == "AR", "Control", NA))

hbca_immune <- subset(hbca_immune, subset = group %in% c("BRCA1", "Control"))


# Myeloid Exhaustion / Immunosuppression
exhaustion_features <- c("CD274", "PDCD1LG2", "IDO1", "LAIR1", "HAVCR2") 

# Interferon Stimulated Genes (ISGs)
isg_features <- c("IFIT1", "IFIT2", "IFIT3", "ISG15", "OAS1", "MX1")

hbca_immune <- AddModuleScore(hbca_immune, features = list(exhaustion_features), name = "Exhaustion_Score")
hbca_immune <- AddModuleScore(hbca_immune, features = list(isg_features), name = "ISG_Score")

VlnPlot(hbca_immune, 
        features = "Exhaustion_Score1", 
        group.by = "predicted.id", 
        split.by = "group", 
        pt.size = 0) + 
  ggtitle("Myeloid Exhaustion: BRCA vs Control")


VlnPlot(hbca_immune, 
        features = "ISG_Score1", 
        group.by = "predicted.id", 
        split.by = "group", 
        pt.size = 0) + 
  ggtitle("Interferon: BRCA vs Control")

# milo analysis
hbca_immune$donor_id <- droplevels(as.factor(hbca_immune$donor_id))
hbca_immune <- subset(hbca_immune, subset = donor_id != "" & !is.na(donor_id))
sce <- as.SingleCellExperiment(hbca_immune, assay = "SCT")
data_milo <- Milo(sce)

colData(data_milo)$cellID <- rownames(colData(data_milo))

data_milo <- buildGraph(data_milo, d = 30, k = 30)
data_milo <- makeNhoods(data_milo, k=60, d=5, prop = 0.15, refined=TRUE)
plotNhoodSizeHist(data_milo, bins=150)

data_milo <- countCells(
  data_milo, 
  meta.data = as.data.frame(colData(data_milo)), 
  sample = "cellID"
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



colData(data_milo)$sample <- colData(data_milo)$donor_id
data_milo <- countCells(
  data_milo,
  meta.data = as.data.frame(colData(data_milo)),
  sample    = "sample"
)

col_df <- as.data.frame(colData(data_milo)) %>%
  dplyr::select(sample, group)

# aggregate per sample: compute fraction of metastatic niche
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


milo_res_endogroups <- groupNhoods(data_milo, milo_res_modified, 
                                   max.lfc.delta = 2, 
                                   overlap = 1)


p1 <- plotNhoodGroups(data_milo, milo_res_endogroups, 
                      size_range=c(1,3)) 

milo_res_endogroups <- annotateNhoods(data_milo, milo_res_endogroups, 'predicted.id')

p2 <- plotDAbeeswarm(milo_res_endogroups, group.by = 'NhoodGroup') +
  facet_grid(predicted.id~., scales="free", space="free")

alpha <- 0.10  # FDR cutoff
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

# reorder cell types by median logFC for nicer reading
df$predicted.id <- reorder(df$predicted.id, df$logFC, FUN = median, na.rm = TRUE)

# size aesthetic = neighborhood size
df$size <- tot_cells[df$Nhood]



# color scale based on fold change
L <- max(2, quantile(abs(df$logFC), 0.99, na.rm = TRUE)) 
p_grad <- ggplot(df, aes(x = logFC, y = predicted.id)) +
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
