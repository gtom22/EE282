library(Seurat)


parent <- "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/"

dirs <- list.dirs(path = parent, full.names = TRUE)
#sample_dirs <- dirs[dir.exists(file.path(dirs, "outs", "filtered_feature_bc_matrix"))]
#matrix_dirs <- file.path(sample_dirs, "outs", "filtered_feature_bc_matrix")
sample_dirs <- dirs[dir.exists(file.path(dirs, "outs", "filtered_feature_bc_matrix"))]
matrix_dirs <- normalizePath(file.path(sample_dirs, "outs", "filtered_feature_bc_matrix"))

matrix_dirs <- matrix_dirs[dir.exists(matrix_dirs)]
matrix_dirs <- matrix_dirs[!grepl("Undetermined", matrix_dirs)]

sample_list <- as.list(matrix_dirs)
print(sample_list)


#sample_list <- c("/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA1140_Ctrl/outs/filtered_feature_bc_matrix", 
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA275_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA278_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA287_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA639_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA719_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CA843_Ctrl/outs/filtered_feature_bc_matrix",
#                 "/share/crsp/lab/dalawson/share/BRCA1_controls/cell_ranger/CB375_Ctrl/outs/filtered_feature_bc_matrix")
sample_seurat_objs <- list()
for (i in seq_along(matrix_dirs)){
  sample_mat <- Read10X(data.dir = matrix_dirs[i])
  #sample_name <- paste0("sample", toString(i))
  sample_name <- basename(dirname(dirname(matrix_dirs[i])))
  
  print(sample_name)
  sample_mat <- CreateSeuratObject(counts = sample_mat, project = sample_name)
  sample_mat$sample <- sample_name
  sample_seurat_objs[[sample_name]] <- sample_mat
}


combined <- Reduce(function(x, y) merge(x, y), sample_seurat_objs)
saveRDS(combined, "/share/crsp/lab/dalawson/share/BRCA1_controls/combined_BRCA.rds")
