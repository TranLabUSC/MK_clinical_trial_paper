##############################################################################
# Consolidated R Script for Heatmaps of Immune Checkpoint Genes
# Ordered by OS.months.
##############################################################################

# -----------------------
# 0. Load Libraries
# -----------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(gridExtra)  # for grid.arrange if desired
library(Matrix)     # for rowMeans on sparse data

# -----------------------
# 1. Define Survivor Groups & Survival Data
# -----------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Adjust to your actual file path:
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

# We only keep "UF" site patients if that is required
survival_df <- survival_df[survival_df$site == "UF", ]

# Create a column "SurvivalGroup"
survival_df$SurvivalGroup <- ifelse(
  survival_df$patient_id %in% short_term_survivor_group, "ShortTerm",
  ifelse(survival_df$patient_id %in% long_term_survivor_group,  "LongTerm", NA)
)

# Filter out patients not in either group
survival_df <- survival_df[!is.na(survival_df$SurvivalGroup), ]

# Make sure 'ShortTerm' is factor level 1
survival_df$SurvivalGroup <- factor(
  survival_df$SurvivalGroup,
  levels = c("ShortTerm", "LongTerm")
)

# We will use the column "OS.months." to order the patients by survival.
# Check if "OS.months." is present:
if(!"OS.months." %in% colnames(survival_df)) {
  stop("ERROR: 'OS.months.' column not found in survival_df.")
}

# -----------------------
# 2. Load Seurat Object & Metadata
# -----------------------
# Adjust to your actual file path:
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

# -----------------------
# 3. Define Cell Type -> Cluster Mapping
# -----------------------
mapping <- list(
  "Proliferating_Effector" = c(14, 16, 17),
  "Effector_CD8" = c(1),
  "Memory_Precursor_Effector_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "CD8" = c(1, 2, 8, 10, 12, 14, 16, 17),
  "All_T_Cells" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 17, 18)
)

celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# -----------------------
# 4. Define Timepoints & Immune Checkpoint Genes
# -----------------------
timepoints <- c("Pre", "C1", "C2", "C4")

# Below is an example set of checkpoint genes by their HGNC symbols
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

# -----------------------
# 5. Helper Function:
#    Compute Average Expression (log-normalized) by Patient
# -----------------------
get_average_expression_matrix <- function(
    seurat_obj,
    metadata,
    genes, 
    cluster_ids,
    cluster_col    = "seurat_clusters",
    timepoint      = "Pre",
    timepoint_col  = "TimePoint",
    patient_col    = "Patient"
) {
  
  # 1) Subset metadata to cells in the desired cluster & timepoint
  meta_filt <- metadata %>%
    filter(.data[[cluster_col]] %in% cluster_ids) %>%
    filter(.data[[timepoint_col]] == timepoint)
  
  if (nrow(meta_filt) == 0) {
    return(NULL)
  }
  
  # 2) Grab the cell barcodes
  cells_use <- rownames(meta_filt)
  
  # 3) Subset expression data from the Seurat "RNA" assay
  #    By default, seurat_obj[["RNA"]]@data is typically log-normalized counts
  expr_mat <- seurat_obj[["RNA"]]@data[genes, cells_use, drop=FALSE]
  
  # 4) Create a list of cells for each patient
  meta_filt$cell_barcode <- rownames(meta_filt)
  patient_list <- split(cells_use, meta_filt[[patient_col]])
  
  # Prepare output matrix: rows=genes, cols=patients
  unique_patients <- names(patient_list)
  out_mat <- matrix(
    NA_real_, 
    nrow = length(genes), 
    ncol = length(unique_patients),
    dimnames = list(genes, unique_patients)
  )
  
  # For each patient, compute the average expression across that patient's cells
  for (pat in unique_patients) {
    cell_subset <- patient_list[[pat]]
    if (length(cell_subset) == 0) next
    if (!all(cell_subset %in% colnames(expr_mat))) next
    
    out_mat[, pat] <- Matrix::rowMeans(expr_mat[, cell_subset, drop=FALSE])
  }
  
  return(out_mat)
}

# -----------------------
# 6. Main Loop: Generate Heatmaps per Cell Type, per Timepoint
# -----------------------
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (celltype in unique(celltype_to_cluster$celltype)) {
  
  # The cluster IDs for the current cell type
  tcell_cluster_list <- as.vector(celltype_to_cluster[
    celltype_to_cluster$celltype == celltype, 
    "cluster"
  ])
  
  # We'll collect a list of matrices (one per timepoint)
  mat_list <- list()
  
  # -----------------------
  # 6A. Build Expression Matrices for Each Timepoint
  # -----------------------
  for (tp in timepoints) {
    avg_expr_mat <- get_average_expression_matrix(
      seurat_obj   = seurat_object_t_cells,
      metadata     = seurat_metadata_t_cells,
      genes        = immune_checkpoint_genes,
      cluster_ids  = tcell_cluster_list,
      cluster_col  = "seurat_clusters",
      timepoint    = tp,
      timepoint_col= "TimePoint",
      patient_col  = "Patient"
    )
    
    if (is.null(avg_expr_mat)) {
      message(paste("No cells found for", celltype, "at timepoint", tp))
      next
    }
    
    # Reorder columns by ascending OS.months.
    patients_in_mat <- colnames(avg_expr_mat)
    surv_sub <- survival_df %>%
      filter(patient_id %in% patients_in_mat) %>%
      arrange(OS.months.)
    
    # Keep only the patients present in avg_expr_mat
    ordered_patients <- intersect(surv_sub$patient_id, patients_in_mat)
    avg_expr_mat <- avg_expr_mat[, ordered_patients, drop=FALSE]
    
    mat_list[[tp]] <- avg_expr_mat
  }
  
  # If we have no data for any timepoint, skip
  if (length(mat_list) == 0) {
    message(paste("Skipping", celltype, "because no data was found."))
    next
  }
  
  # -----------------------
  # 6B. Create pheatmaps for each timepoint and combine
  # -----------------------
  pheatmap_list <- list()
  
  for (tp in timepoints) {
    if (!tp %in% names(mat_list)) {
      # Means we had no data for that timepoint
      pheatmap_list[[tp]] <- NULL
      next
    }
    
    expr_mat_tp <- mat_list[[tp]]
    
    # 1) Identify rows that have zero variance
    row_sd <- apply(expr_mat_tp, 1, sd, na.rm = TRUE)
    zero_var_rows <- which(row_sd == 0)
    
    # 2) Remove those rows
    if (length(zero_var_rows) > 0) {
      expr_mat_tp <- expr_mat_tp[-zero_var_rows, ]
    }
    
    # If you'd like a custom color scale:
    # my_palette <- colorRampPalette(c("blue","white","red"))(100)
    
    # We keep columns unclustered (cluster_cols = FALSE) to preserve survival-based ordering
    p_obj <- pheatmap(
      expr_mat_tp,
      scale = "row",           # often helps see variation across patients
      cluster_rows = FALSE,     # or FALSE if you want your genes in a given order
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = paste(celltype, "-", tp),
      fontsize = 10,
      border_color = NA
      # color = my_palette
    )
    
    pheatmap_list[[tp]] <- p_obj
  }
  
  # -----------------------
  # 6C. Arrange 4 heatmaps (Pre, C1, C2, C4) in a 2x2 grid
  # -----------------------
  pdf_filename <- file.path(output_dir, paste0(celltype, "_ImmuneCheckpoints_Heatmaps.pdf"))
  pdf(pdf_filename, width=18, height=10)
  
  # Extract the gtable from each pheatmap object
  heatmaps_to_plot <- lapply(pheatmap_list, function(x) {
    if (!is.null(x)) x$gtable else NULL
  })
  
  # We'll place them in the layout:
  #   Pre (top-left), C1 (top-right), 
  #   C2 (bottom-left), C4 (bottom-right).
  grid.arrange(
    heatmaps_to_plot[["Pre"]],
    heatmaps_to_plot[["C1"]],
    heatmaps_to_plot[["C2"]],
    heatmaps_to_plot[["C4"]],
    ncol = 2, nrow = 2
  )
  
  dev.off()
  message(paste("Saved heatmap PDF for celltype:", celltype, "to:", pdf_filename))
}

##############################################################################
# End of Consolidated Script
##############################################################################











##############################################################################
# Generate a grid of FeaturePlots (rows = genes, columns = patients),
# ordered by ascending OS.months.
##############################################################################

# -----------------------
# 0. Load Libraries
# -----------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)


short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 1. Load Survival Data & Order Patients
# -----------------------
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

# We only keep "UF" site patients if that is required
survival_df <- survival_df[survival_df$site == "UF", ]

# Create a column "SurvivalGroup"
survival_df$SurvivalGroup <- ifelse(
  survival_df$patient_id %in% short_term_survivor_group, "ShortTerm",
  ifelse(survival_df$patient_id %in% long_term_survivor_group,  "LongTerm", NA)
)

# Filter out patients not in either group
survival_df <- survival_df[!is.na(survival_df$SurvivalGroup), ]

# Check that OS.months. and patient_id exist in your data
if (!"OS.months." %in% colnames(survival_df)) {
  stop("ERROR: 'OS.months.' column not found in survival_df.")
}
if (!"patient_id" %in% colnames(survival_df)) {
  stop("ERROR: 'patient_id' column not found in survival_df.")
}

# Sort patients by ascending OS.months.
survival_df <- survival_df %>%
  arrange(OS.months.)

# Create a vector of patient IDs in ascending survival order
ordered_patients <- survival_df$patient_id

# -----------------------
# 2. Load Seurat Object
# -----------------------
# This should contain your UMAP reduction and all cells you want to plot
# (e.g., T cells, or entire dataset).
# Adjust the file path as appropriate:
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# -----------------------
# 3. Reorder the Factor Levels for "Patient" in the Seurat Metadata
#    so that FeaturePlot(..., split.by="Patient") respects your custom order
# -----------------------
# We'll assume the Seurat metadata column is named "Patient"
if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

# Convert to factor and reorder levels
seurat_object_t_cells@meta.data$Patient <- factor(
  seurat_object_t_cells@meta.data$Patient,
  levels = ordered_patients
)

# -----------------------
# 4. Define Immune Checkpoint Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

# -----------------------
# 4. Build a Top Row of Patient Labels
#    (One blank cell + one label for each patient)
# -----------------------
# This ensures columns are labeled only once.

top_labels <- c("", ordered_patients)  # first blank, then all patients

top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() +
    theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=8)
})

# One row, (P+1) columns
top_row <- plot_grid(plotlist = top_label_plots, nrow = 1)

# -----------------------
# 5. Create a Row of Plots for Each Gene
# -----------------------
# For each gene:
#   - We call FeaturePlot(..., split.by="Patient", combine=FALSE)
#   - That returns a list of subplots (one per patient)
#   - We remove all axis/legend/title from each subplot
#   - We apply coord_flip()
#   - We combine them horizontally
#   - Then we prepend a "gene label" column on the left

all_rows_list <- list()

for (gene in immune_checkpoint_genes) {
  
  # (A) Get subplots for this gene across all patients
  p_list <- FeaturePlot(
    object   = seurat_object_t_cells,
    features = gene,
    split.by = "Patient",
    cols     = c("lightgrey", "red"),
    combine  = FALSE,
    raster = FALSE,
    pt.size = 0.6,      # Slightly bigger points to highlight red cells
    order = TRUE,        # Ensures red (high expression) points are plotted last
    min.cutoff = 'q10',
    max.cutoff = 'q90'        # Ensures red (high expression) points are plotted last
  )
  # p_list should have length = #patients
  
  # (B) Remove all labels, legends, titles, add coord_flip()
  p_list <- lapply(p_list, function(x) {
    x + 
      coord_flip() +
      ggtitle(NULL) +
      NoLegend() +
      theme_void() # removes axis text/ticks/lines
  })
  
  # (C) Combine these subplots horizontally (i.e., 1 row, #patients columns)
  gene_row_plots <- plot_grid(plotlist = p_list, nrow = 1)
  
  # (D) Create a "gene label" plot (first column)
  gene_label_plot <- ggplot() +
    theme_void() +
    # rotate text 90 degrees for vertical label on the left
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=8)
  
  # (E) Combine the gene label (first column) + the splitted subplots
  #     total columns = 1 + #patients
  row_for_gene <- plot_grid(
    gene_label_plot,
    gene_row_plots,
    nrow       = 1,
    rel_widths = c(0.05, 1)  # make label column narrower
  )
  
  all_rows_list[[gene]] <- row_for_gene
}

# -----------------------
# 6. Stack All Gene-Rows Vertically
#    Then place the top_row of patient labels above them
# -----------------------
main_matrix <- plot_grid(plotlist = all_rows_list, ncol = 1)

final_plot <- plot_grid(
  top_row,
  main_matrix,
  ncol = 1,
  rel_heights = c(0.1, 1)  # adjust as desired
)

# -----------------------
# 6. Combine the Plots into One Grid & Save as PDF
# -----------------------
output_pdf <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoints_FeaturePlots.pdf"
pdf(output_pdf, width = 4 * length(ordered_patients), height = 3 * length(immune_checkpoint_genes))

# We'll rely on patchwork to arrange the entire list in a single layout.
# ncol = length(ordered_patients) so that each row = 1 gene.
print(final_plot)

dev.off()

message("Saved FeaturePlot grid PDF to: ", output_pdf)

##############################################################################
# End of Script
##############################################################################


# 
# 
# 
# 
# 
# 
# # Load necessary libraries
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(cowplot)
# library(Matrix)  # For handling sparse matrices
# 
# # Define survivor groups
# short_term_survivor_group <- c(7, 10, 12, 14, 18)
# long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)
# 
# # -----------------------
# # 1. Load Survival Data & Order Patients
# # -----------------------
# survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
# survival_df <- read.csv(survival_data)
# 
# # Filter for "UF" site
# survival_df <- survival_df[survival_df$site == "UF", ]
# 
# # Create "SurvivalGroup" column
# survival_df$SurvivalGroup <- ifelse(
#   survival_df$patient_id %in% short_term_survivor_group, "ShortTerm",
#   ifelse(survival_df$patient_id %in% long_term_survivor_group, "LongTerm", NA)
# )
# 
# # Filter out patients not in either group
# survival_df <- survival_df[!is.na(survival_df$SurvivalGroup), ]
# 
# # Check essential columns
# if (!"OS.months." %in% colnames(survival_df)) {
#   stop("ERROR: 'OS.months.' column not found in survival_df.")
# }
# if (!"patient_id" %in% colnames(survival_df)) {
#   stop("ERROR: 'patient_id' column not found in survival_df.")
# }
# 
# # Sort by survival time
# survival_df <- survival_df %>%
#   arrange(OS.months.)
# 
# # Ordered patients
# ordered_patients <- survival_df$patient_id
# 
# # -----------------------
# # 2. Load Seurat Object
# # -----------------------
# seurat_object_t_cells <- readRDS(
#   "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
# )
# 
# # -----------------------
# # 3. Merge SurvivalGroup into Seurat Metadata
# # -----------------------
# # Ensure 'Patient' column exists
# if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
#   stop("ERROR: 'Patient' column not found in Seurat object metadata.")
# }
# 
# # Ensure both keys are character
# seurat_object_t_cells@meta.data$Patient <- as.character(seurat_object_t_cells@meta.data$Patient)
# survival_df$patient_id <- as.character(survival_df$patient_id)
# 
# # Perform the left join
# seurat_object_t_cells@meta.data <- seurat_object_t_cells@meta.data %>%
#   left_join(survival_df[, c("patient_id", "SurvivalGroup")], by = c("Patient" = "patient_id"))
# 
# # Verify merge
# if (!"SurvivalGroup" %in% colnames(seurat_object_t_cells@meta.data)) {
#   stop("ERROR: 'SurvivalGroup' column not found in Seurat object metadata after merging.")
# }
# 
# # -----------------------
# # 4. Reorder the Factor Levels for "Patient" in the Seurat Metadata
# # -----------------------
# seurat_object_t_cells@meta.data$Patient <- factor(
#   seurat_object_t_cells@meta.data$Patient,
#   levels = ordered_patients
# )
# 
# # -----------------------
# # 5. Define Immune Checkpoint Genes
# # -----------------------
# immune_checkpoint_genes <- c(
#   "PDCD1",    # PD-1
#   "CD274",    # PD-L1
#   "PDCD1LG2", # PD-L2
#   "CTLA4",
#   "LAG3",
#   "HAVCR2",   # TIM-3
#   "TIGIT",
#   "BTLA",
#   "VSIR",     # VISTA
#   "CD276",    # B7-H3
#   "VTCN1",    # B7-H4
#   "IDO1",
#   "IDO2",
#   "CD47",
#   "CD28",
#   "ICOS",
#   "TNFRSF4",  # OX40
#   "TNFRSF9",  # 4-1BB
#   "TNFRSF18", # GITR
#   "NCR3LG1",  # B7-H6
#   "HHLA2"     # B7-H7
# )
# 
# # -----------------------
# # 6. Define the Enhanced Plotting Function
# # -----------------------
# create_plot_grid <- function(seurat_obj, genes, split_by, group_levels, top_labels, output_pdf, plot_width, plot_height) {
#   
#   # Check if split_by column exists
#   if (!split_by %in% colnames(seurat_obj@meta.data)) {
#     stop(paste("ERROR: '", split_by, "' column not found in Seurat object metadata.", sep = ""))
#   }
#   
#   # Remove cells with NA in split_by column
#   initial_cell_count <- ncol(seurat_obj)
#   seurat_obj <- subset(seurat_obj, !is.na(seurat_obj@meta.data[[split_by]]))
#   final_cell_count <- ncol(seurat_obj)
#   
#   if (final_cell_count < initial_cell_count) {
#     warning(initial_cell_count - final_cell_count, " cells were removed due to NA in '", split_by, "'.")
#   }
#   
#   # Check if all genes are present
#   all_genes <- rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
#   
#   missing_genes <- setdiff(genes, all_genes)
#   if (length(missing_genes) > 0) {
#     warning("The following genes are not present in the Seurat object and will be skipped: ", 
#             paste(missing_genes, collapse = ", "))
#     genes <- setdiff(genes, missing_genes)
#   }
#   
#   # Check for zero-expression genes
#   if (length(genes) > 0) {
#     zero_expression_genes <- genes[Matrix::rowSums(seurat_obj@assays[[DefaultAssay(seurat_obj)]]@data[genes, ]) == 0]
#     if (length(zero_expression_genes) > 0) {
#       warning("The following genes have zero expression across all cells and will be skipped: ", 
#               paste(zero_expression_genes, collapse = ", "))
#       genes <- setdiff(genes, zero_expression_genes)
#     }
#   }
#   
#   if (length(genes) == 0) {
#     stop("No genes left to plot after filtering missing and zero-expression genes.")
#   }
#   
#   # Create top labels
#   label_plots <- lapply(top_labels, function(txt) {
#     ggplot() +
#       theme_void() +
#       geom_text(aes(x=0.5, y=0.5, label=txt), size=8)
#   })
#   
#   top_row <- plot_grid(plotlist = label_plots, nrow = 1)
#   
#   all_rows_list <- list()
#   
#   for (gene in genes) {
#     
#     # Attempt to get FeaturePlot subplots
#     p_list <- tryCatch({
#       FeaturePlot(
#         object   = seurat_obj,
#         features = gene,
#         split.by = split_by,
#         cols     = c("lightgrey", "red"),
#         combine  = FALSE
#       )
#     }, error = function(e) {
#       warning("Error generating FeaturePlot for gene ", gene, ": ", e$message)
#       return(NULL)
#     })
#     
#     if (is.null(p_list)) {
#       next  # Skip to next gene if FeaturePlot failed
#     }
#     
#     # Check if number of plots matches group levels
#     if (length(p_list) != length(group_levels)) {
#       warning(paste("Number of plots for gene", gene, "does not match number of groups. Skipping this gene."))
#       next
#     }
#     
#     # Modify plots: remove labels, legends, titles, add coord_flip
#     p_list <- lapply(p_list, function(x) {
#       x + 
#         coord_flip() +
#         ggtitle(NULL) +
#         NoLegend() +
#         theme_void()
#     })
#     
#     # Combine gene plots horizontally
#     gene_row_plots <- plot_grid(plotlist = p_list, nrow = 1)
#     
#     # Create gene label plot
#     gene_label_plot <- ggplot() +
#       theme_void() +
#       geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=8)
#     
#     # Combine gene label and gene plots
#     row_for_gene <- plot_grid(
#       gene_label_plot,
#       gene_row_plots,
#       nrow = 1,
#       rel_widths = c(0.05, 1)
#     )
#     
#     all_rows_list[[gene]] <- row_for_gene
#   }
#   
#   if (length(all_rows_list) == 0) {
#     stop("No genes left to plot after processing.")
#   }
#   
#   # Stack all gene rows vertically
#   main_matrix <- plot_grid(plotlist = all_rows_list, ncol = 1)
#   
#   # Combine top row and main matrix
#   final_plot <- plot_grid(
#     top_row,
#     main_matrix,
#     ncol = 1,
#     rel_heights = c(0.1, 1)
#   )
#   
#   # Save to PDF
#   pdf(output_pdf, width = plot_width, height = plot_height)
#   print(final_plot)
#   dev.off()
#   
#   message("Saved FeaturePlot grid PDF to: ", output_pdf)
# }
# 
# # -----------------------
# # 7. Check for Missing Genes and Metadata Integrity
# # -----------------------
# # Check which genes are missing and handle accordingly
# all_genes <- rownames(seurat_object_t_cells[[DefaultAssay(seurat_object_t_cells)]])
# 
# missing_genes <- setdiff(immune_checkpoint_genes, all_genes)
# if (length(missing_genes) > 0) {
#   warning("The following genes are not present in the Seurat object and will be skipped: ", 
#           paste(missing_genes, collapse = ", "))
#   immune_checkpoint_genes <- setdiff(immune_checkpoint_genes, missing_genes)
# }
# 
# # Check for zero-expression genes
# zero_expression_genes <- immune_checkpoint_genes[Matrix::rowSums(seurat_object_t_cells@assays[[DefaultAssay(seurat_object_t_cells)]]@data[immune_checkpoint_genes, ]) == 0]
# if (length(zero_expression_genes) > 0) {
#   warning("The following genes have zero expression across all cells and will be skipped: ", 
#           paste(zero_expression_genes, collapse = ", "))
#   immune_checkpoint_genes <- setdiff(immune_checkpoint_genes, zero_expression_genes)
# }
# 
# # Check for NA in 'Patient' and 'SurvivalGroup'
# na_patients <- sum(is.na(seurat_object_t_cells@meta.data$Patient))
# na_survival <- sum(is.na(seurat_object_t_cells@meta.data$SurvivalGroup))
# 
# cat("Number of cells with NA in 'Patient':", na_patients, "\n")
# cat("Number of cells with NA in 'SurvivalGroup':", na_survival, "\n")
# 
# # Optionally, remove cells with NA
# seurat_object_t_cells <- subset(seurat_object_t_cells, !is.na(Patient) & !is.na(SurvivalGroup))
# 
# # -----------------------
# # 8. Create and Save the Original Patient-Based Plot
# # -----------------------
# output_pdf_patient <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoints_FeaturePlots_Patients.pdf"
# 
# create_plot_grid(
#   seurat_obj   = seurat_object_t_cells,
#   genes        = immune_checkpoint_genes,
#   split_by     = "Patient",
#   group_levels = as.character(ordered_patients),
#   top_labels   = c("", as.character(ordered_patients)),
#   output_pdf   = output_pdf_patient,
#   plot_width   = 4 * length(ordered_patients),  # Adjust as needed
#   plot_height  = 3 * length(immune_checkpoint_genes)  # Adjust as needed
# )
# 
# # -----------------------
# # 9. Create and Save the Survival Group-Based Plot
# # -----------------------
# # Define group levels and labels
# survival_groups <- c("ShortTerm", "LongTerm")
# top_labels_survival <- c("", survival_groups)
# 
# output_pdf_survival <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoints_FeaturePlots_SurvivalGroups.pdf"
# 
# create_plot_grid(
#   seurat_obj   = seurat_object_t_cells,
#   genes        = immune_checkpoint_genes,
#   split_by     = "SurvivalGroup",
#   group_levels = survival_groups,
#   top_labels   = top_labels_survival,
#   output_pdf   = output_pdf_survival,
#   plot_width   = 4 * length(survival_groups),  # Adjust as needed
#   plot_height  = 3 * length(immune_checkpoint_genes)  # Adjust as needed
# )
# 
# ##############################################################################
# # End of Script
# ##############################################################################
# 









# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 1. Define Survivor Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 2. Load and Prepare Survival Data
# -----------------------
# Load survival data
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df$patient_id <- as.integer(survival_df$patient_id)

# Filter for 'UF' site patients and assign SurvivalGroup
survival_df <- survival_df %>%
  filter(site == "UF") %>%
  mutate(SurvivalGroup = case_when(
    patient_id %in% short_term_survivor_group ~ "ShortTerm",
    patient_id %in% long_term_survivor_group ~ "LongTerm",
    patient_id %in% control_group ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(SurvivalGroup))

# Verify necessary columns
required_columns <- c("OS.months.", "patient_id")
missing_columns <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: Missing columns in survival_df:", paste(missing_columns, collapse = ", ")))
}

# Sort patients by ascending OS.months.
survival_df <- survival_df %>%
  arrange(OS.months.)

# -----------------------
# 3. Load and Subset Seurat Object
# -----------------------
# Load Seurat Object
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

# Verify 'Patient' column exists in metadata
if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

# Check overlapping patients
unique_patients_seurat <- unique(seurat_object_t_cells@meta.data$Patient)
unique_patients_survival <- unique(survival_df$patient_id)
common_patients <- intersect(unique_patients_seurat, unique_patients_survival)

if (length(common_patients) == 0) {
  stop("ERROR: No overlapping patients found between Seurat object and survival data.")
} else {
  message(paste("Number of overlapping patients:", length(common_patients)))
}

# Subset Seurat object to include only cells from patients in the survivor groups
patients_to_keep <- survival_df$patient_id
seurat_subset <- subset(seurat_object_t_cells, 
                        cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep))
message(paste("Number of patients after subsetting:", length(unique(seurat_subset@meta.data$Patient))))

# -----------------------
# 4. Add SurvivalGroup Metadata
# -----------------------
# Preserve cell barcodes by converting row names to a column
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  mutate(CellBarcode = rownames(seurat_subset@meta.data))

# Create a mapping dataframe
mapping_df <- survival_df %>%
  select(patient_id, SurvivalGroup) %>%
  rename(Patient = patient_id)

# Join the mapping to the Seurat metadata using 'Patient'
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  left_join(mapping_df, by = "Patient")

# Restore row names from 'CellBarcode' column
rownames(seurat_subset@meta.data) <- seurat_subset@meta.data$CellBarcode

# Remove the temporary 'CellBarcode' column
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  select(-CellBarcode)

# Verify the new metadata
print(table(seurat_subset@meta.data$SurvivalGroup, useNA = "ifany"))

# Check for any NA values in SurvivalGroup
num_na <- sum(is.na(seurat_subset@meta.data$SurvivalGroup))
if (num_na > 0) {
  warning(paste(num_na, "cells have NA for SurvivalGroup."))
}

# Set SurvivalGroup as a factor with desired levels
seurat_subset@meta.data$SurvivalGroup <- factor(
  seurat_subset@meta.data$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)
# -----------------------
# 5. Downsample the Larger Group to Match the Smaller Group
# -----------------------

# Set seed for reproducibility
set.seed(123)

# Get cell names for each SurvivalGroup
control_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == "Control")
short_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == "ShortTerm")
long_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == "LongTerm")

# Determine the number of cells in each group
num_control <- length(control_cells)
num_short <- length(short_cells)
num_long <- length(long_cells)

cat("Number of Control cells:", num_control, "\n")
cat("Number of ShortTerm cells:", num_short, "\n")
cat("Number of LongTerm cells:", num_long, "\n")

# Determine the minimum number of cells
min_cells <- min(num_control, num_short, num_long)
cat("Downsampling to", min_cells, "cells per group.\n")

# Randomly sample cells from each group
control_cells_downsampled <- sample(control_cells, min_cells)
short_cells_downsampled <- sample(short_cells, min_cells)
long_cells_downsampled <- sample(long_cells, min_cells)

# Combine the downsampled cells
downsampled_cells <- c(control_cells_downsampled, short_cells_downsampled, long_cells_downsampled)

# Create the downsampled Seurat object
seurat_downsampled <- subset(seurat_subset, cells = downsampled_cells)


# Verify the new cell counts
cell_counts <- table(seurat_downsampled@meta.data$SurvivalGroup)
print(cell_counts)

# -----------------------
# 6. Define Immune Checkpoint Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

# -----------------------
# 7. Build a Top Row of Survival Group Labels
#    (One blank cell + one label for each group)
# -----------------------
# This ensures columns are labeled only once.

# Adjust top_labels to match the number of groups
top_labels <- c("", levels(seurat_downsampled@meta.data$SurvivalGroup))  # first blank, then all groups

top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() +
    theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=8)
})

# One row, (Group + 1) columns
top_row <- plot_grid(plotlist = top_label_plots, nrow = 1)

# -----------------------
# 8. Create a Row of Plots for Each Gene
# -----------------------
# For each gene:
#   - Call FeaturePlot(..., split.by="SurvivalGroup", combine=FALSE)
#   - Remove axis/legend/title from each subplot
#   - Apply coord_flip()
#   - Combine horizontally
#   - Prepend a "gene label" column on the left

all_rows_list <- list()

for (gene in immune_checkpoint_genes) {
  
  # (A) Get subplots for this gene across all groups
  p_list <- FeaturePlot(
    object   = seurat_downsampled,
    features = gene,
    split.by = "SurvivalGroup",
    cols     = c("lightgrey", "red"),
    combine  = FALSE,
    raster = FALSE,
    pt.size = 0.8,      # Slightly bigger points to highlight red cells
    order = TRUE,        # Ensures red (high expression) points are plotted last
    min.cutoff = 'q10',
    max.cutoff = 'q90'
  )
  # p_list should have length = number of groups (2)
  
  # (B) Remove all labels, legends, titles, add coord_flip()
  p_list <- lapply(p_list, function(x) {
    x + 
      coord_flip() +
      ggtitle(NULL) +
      NoLegend() +
      theme_void() # removes axis text/ticks/lines
  })
  
  # (C) Combine these subplots horizontally (i.e., 1 row, #groups columns)
  gene_row_plots <- plot_grid(plotlist = p_list, nrow = 1)
  
  # (D) Create a "gene label" plot (first column)
  gene_label_plot <- ggplot() +
    theme_void() +
    # rotate text 90 degrees for vertical label on the left
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=8)
  
  # (E) Combine the gene label (first column) + the splitted subplots
  #     total columns = 1 + #groups
  row_for_gene <- plot_grid(
    gene_label_plot,
    gene_row_plots,
    nrow       = 1,
    rel_widths = c(0.05, 1)  # make label column narrower
  )
  
  all_rows_list[[gene]] <- row_for_gene
}

# -----------------------
# 9. Stack All Gene-Rows Vertically
#    Then place the top_row of group labels above them
# -----------------------
main_matrix <- plot_grid(plotlist = all_rows_list, ncol = 1)

final_plot <- plot_grid(
  top_row,
  main_matrix,
  ncol = 1,
  rel_heights = c(0.05, 1)  # adjust as desired
)

# -----------------------
# 10. Combine the Plots into One Grid & Save as PDF
# -----------------------
output_pdf <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoints_FeaturePlots_SurvivalGroup_with_control.pdf"

# Determine the number of plots to set appropriate PDF size
num_genes <- length(immune_checkpoint_genes)
# Estimate width: (number of groups + gene label) * per group width
# height: number of genes * per gene height

# Adjust these values based on your specific needs and the number of groups
pdf_width <- 4 * 3  # e.g., 8 inches
pdf_height <- 3 * num_genes  # e.g., 60 inches for 20 genes; adjust as needed

# Open PDF device
pdf(output_pdf, width = pdf_width, height = pdf_height)

# Print the final plot
print(final_plot)

# Close the PDF device
dev.off()

message("Saved FeaturePlot grid PDF to: ", output_pdf)

##############################################################################
# End of Script
##############################################################################







##############################################################################
# this one is for a new kind of visualization
##############################################################################
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 1. Define Survivor Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 2. Load and Prepare Survival Data
# -----------------------
# Load survival data
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df$patient_id <- as.integer(survival_df$patient_id)

# Filter for 'UF' site patients and assign SurvivalGroup
survival_df <- survival_df %>%
  filter(site == "UF") %>%
  mutate(SurvivalGroup = case_when(
    patient_id %in% short_term_survivor_group ~ "ShortTerm",
    patient_id %in% long_term_survivor_group ~ "LongTerm",
    patient_id %in% control_group ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(SurvivalGroup))

# Verify necessary columns
required_columns <- c("OS.months.", "patient_id")
missing_columns <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: Missing columns in survival_df:", paste(missing_columns, collapse = ", ")))
}

# Sort patients by ascending OS.months.
survival_df <- survival_df %>%
  arrange(OS.months.)

# -----------------------
# 3. Load and Subset Seurat Object
# -----------------------
# Load Seurat Object
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

# Verify 'Patient' column exists in metadata
if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

# Check overlapping patients
unique_patients_seurat <- unique(seurat_object_t_cells@meta.data$Patient)
unique_patients_survival <- unique(survival_df$patient_id)
common_patients <- intersect(unique_patients_seurat, unique_patients_survival)

if (length(common_patients) == 0) {
  stop("ERROR: No overlapping patients found between Seurat object and survival data.")
} else {
  message(paste("Number of overlapping patients:", length(common_patients)))
}

# Subset Seurat object to include only cells from patients in the survivor groups
patients_to_keep <- survival_df$patient_id
seurat_subset <- subset(
  seurat_object_t_cells, 
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)
message(paste("Number of patients after subsetting:", length(unique(seurat_subset@meta.data$Patient))))

# -----------------------
# 4. Add SurvivalGroup Metadata
# -----------------------
# Preserve cell barcodes by converting row names to a column
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  mutate(CellBarcode = rownames(seurat_subset@meta.data))

# Create a mapping dataframe
mapping_df <- survival_df %>%
  select(patient_id, SurvivalGroup) %>%
  rename(Patient = patient_id)

# Join the mapping to the Seurat metadata using 'Patient'
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  left_join(mapping_df, by = "Patient")

# Restore row names from 'CellBarcode' column
rownames(seurat_subset@meta.data) <- seurat_subset@meta.data$CellBarcode

# Remove the temporary 'CellBarcode' column
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  select(-CellBarcode)

# Verify the new metadata
print(table(seurat_subset@meta.data$SurvivalGroup, useNA = "ifany"))

# Check for any NA values in SurvivalGroup
num_na <- sum(is.na(seurat_subset@meta.data$SurvivalGroup))
if (num_na > 0) {
  warning(paste(num_na, "cells have NA for SurvivalGroup."))
}

# Set SurvivalGroup as a factor with desired levels
seurat_subset@meta.data$SurvivalGroup <- factor(
  seurat_subset@meta.data$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------
# 5. (Optionally) Subset to T-cells only, or any other pre-processing
#    Already done above by removing cluster 13, etc.
# -----------------------
# (Nothing additional here, unless desired.)

# -----------------------
# 6. **NEW**: Subset to TimePoints of Interest (Pre, C1, C2)
# -----------------------
valid_timepoints <- c("Pre","C1","C2")
seurat_subset <- subset(
  seurat_subset, 
  subset = TimePoint %in% valid_timepoints
)
message("Cells retained for timepoints Pre, C1, C2: ", ncol(seurat_subset))

# -----------------------
# 7. **NEW**: Downsample *within each SurvivalGroup* so that
#    each of the 3 timepoints in that group has the same number of cells.
# -----------------------

set.seed(123)  # for reproducibility

final_cells <- c()

# Loop over each group
for(grp in levels(seurat_subset@meta.data$SurvivalGroup)) {
  
  # Within this group, gather cells for each timepoint
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  # For safety, check only the timepoints we want
  pre_cells <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "Pre"))
  c1_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C1"))
  c2_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C2"))
  
  n_pre <- length(pre_cells)
  n_c1  <- length(c1_cells)
  n_c2  <- length(c2_cells)
  
  # Find the smallest among the three timepoints in this group
  n_min <- min(n_pre, n_c1, n_c2)
  
  if(n_min == 0) {
    # If one timepoint has zero cells, this group effectively can't be balanced.
    # You can decide to skip or keep as is. Here, we skip the group if any timepoint is zero.
    warning(paste("Skipping group", grp, "because at least one timepoint has 0 cells."))
    next
  }
  
  # Randomly downsample each timepoint to match n_min
  pre_cells_sub <- sample(pre_cells, n_min)
  c1_cells_sub  <- sample(c1_cells, n_min)
  c2_cells_sub  <- sample(c2_cells, n_min)
  
  # Combine them back for this group
  final_cells <- c(final_cells, pre_cells_sub, c1_cells_sub, c2_cells_sub)
}

# Now subset to only the final downsampled cells
seurat_timepoint_downsampled <- subset(seurat_subset, cells = final_cells)

cat("Final cell counts after timepoint-downsampling within each group:\n")
table(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint
)

# -----------------------
# 8. **NEW**: Create a combined "Group_TimePoint" metadata column
#    so that we can do a single FeaturePlot call with split.by = "Group_TimePoint".
# -----------------------
seurat_timepoint_downsampled@meta.data$Group_TimePoint <- paste(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint,
  sep = "_"
)

# Define factor levels so subplots appear in a nice order:
desired_order <- c(
  "Control_Pre", "Control_C1", "Control_C2",
  "ShortTerm_Pre", "ShortTerm_C1", "ShortTerm_C2",
  "LongTerm_Pre", "LongTerm_C1", "LongTerm_C2"
)
seurat_timepoint_downsampled@meta.data$Group_TimePoint <- factor(
  seurat_timepoint_downsampled@meta.data$Group_TimePoint,
  levels = desired_order
)

# -----------------------
# 9. Define Immune Checkpoint Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

# -----------------------
# 10. Create a Top Row of Labels for the 9 Columns
# -----------------------
top_labels <- levels(seurat_timepoint_downsampled@meta.data$Group_TimePoint)
# Or you could make them more human-readable:
# top_labels <- c("Control\nPre","Control\nC1",...)

# Each label as its own small plot
top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() +
    theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=4)
})

# Put them all in one row
top_row <- plot_grid(plotlist = top_label_plots, nrow = 1)

# -----------------------
# 11. Build the main FeaturePlot rows
# -----------------------
all_rows_list <- list()

for (gene in immune_checkpoint_genes) {
  
  # (A) Run FeaturePlot with split.by = "Group_TimePoint" so we get 9 subplots for each gene
  p_list <- FeaturePlot(
    object   = seurat_timepoint_downsampled,
    features = gene,
    split.by = "Group_TimePoint",
    cols     = c("lightgrey", "red"),
    combine  = FALSE,
    pt.size  = 0.6,
    order    = TRUE,     # to ensure higher-expressing cells are on top
    min.cutoff = "q10",
    max.cutoff = "q90"
  )
  # p_list should be a list of 9 subplots in the order of levels(Group_TimePoint)
  
  # (B) Remove legends, axis ticks, etc., from each subplot
  p_list <- lapply(p_list, function(x) {
    x + 
      coord_flip() +
      ggtitle(NULL) +
      NoLegend() +
      theme_void()  # removes axis text/ticks/lines
  })
  
  # (C) Combine these 9 subplots horizontally
  gene_row_plots <- plot_grid(plotlist = p_list, nrow = 1)
  
  # (D) Create a "gene label" plot (the first column)
  gene_label_plot <- ggplot() +
    theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=5)
  
  # (E) Combine the gene label (first col) + the 9 subplots (next cols)
  row_for_gene <- plot_grid(
    gene_label_plot,
    gene_row_plots,
    nrow       = 1,
    rel_widths = c(0.04, 1)  # narrower for the label
  )
  
  all_rows_list[[gene]] <- row_for_gene
}

# -----------------------
# 12. Stack All Gene-Rows Vertically and put the top_row above
# -----------------------
main_matrix <- plot_grid(plotlist = all_rows_list, ncol = 1)

final_plot <- plot_grid(
  top_row,
  main_matrix,
  ncol = 1,
  rel_heights = c(0.05, 1)  # adjust as needed
)

# -----------------------
# 13. Save as PDF
# -----------------------
output_pdf <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_splitted_featureplots.pdf"

num_genes <- length(immune_checkpoint_genes)
# We have 9 columns of subplots + 1 label column = 10 columns total visually
# Adjust width/height to suit your preferences

pdf_width <- 2 * 10  # e.g., 20 inches wide
pdf_height <- 2 * num_genes  # e.g., 2 inches per gene row
pdf(output_pdf, width = pdf_width, height = pdf_height)
print(final_plot)
dev.off()

message("Saved timepoint-splitted FeaturePlot PDF to: ", output_pdf)










##########################################################################################
# with normalized expression
##########################################################################################
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 1. Define Survivor Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 2. Load and Prepare Survival Data
# -----------------------
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df$patient_id <- as.integer(survival_df$patient_id)

survival_df <- survival_df %>%
  filter(site == "UF") %>%
  mutate(SurvivalGroup = case_when(
    patient_id %in% short_term_survivor_group ~ "ShortTerm",
    patient_id %in% long_term_survivor_group ~ "LongTerm",
    patient_id %in% control_group ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(SurvivalGroup))

required_columns <- c("OS.months.", "patient_id")
missing_columns <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: Missing columns in survival_df:", paste(missing_columns, collapse = ", ")))
}

survival_df <- survival_df %>%
  arrange(OS.months.)

# -----------------------
# 3. Load and Subset Seurat Object
# -----------------------
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# Example: remove cluster 13 if desired
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

common_patients <- intersect(
  unique(seurat_object_t_cells@meta.data$Patient),
  unique(survival_df$patient_id)
)
if (length(common_patients) == 0) {
  stop("ERROR: No overlapping patients found between Seurat object and survival data.")
} else {
  message(paste("Number of overlapping patients:", length(common_patients)))
}

patients_to_keep <- survival_df$patient_id
seurat_subset <- subset(
  seurat_object_t_cells, 
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)

message(paste("Number of patients after subsetting:", length(unique(seurat_subset@meta.data$Patient))))

# -----------------------
# 4. Add SurvivalGroup Metadata
# -----------------------
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  mutate(CellBarcode = rownames(seurat_subset@meta.data))

mapping_df <- survival_df %>%
  select(patient_id, SurvivalGroup) %>%
  rename(Patient = patient_id)

seurat_subset@meta.data <- seurat_subset@meta.data %>%
  left_join(mapping_df, by = "Patient")

rownames(seurat_subset@meta.data) <- seurat_subset@meta.data$CellBarcode

seurat_subset@meta.data <- seurat_subset@meta.data %>%
  select(-CellBarcode)

# Ensure factor levels
seurat_subset@meta.data$SurvivalGroup <- factor(
  seurat_subset@meta.data$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------
# 5. Subset to TimePoints of Interest and Downsample *within group*
#    so each timepoint has the same number of cells
# -----------------------
valid_timepoints <- c("Pre","C1","C2")
seurat_subset <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)

set.seed(123)
final_cells <- c()

for(grp in levels(seurat_subset@meta.data$SurvivalGroup)) {
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  pre_cells <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "Pre"))
  c1_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C1"))
  c2_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C2"))
  
  n_pre <- length(pre_cells)
  n_c1  <- length(c1_cells)
  n_c2  <- length(c2_cells)
  
  n_min <- min(n_pre, n_c1, n_c2)
  if(n_min == 0) {
    warning(paste("Skipping group", grp, " - zero cells in at least one timepoint."))
    next
  }
  
  pre_sub <- sample(pre_cells, n_min)
  c1_sub  <- sample(c1_cells, n_min)
  c2_sub  <- sample(c2_cells, n_min)
  
  final_cells <- c(final_cells, pre_sub, c1_sub, c2_sub)
}

seurat_timepoint_downsampled <- subset(seurat_subset, cells = final_cells)

cat("Final cell counts after timepoint-downsampling:\n")
print(table(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint
))

# -----------------------
# 6. Create "Group_TimePoint" metadata
# -----------------------
seurat_timepoint_downsampled@meta.data$Group_TimePoint <- paste(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint,
  sep = "_"
)

desired_order <- c(
  "Control_Pre", "Control_C1", "Control_C2",
  "ShortTerm_Pre", "ShortTerm_C1", "ShortTerm_C2",
  "LongTerm_Pre", "LongTerm_C1", "LongTerm_C2"
)
seurat_timepoint_downsampled@meta.data$Group_TimePoint <- factor(
  seurat_timepoint_downsampled@meta.data$Group_TimePoint,
  levels = desired_order
)

# -----------------------
# 7. Define Immune Checkpoint Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    
  "CD274",    
  "PDCD1LG2", 
  "CTLA4",
  "LAG3",
  "HAVCR2",   
  "TIGIT",
  "BTLA",
  "VSIR",     
  "CD276",    
  "VTCN1",    
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  
  "TNFRSF9",  
  "TNFRSF18", 
  "NCR3LG1",  
  "HHLA2"     
)

# -----------------------
# 8A. Make the "Original" FeaturePlot (unmodified expression)
# -----------------------
top_labels <- levels(seurat_timepoint_downsampled@meta.data$Group_TimePoint)
top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=4)
})
top_row <- plot_grid(plotlist = top_label_plots, nrow = 1)

all_rows_list <- list()
for (gene in immune_checkpoint_genes) {
  cat("Default Assay is:", DefaultAssay(seurat_timepoint_downsampled), "\n")
  p_list <- FeaturePlot(
    object   = seurat_timepoint_downsampled,
    features = gene,
    split.by = "Group_TimePoint",
    cols     = c("lightgrey", "red"),
    combine  = FALSE,
    pt.size  = 0.6,
    order    = TRUE,
    min.cutoff = "q10",
    max.cutoff = "q90"
  )
  p_list <- lapply(p_list, function(x) {
    x + coord_flip() + ggtitle(NULL) + NoLegend() + theme_void()
  })
  
  gene_row_plots <- plot_grid(plotlist = p_list, nrow = 1)
  
  gene_label_plot <- ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=5)
  
  row_for_gene <- plot_grid(
    gene_label_plot,
    gene_row_plots,
    nrow = 1,
    rel_widths = c(0.04, 1)
  )
  all_rows_list[[gene]] <- row_for_gene
}

main_matrix <- plot_grid(plotlist = all_rows_list, ncol = 1)
final_plot <- plot_grid(top_row, main_matrix, ncol = 1, rel_heights = c(0.05, 1))

output_pdf_original <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_splitted_featureplots.pdf"
num_genes <- length(immune_checkpoint_genes)
pdf_width  <- 2 * 10  
pdf_height <- 2 * num_genes
pdf(output_pdf_original, width = pdf_width, height = pdf_height)
print(final_plot)
dev.off()

message("Saved original FeaturePlot PDF to: ", output_pdf_original)


# -----------------------
# 8B. Make the Z-score "Normalized" FeaturePlot using Pre as baseline (mean & sd)
# -----------------------
# Intersect with the actual genes present in your data:
immune_genes_present <- intersect(immune_checkpoint_genes, rownames(seurat_timepoint_downsampled))

# If for some reason none overlap, handle that scenario:
if (length(immune_genes_present) == 0) {
  stop("No overlap between immune_checkpoint_genes and rownames of object!")
}
# (1) Extract existing log data from the 'RNA' assay:
orig_data <- GetAssayData(seurat_timepoint_downsampled, slot = "counts", assay = "RNA")

subset_data <- orig_data[immune_genes_present, , drop = FALSE]

# (2) We'll create a new matrix, z-scoring each group/timepoint relative to that group's Pre
normalized_data <- subset_data # copy

unique_groups <- levels(seurat_timepoint_downsampled@meta.data$SurvivalGroup)
for (grp in unique_groups) {
  print(grp)
  # Identify Pre cells in this group
  cells_in_pre <- WhichCells(
    seurat_timepoint_downsampled,
    expression = SurvivalGroup == grp & TimePoint == "Pre"
  )
  if (length(cells_in_pre) == 0) next  # skip if no Pre cells for this group
  
  # Mean & SD of each gene among Pre cells
  pre_mean_vec <- Matrix::rowMeans(subset_data[, cells_in_pre, drop = FALSE])
  pre_sd_vec   <- apply(subset_data[, cells_in_pre, drop = FALSE], 1, sd)
  
  # Avoid division by zero
  pre_sd_vec[pre_sd_vec == 0] <- 1
  
  # Identify all cells in this group (Pre, C1, C2)
  cells_in_group <- WhichCells(
    seurat_timepoint_downsampled,
    expression = SurvivalGroup == grp
  )
  if (length(cells_in_group) == 0) next
  
  # Z-score transformation: (value - mean) / sd
  zscore_mat <- sweep(subset_data[, cells_in_group, drop = FALSE], 1, pre_mean_vec, FUN = "-")
  zscore_mat <- sweep(zscore_mat, 1, pre_sd_vec, FUN = "/")
  
  normalized_data[, cells_in_group] <- zscore_mat
}

# (3) Store z-score data in a new assay
normalized_assay <- CreateAssayObject(data = normalized_data)
seurat_timepoint_downsampled[["NormalizedRNA"]] <- normalized_assay

# (4) Plot again using the "NormalizedRNA" assay
all_rows_list_norm <- list()
for (gene in immune_checkpoint_genes) {
  print(gene)
  DefaultAssay(seurat_timepoint_downsampled) <- "NormalizedRNA"
  p_list_norm <- FeaturePlot(
    seurat_timepoint_downsampled,
    features = gene,
    split.by = "Group_TimePoint",
    cols     = c("blue", "white", "red"), # optional palette for negative->positive
    combine  = FALSE,
    pt.size  = 0.6,
    order    = TRUE
    # Because these are z-scores (which can be negative),
    # you might want a numeric cutoff or extended range
    # For example:
    # min.cutoff = -2,
    # max.cutoff = 2
  )
  p_list_norm <- lapply(p_list_norm, function(x) {
    x + coord_flip() + ggtitle(NULL) + NoLegend() + theme_void()
  })
  
  gene_row_plots_norm <- plot_grid(plotlist = p_list_norm, nrow = 1)
  
  gene_label_plot_norm <- ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=5)
  
  row_for_gene_norm <- plot_grid(
    gene_label_plot_norm,
    gene_row_plots_norm,
    nrow = 1,
    rel_widths = c(0.04, 1)
  )
  all_rows_list_norm[[gene]] <- row_for_gene_norm
}

main_matrix_norm <- plot_grid(plotlist = all_rows_list_norm, ncol = 1)
final_plot_norm <- plot_grid(top_row, main_matrix_norm, ncol = 1, rel_heights = c(0.05, 1))

# (5) Save the z-score "normalized" version
output_pdf_norm <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_splitted_featureplots_normalized.pdf"
pdf(output_pdf_norm, width = pdf_width, height = pdf_height)
print(final_plot_norm)
dev.off()

message("Saved z-score normalized FeaturePlot PDF to: ", output_pdf_norm)




















##########################################################################################
# density UMAP for short term and long term survivors
##########################################################################################
###############################################################################
# SCRIPT: 2D Density of Immune Checkpoint Expression in Short vs. Long
#         - All cells in gray background
#         - ShortTerm (blue contours) & LongTerm (red contours) if above threshold
#         - Two variants: no downsampling vs. downsampling
###############################################################################

# -----------------------
# 1. Load Libraries
# -----------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# -----------------------
# 2. Define Short/Long-Term Groups
# -----------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 3. File Paths (Adjust as Needed)
# -----------------------
seurat_obj_path    <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
pdf_no_ds_path     <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_Expression_NoDownsample.pdf"
pdf_ds_path        <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_Expression_Downsample.pdf"

# -----------------------
# 4. Define Immune Checkpoint Genes
#    We'll plot each gene on a separate page
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

# -----------------------
# 5. Load the Seurat Object
# -----------------------
seurat_object <- readRDS(seurat_obj_path)

if (!"Patient" %in% colnames(seurat_object@meta.data)) {
  stop("ERROR: 'Patient' column not found in the Seurat object's metadata.")
}

# Make sure the assay is set so we can fetch gene expression
DefaultAssay(seurat_object) <- "RNA"

# Ensure UMAP is present (or run it if needed)
if (!"umap" %in% Reductions(seurat_object)) {
  seurat_object <- RunUMAP(seurat_object, dims = 1:20)
}

# -----------------------
# 6. Label Which Cells Are Short/Long
#    (But we won't define a "Background" label; all unlabeled are just 'NA' in metadata)
# -----------------------
seurat_object$SurvivalGroup <- NA
seurat_object$SurvivalGroup[
  seurat_object$Patient %in% short_term_survivor_group
] <- "ShortTerm"
seurat_object$SurvivalGroup[
  seurat_object$Patient %in% long_term_survivor_group
] <- "LongTerm"

# -----------------------
# 7. Optional: Downsampling Function
#    Balances short vs. long, keeps all 'other' cells too
# -----------------------
downsample_short_long <- function(sobj) {
  short_cells <- WhichCells(sobj, expression = SurvivalGroup == "ShortTerm")
  long_cells  <- WhichCells(sobj, expression = SurvivalGroup == "LongTerm")
  
  num_short <- length(short_cells)
  num_long  <- length(long_cells)
  
  message("ShortTerm cells: ", num_short, " | LongTerm cells: ", num_long)
  if (num_short == 0 || num_long == 0) {
    warning("One group is zero; cannot downsample. Returning original object.")
    return(sobj)
  }
  
  min_cells <- min(num_short, num_long)
  set.seed(123)
  short_down <- sample(short_cells, min_cells)
  long_down  <- sample(long_cells, min_cells)
  
  keep_cells <- c(short_down, long_down)
  sobj_ds <- subset(sobj, cells = keep_cells)
  
  message("Downsampled to ", ncol(sobj_ds), " total cells.")
  return(sobj_ds)
}

# -----------------------
# 8. 2D Density Plot Function (Expression Threshold)
#    Loops over immune_checkpoint_genes, one page per gene
# -----------------------
create_expression_density_plots <- function(
    sobj,
    gene_list      = immune_checkpoint_genes,
    pdf_file       = "Expression_Density.pdf",
    expression_thr = 1,   # threshold for "high expression"
    flip_coords    = FALSE
) {
  # Extract UMAP for all cells => "all_df"
  umap_coords <- Embeddings(sobj, "umap")
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  all_df <- data.frame(
    cell   = rownames(umap_coords),
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    group  = sobj$SurvivalGroup  # "ShortTerm", "LongTerm", or NA
  )
  
  x_limits <- range(all_df$UMAP_1)
  y_limits <- range(all_df$UMAP_2)
  
  pdf(pdf_file, width = 7, height = 6)
  
  for (gene in gene_list) {
    print(gene)
    
    # 1) Fetch expression for this gene
    expr_values <- FetchData(sobj, vars = gene)[, 1]
    
    # Merge with all_df
    plot_df <- all_df
    plot_df$expression <- expr_values
    
    # 2) Gray background => all cells
    #    We'll just keep that as one data frame
    gray_df <- plot_df
    
    # 3) Subset short-term cells that pass expression threshold
    short_df <- plot_df %>%
      filter(group == "ShortTerm" & expression > expression_thr) %>%
      filter(complete.cases(UMAP_1, UMAP_2))
    
    # 4) Subset long-term cells that pass expression threshold
    long_df <- plot_df %>%
      filter(group == "LongTerm" & expression > expression_thr) %>%
      filter(complete.cases(UMAP_1, UMAP_2))
    
    print(dim(short_df))
    print(dim(long_df))
    
    if(gene == "IDO1") {
      print(short_df)
      print(long_df)
    }
    
    # 5) Build the ggplot
    p <- ggplot() +
      # Gray points for all cells
      geom_point(
        data  = gray_df,
        aes(x = UMAP_1, y = UMAP_2),
        color = "gray",
        size  = 0.5,
        alpha = 0.5
      )
    
    # 5A) Add short-term density if enough data
    if (
      nrow(short_df) >= 2 &&
      length(unique(short_df$UMAP_1)) > 1 &&
      length(unique(short_df$UMAP_2)) > 1
    ) {
      p <- p + geom_density_2d(
        data  = short_df,
        aes(x = UMAP_1, y = UMAP_2),
        color = "blue",
        size  = 0.7
      )
    }
    
    # 5B) Add long-term density if enough data
    if (
      nrow(long_df) >= 2 &&
      length(unique(long_df$UMAP_1)) > 1 &&
      length(unique(long_df$UMAP_2)) > 1
    ) {
      p <- p + geom_density_2d(
        data  = long_df,
        aes(x = UMAP_1, y = UMAP_2),
        color = "red",
        size  = 0.7
      )
    }
    
    # 5C) Coordinate settings
    p <- p +
      coord_fixed() +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      theme_void() +
      theme(
        panel.border     = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_blank(),
        axis.line        = element_line(colour = "black")
      ) +
      ggtitle(paste0("Gene: ", gene, " (Threshold > ", expression_thr, ")"))
    
    if (flip_coords) {
      p <- p + coord_flip()
    }
    
    print(p)
  }
  
  dev.off()
  message("Saved expression-based density plots to: ", pdf_file)
}

# -----------------------
# 9. Make Two Variants
# -----------------------

# (A) No Downsampling
pdf_no_ds <- pdf_no_ds_path
create_expression_density_plots(
  sobj           = seurat_object,
  pdf_file       = pdf_no_ds,
  expression_thr = 0.1,  # adjust the threshold as needed
  flip_coords    = TRUE
)

# (B) Downsampling (Short=Long balanced, keep others)
seurat_ds <- downsample_short_long(seurat_object)
pdf_ds    <- pdf_ds_path
create_expression_density_plots(
  sobj           = seurat_ds,
  pdf_file       = pdf_ds,
  expression_thr = 0.1,
  flip_coords    = TRUE
)

###############################################################################
# End of Script
###############################################################################

























# ################################################################################
# # -----------------------
# # 1. Load Libraries
# # -----------------------
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(cowplot)
# 
# # -----------------------
# # 2. Define Groups
# # -----------------------
# control_group <- c(1, 4, 8, 9, 11)
# short_term_survivor_group <- c(7, 10, 12, 14, 18)
# long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)
# 
# # -----------------------
# # 3. File Paths
# # -----------------------
# seurat_obj_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
# 
# # Output PDFs
# pdf_no_ds_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_timepoint_splitted_NoDownsample.pdf"
# pdf_ds_path    <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_timepoint_splitted_Downsample.pdf"
# 
# # -----------------------
# # 4. Define Immune Checkpoint Genes
# # -----------------------
# immune_checkpoint_genes <- c(
#   "PDCD1",    # PD-1
#   "CD274",    # PD-L1
#   "PDCD1LG2", # PD-L2
#   "CTLA4",
#   "LAG3",
#   "HAVCR2",   # TIM-3
#   "TIGIT",
#   "BTLA",
#   "VSIR",     # VISTA
#   "CD276",    # B7-H3
#   "VTCN1",    # B7-H4
#   "IDO1",
#   "IDO2",
#   "CD47",
#   "CD28",
#   "ICOS",
#   "TNFRSF4",  # OX40
#   "TNFRSF9",  # 4-1BB
#   "TNFRSF18", # GITR
#   "NCR3LG1",  # B7-H6
#   "HHLA2"     # B7-H7
# )
# 
# # -----------------------
# # 5. Load and Subset Seurat Object
# # -----------------------
# seurat_object <- readRDS(seurat_obj_path)
# 
# # Example: remove a cluster if desired
# seurat_object <- subset(seurat_object, subset = seurat_clusters != 13)
# 
# if (!"Patient" %in% colnames(seurat_object@meta.data)) {
#   stop("ERROR: 'Patient' column not found in the Seurat object's metadata.")
# }
# 
# DefaultAssay(seurat_object) <- "RNA"
# 
# if (!"umap" %in% names(seurat_object@reductions)) {
#   seurat_object <- RunUMAP(seurat_object, dims = 1:20)
# }
# 
# # -----------------------
# # 6. Assign SurvivalGroup
# # -----------------------
# seurat_object$SurvivalGroup <- NA
# seurat_object$SurvivalGroup[seurat_object$Patient %in% control_group] <- "Control"
# seurat_object$SurvivalGroup[seurat_object$Patient %in% short_term_survivor_group] <- "ShortTerm"
# seurat_object$SurvivalGroup[seurat_object$Patient %in% long_term_survivor_group]  <- "LongTerm"
# 
# # Factorize in desired order
# seurat_object$SurvivalGroup <- factor(
#   seurat_object$SurvivalGroup,
#   levels = c("Control", "ShortTerm", "LongTerm")
# )
# 
# # -----------------------
# # 7. Subset to TimePoints: Pre, C1, C2
# # -----------------------
# valid_timepoints <- c("Pre", "C1", "C2")
# seurat_tp <- subset(seurat_object, subset = TimePoint %in% valid_timepoints)
# 
# # -----------------------
# # 8. Timepoint-based Downsampling (Optional)
# # -----------------------
# downsample_timepoints <- function(sobj, groups = c("Control","ShortTerm","LongTerm"), tps = c("Pre","C1","C2")) {
#   set.seed(123)
#   final_cells <- c()
#   
#   for (grp in groups) {
#     grp_cells <- WhichCells(sobj, expression = SurvivalGroup == grp)
#     
#     for (tp in tps) {
#       tp_cells <- intersect(grp_cells, WhichCells(sobj, expression = TimePoint == tp))
#       assign(paste0(grp, "_", tp), tp_cells, envir = .GlobalEnv)
#     }
#   }
#   
#   downsampled_cells <- c()
#   for (tp in tps) {
#     grp_counts <- sapply(groups, function(g) length(get(paste0(g, "_", tp), envir = .GlobalEnv)))
#     n_min <- min(grp_counts)
#     if (n_min == 0) {
#       warning(paste("Skipping timepoint:", tp, "- at least one group has zero cells."))
#       next
#     }
#     for (g in groups) {
#       all_cells_g_tp <- get(paste0(g, "_", tp), envir = .GlobalEnv)
#       downsampled_cells <- c(downsampled_cells, sample(all_cells_g_tp, n_min))
#     }
#   }
#   
#   sobj_ds <- subset(sobj, cells = downsampled_cells)
#   return(sobj_ds)
# }
# 
# # -----------------------
# # 9. Create the 2D Density Plot Function
# #    - 3 subplots (Pre, C1, C2) per gene
# #    - Overlays density contours for Control=yellow, ShortTerm=blue, LongTerm=red
# #    - coord_flip() is used
# #    - Capping IDO1 expression at 0.025
# # -----------------------
# plot_2d_density_per_gene <- function(sobj, gene) {
#   timepoints   <- c("Pre", "C1", "C2")
#   group_colors <- c("Control"="yellow", "ShortTerm"="blue", "LongTerm"="red")
#   
#   # Extract UMAP
#   umap_coords <- Embeddings(sobj, "umap")[, 1:2]
#   colnames(umap_coords) <- c("UMAP_1","UMAP_2")
#   
#   # Fetch expression (cap IDO1 if desired)
#   expr_vals <- FetchData(sobj, vars = gene)[,1]
#   # if (gene == "IDO1") {
#   #   expr_vals[expr_vals > 0.025] <- 0.025
#   # }
#   
#   # Build df
#   plot_df <- data.frame(
#     cell          = rownames(umap_coords),
#     UMAP_1        = umap_coords[,1],
#     UMAP_2        = umap_coords[,2],
#     TimePoint     = sobj$TimePoint,
#     SurvivalGroup = sobj$SurvivalGroup,
#     expression    = expr_vals
#   )
#   
#   x_limits <- range(plot_df$UMAP_1)
#   y_limits <- range(plot_df$UMAP_2)
#   
#   p_list <- list()
#   for (tp in timepoints) {
#     df_tp <- dplyr::filter(plot_df, TimePoint == tp)
#     
#     base_plot <- ggplot(df_tp, aes(x=UMAP_1,y=UMAP_2)) +
#       geom_point(color="gray70",size=0.4,alpha=0.5) +
#       coord_fixed() + coord_flip() +  # your requested flip
#       xlim(x_limits) + ylim(y_limits) +
#       theme_void() +
#       ggtitle(paste0(gene," – ",tp))
#     
#     # For each group, draw weighted density contours
#     for (grp in c("Control","ShortTerm","LongTerm")) {
#       df_grp <- df_tp %>% filter(SurvivalGroup == grp & expression > 0)
#       if(nrow(df_grp) > 1) {
#         base_plot <- base_plot + 
#           geom_density_2d(
#             data  = df_grp,
#             aes(x=UMAP_1, y=UMAP_2, weight=expression),
#             color = group_colors[grp],
#             size  = 0.7
#             # h = c(0.3, 0.3)
#             # adjust = 0.3
#           )
#       }
#     }
#     
#     p_list[[tp]] <- base_plot
#   }
#   return(cowplot::plot_grid(plotlist=p_list, nrow=1))
# }
# 
# 
# # Wrapper to loop over genes and save to PDF
# create_2d_density_plots_for_genes <- function(sobj, genes, output_pdf, expression_thr = 0.1) {
#   pdf(output_pdf, width=15, height=5)  # wide enough for 3 subplots side-by-side
#   for (g in genes) {
#     cat("Plotting:", g, "\n")
#     p_gene <- plot_2d_density_per_gene(sobj, gene=g)
#     print(p_gene)
#   }
#   dev.off()
#   message("Saved plots to: ", output_pdf)
# }
# 
# # -----------------------
# # 10A. Generate the Plots Without Downsampling
# # -----------------------
# create_2d_density_plots_for_genes(
#   sobj           = seurat_tp,
#   genes          = immune_checkpoint_genes,
#   output_pdf     = pdf_no_ds_path,
#   expression_thr = 0.1
# )
# 
# # -----------------------
# # 10B. Generate the Plots With Downsampling
# # -----------------------
# seurat_tp_ds <- downsample_timepoints(seurat_tp)
# create_2d_density_plots_for_genes(
#   sobj           = seurat_tp_ds,
#   genes          = immune_checkpoint_genes,
#   output_pdf     = pdf_ds_path,
#   expression_thr = 0.1
# )
# 
# message("Done!")










################################################################################
# -----------------------
# 1. Load Libraries
# -----------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 2. Define Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 3. File Paths
# -----------------------
seurat_obj_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"

# Output PDFs
pdf_no_ds_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_timepoint_splitted_NoDownsample.pdf"
pdf_ds_path    <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_timepoint_splitted_Downsample.pdf"

# -----------------------
# 4. Define Immune Checkpoint Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",   # TIM-3
  "TIGIT",
  "BTLA",
  "VSIR",     # VISTA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",
  "IDO2",
  "CD47",
  "CD28",
  "ICOS",
  "TNFRSF4",  # OX40
  "TNFRSF9",  # 4-1BB
  "TNFRSF18", # GITR
  "NCR3LG1",  # B7-H6
  "HHLA2"     # B7-H7
)

################################################################################
# 1. Load Seurat Object, Assign Groups & Subset
################################################################################
seurat_object <- readRDS(seurat_obj_path)

# Example: remove an unwanted cluster
seurat_object <- subset(seurat_object, subset = seurat_clusters != 13)

# Must have "Patient" in metadata
if (!"Patient" %in% colnames(seurat_object@meta.data)) {
  stop("ERROR: 'Patient' column not found in the metadata.")
}

DefaultAssay(seurat_object) <- "RNA"
if (!"umap" %in% Reductions(seurat_object)) {
  seurat_object <- RunUMAP(seurat_object, dims=1:20)
}

# Mark groups
seurat_object$SurvivalGroup <- NA
seurat_object$SurvivalGroup[ seurat_object$Patient %in% control_group ] <- "Control"
seurat_object$SurvivalGroup[ seurat_object$Patient %in% short_term_survivor_group ] <- "ShortTerm"
seurat_object$SurvivalGroup[ seurat_object$Patient %in% long_term_survivor_group  ] <- "LongTerm"
seurat_object$SurvivalGroup <- factor(
  seurat_object$SurvivalGroup,
  levels = c("Control","ShortTerm","LongTerm")
)

# Keep only timepoints of interest
valid_tps <- c("Pre","C1","C2")
seurat_tp <- subset(seurat_object, subset = TimePoint %in% valid_tps)

################################################################################
# 2. Downsampling Function (Optional)
################################################################################
downsample_timepoints <- function(sobj, groups = c("Control","ShortTerm","LongTerm"), tps=c("Pre","C1","C2")) {
  set.seed(123)
  final_cells <- c()
  for (grp in groups) {
    grp_cells <- WhichCells(sobj, expression = SurvivalGroup == grp)
    for (tp in tps) {
      tp_cells <- intersect(grp_cells, WhichCells(sobj, expression=TimePoint==tp))
      assign(paste0(grp,"_",tp), tp_cells, envir=.GlobalEnv)
    }
  }
  ds_cells <- c()
  for (tp in tps) {
    counts_vec <- sapply(groups, function(g) length(get(paste0(g,"_",tp), envir=.GlobalEnv)))
    n_min <- min(counts_vec)
    if (n_min == 0) {
      warning("Skipping timepoint:", tp, "because at least one group has zero cells.")
      next
    }
    for (g in groups) {
      gp_tp_cells <- get(paste0(g,"_",tp), envir=.GlobalEnv)
      ds_cells <- c(ds_cells, sample(gp_tp_cells, n_min))
    }
  }
  return(subset(sobj, cells = ds_cells))
}

################################################################################
# 3. Subset-Specific q10/q90 Clamping
#    exactly like FeaturePlot does per split panel
################################################################################
# For each gene, we replicate the logic of:
#    FeaturePlot(..., split.by="Group_TimePoint", min.cutoff="q10", max.cutoff="q90")
# i.e. each (TimePoint,SurvivalGroup) subset gets its own local q10, q90.
#
# We'll do:
#   1) Build a data frame with all cells' raw expression from the 'data' slot
#   2) For each group-timepoint subset, compute q10,q90 & clamp
#   3) Return a vector of final expression values, aligned to the full set of cells

get_split_specific_expr <- function(sobj, gene) {
  # For integrated assay, usually scale.data stores the expression used in FeaturePlot
  if(gene %in% rownames(GetAssayData(sobj, assay="integrated", slot="scale.data"))){
    print(gene)
    expr_raw <- GetAssayData(sobj, assay = "integrated", slot = "scale.data")[gene, ]
  } else {
    expr_raw <- GetAssayData(sobj, assay = "RNA", slot = "data")[gene, ]
  }
  
  
  df_all <- data.frame(
    cell          = names(expr_raw),
    raw_expr      = as.numeric(expr_raw),
    TimePoint     = sobj$TimePoint,
    SurvivalGroup = sobj$SurvivalGroup
  )
  df_all$Group_TimePoint <- paste0(df_all$SurvivalGroup, "_", df_all$TimePoint)
  
  # Subset-specific q10/q90
  df_all$expr_clamped <- df_all$raw_expr
  for (sp in unique(df_all$Group_TimePoint)) {
    idx <- which(df_all$Group_TimePoint == sp)
    if (length(idx) < 1) next
    
    raw_vals <- df_all$raw_expr[idx]
    q10 <- as.numeric(quantile(raw_vals, probs=0.1, na.rm=TRUE))
    q90 <- as.numeric(quantile(raw_vals, probs=0.9, na.rm=TRUE))
    
    print(sp)
    print(q10)
    print(q90)
    
    clamped_vals <- pmin( pmax(raw_vals, q10), q90 )
    df_all$expr_clamped[idx] <- clamped_vals
  }
  
  # Return in same order as the original cell names
  out_vec <- df_all$expr_clamped
  names(out_vec) <- df_all$cell
  return(out_vec)
}


################################################################################
# 4. Plotting 2D Density (3 subplots per gene: Pre, C1, C2)
#    Weighted by the subset-specific clamped expression
################################################################################
plot_2d_density_per_gene <- function(sobj, gene) {
  # color map for groups
  group_colors <- c("Control"="yellow", "ShortTerm"="blue", "LongTerm"="red")
  timepoints   <- c("Pre","C1","C2")
  
  # get UMAP
  umap_mat <- Embeddings(sobj, "umap")[,1:2]
  colnames(umap_mat) <- c("UMAP_1","UMAP_2")
  
  # get expression clamped subset-wise
  expr_clamped <- get_split_specific_expr(sobj, gene)
  
  # build data frame
  plot_df <- data.frame(
    cell          = rownames(umap_mat),
    UMAP_1        = umap_mat[,1],
    UMAP_2        = umap_mat[,2],
    TimePoint     = sobj$TimePoint,
    SurvivalGroup = sobj$SurvivalGroup,
    expression    = expr_clamped
  )
  
  x_rng <- range(plot_df$UMAP_1, na.rm=TRUE)
  y_rng <- range(plot_df$UMAP_2, na.rm=TRUE)
  
  # one subplot per timepoint
  p_list <- list()
  for (tp in timepoints) {
    df_tp <- filter(plot_df, TimePoint == tp)
    plt <- ggplot(df_tp, aes(x=UMAP_1, y=UMAP_2)) +
      geom_point(color="gray70", size=0.4, alpha=0.5) +
      coord_fixed() +
      coord_flip() +
      xlim(x_rng) + ylim(y_rng) +
      theme_void() +
      ggtitle(paste0(gene, " – ", tp))
    
    # add group density contours
    for (grp in c("Control","ShortTerm","LongTerm")) {
      df_grp <- filter(df_tp, SurvivalGroup==grp)
      # if at least a couple cells
      if (nrow(df_grp) >= 2) {
        plt <- plt + geom_density_2d(
          data = df_grp,
          aes(weight=expression),
          color = group_colors[grp],
          size  = 0.7
        )
      }
    }
    p_list[[tp]] <- plt
  }
  
  return(plot_grid(plotlist=p_list, nrow=1))
}

################################################################################
# 5. Wrapper to Plot All Genes
################################################################################
create_2d_density_plots_for_genes <- function(sobj, genes, pdf_file) {
  pdf(pdf_file, width=15, height=5)
  for (g in genes) {
    cat("Plotting:", g, "\n")
    p_gene <- plot_2d_density_per_gene(sobj, g)
    print(p_gene)
  }
  dev.off()
  message("Saved plots to: ", pdf_file)
}

################################################################################
# 6. Try Without and With Downsampling
################################################################################
# (A) No Downsampling
create_2d_density_plots_for_genes(
  sobj     = seurat_tp,
  genes    = immune_checkpoint_genes,
  pdf_file = pdf_no_ds_path
)

# (B) With Downsampling
seurat_tp_ds <- downsample_timepoints(seurat_tp)
create_2d_density_plots_for_genes(
  sobj     = seurat_tp_ds,
  genes    = immune_checkpoint_genes,
  pdf_file = pdf_ds_path
)

message("Done!")



















#########################################################################################
#########################################################################################
# Load necessary libraries
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 1. Define Survivor Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 2. Load & Prepare Survival Data
# -----------------------
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df$patient_id <- as.integer(survival_df$patient_id)

survival_df <- survival_df %>%
  filter(site == "UF") %>%
  mutate(SurvivalGroup = case_when(
    patient_id %in% short_term_survivor_group ~ "ShortTerm",
    patient_id %in% long_term_survivor_group  ~ "LongTerm",
    patient_id %in% control_group             ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(SurvivalGroup))

required_columns <- c("OS.months.", "patient_id")
missing_columns <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: Missing columns in survival_df:", paste(missing_columns, collapse = ", ")))
}
survival_df <- survival_df %>% arrange(OS.months.)

# -----------------------
# 3. Load & Subset Seurat Object
# -----------------------
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# Remove cluster 13 if desired
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

common_patients <- intersect(
  unique(seurat_object_t_cells@meta.data$Patient),
  unique(survival_df$patient_id)
)
if (length(common_patients) == 0) {
  stop("ERROR: No overlapping patients found.")
} else {
  message(paste("Number of overlapping patients:", length(common_patients)))
}

patients_to_keep <- survival_df$patient_id
seurat_subset <- subset(
  seurat_object_t_cells,
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)

message(paste("Number of patients after subsetting:", length(unique(seurat_subset@meta.data$Patient))))
        
# -----------------------
# 4. Add SurvivalGroup Metadata
# -----------------------
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  mutate(CellBarcode = rownames(.))

mapping_df <- survival_df %>%
  select(patient_id, SurvivalGroup) %>%
  rename(Patient = patient_id)

seurat_subset@meta.data <- seurat_subset@meta.data %>%
  left_join(mapping_df, by = "Patient")

rownames(seurat_subset@meta.data) <- seurat_subset@meta.data$CellBarcode
seurat_subset@meta.data$CellBarcode <- NULL

# Factor levels
seurat_subset@meta.data$SurvivalGroup <- factor(
  seurat_subset@meta.data$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------
# 5. Subset to TimePoints of Interest & Downsample
# -----------------------
valid_timepoints <- c("Pre","C1","C2")
seurat_subset <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)

set.seed(123)
final_cells <- c()
for(grp in levels(seurat_subset@meta.data$SurvivalGroup)) {
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  pre_cells <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "Pre"))
  c1_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C1"))
  c2_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C2"))
  
  n_pre <- length(pre_cells)
  n_c1  <- length(c1_cells)
  n_c2  <- length(c2_cells)
  
  n_min <- min(n_pre, n_c1, n_c2)
  if(n_min == 0) {
    warning(paste("Skipping group", grp, " - zero cells in at least one timepoint."))
    next
  }
  
  final_cells <- c(
    final_cells,
    sample(pre_cells, n_min),
    sample(c1_cells, n_min),
    sample(c2_cells, n_min)
  )
}
seurat_timepoint_downsampled <- subset(seurat_subset, cells = final_cells)

cat("Final cell counts after timepoint-downsampling:\n")
print(table(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint
))

# -----------------------
# 6. Create "Group_TimePoint" metadata
# -----------------------
seurat_timepoint_downsampled$Group_TimePoint <- paste(
  seurat_timepoint_downsampled$SurvivalGroup,
  seurat_timepoint_downsampled$TimePoint,
  sep="_"
)
desired_order <- c(
  "Control_Pre","Control_C1","Control_C2",
  "ShortTerm_Pre","ShortTerm_C1","ShortTerm_C2",
  "LongTerm_Pre","LongTerm_C1","LongTerm_C2"
)
seurat_timepoint_downsampled$Group_TimePoint <- factor(
  seurat_timepoint_downsampled$Group_TimePoint,
  levels=desired_order
)

# -----------------------
# 7. Define Genes
# -----------------------
immune_checkpoint_genes <- c(
  "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
  "BTLA","VSIR","CD276","VTCN1","IDO1","IDO2","CD47","CD28",
  "ICOS","TNFRSF4","TNFRSF9","TNFRSF18","NCR3LG1","HHLA2"
)

# -----------------------
# 8. Make FeaturePlot & SAVE the underlying data (manually determine split labels)
# -----------------------
top_labels <- levels(seurat_timepoint_downsampled$Group_TimePoint)
top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=4)
})
top_row <- plot_grid(plotlist=top_label_plots, nrow=1)

all_rows_list  <- list()  # for final PDF
all_feature_df <- list()  # to store data from each splitted panel

for (gene in immune_checkpoint_genes) {
  cat("Default Assay is:", DefaultAssay(seurat_timepoint_downsampled), "\n")
  
  # combine=FALSE => returns a list of ggplots (one per split panel)
  p_list <- FeaturePlot(
    object    = seurat_timepoint_downsampled,
    features  = gene,
    split.by  = "Group_TimePoint",
    cols      = c("lightgrey","red"),
    combine   = FALSE,
    pt.size   = 0.6,
    order     = TRUE,
    min.cutoff= "q10",
    max.cutoff= "q90"
  )
  
  # We'll figure out the actual "Group_TimePoint" subset for each p_list[[i]]
  splitted_data_dfs <- lapply(seq_along(p_list), function(i) {
    df_i <- p_list[[i]]$data  # x=UMAP_1, y=UMAP_2, <gene>, etc.
    df_i$panel_index <- i
    df_i$gene        <- gene
    
    # Here: rownames(df_i) are the cell names of that panel.
    # We can see which "Group_TimePoint" is in seurat_timepoint_downsampled:
    these_cells <- rownames(df_i)
    # Typically, they should all have the same Group_TimePoint if splitted properly:
    gtp <- unique(seurat_timepoint_downsampled$Group_TimePoint[these_cells])
    if (length(gtp) == 1) {
      df_i$split_label <- as.character(gtp)  # e.g. "Control_Pre"
    } else {
      # If for some reason it's multiple or none, store NA
      df_i$split_label <- NA
    }
    df_i
  })
  
  # combine them for this gene
  splitted_data_combined <- dplyr::bind_rows(splitted_data_dfs)
  all_feature_df[[gene]] <- splitted_data_combined
  
  # Cosmetic modifications for the PDF
  p_list <- lapply(p_list, function(x) {
    x + coord_flip() + ggtitle(NULL) + NoLegend() + theme_void()
  })
  gene_row_plots <- plot_grid(plotlist=p_list, nrow=1)
  
  gene_label_plot <- ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=5)
  
  row_for_gene <- plot_grid(
    gene_label_plot,
    gene_row_plots,
    nrow=1, rel_widths=c(0.04,1)
  )
  all_rows_list[[gene]] <- row_for_gene
}

# Combine all rows into one big figure
main_matrix <- plot_grid(plotlist=all_rows_list, ncol=1)
final_plot <- plot_grid(top_row, main_matrix, ncol=1, rel_heights=c(0.05,1))

# Save the PDF
output_pdf_original <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_splitted_featureplots.pdf"
num_genes <- length(immune_checkpoint_genes)
pdf_width  <- 2*10
pdf_height <- 2*num_genes
pdf(output_pdf_original, width=pdf_width, height=pdf_height)
print(final_plot)
dev.off()
cat("Saved FeaturePlots to:", output_pdf_original, "\n")

# Also save the combined FeaturePlot data as RDS
featureplot_data <- dplyr::bind_rows(all_feature_df)
saveRDS(featureplot_data,
        file="/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds"
)
cat("Saved underlying FeaturePlot data to: /project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds\n")
        
        


################################################################################
# 1. Libraries & Input
################################################################################
library(dplyr)
library(ggplot2)
library(cowplot)

df_all <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds")
# This data was generated by the modified Code 1.
# It has columns like x=UMAP_1, y=UMAP_2, the expression columns, plus
# "split_label" = e.g. "Control_Pre"

pdf_output <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_from_FeaturePlotData.pdf"

# ################################################################################
# # A function to create 3 subplots (Pre, C1, C2) from the FeaturePlot data
# ################################################################################
# plot_2d_density <- function(df_all, target_gene) {
#   # 1) Subset the rows for just this gene
#   df_gene <- df_all %>% filter(gene == target_gene)
#   
#   # 2) Figure out which column holds the expression for this gene
#   #    Many genes have "rna_<GENE>", IDO1 might be just "IDO1".
#   if (target_gene == "IDO1") {
#     expr_col <- "IDO1"
#   } else {
#     expr_col <- paste0("rna_", target_gene)
#   }
#   if (! expr_col %in% colnames(df_gene)) {
#     stop("No column named ", expr_col, " for gene=", target_gene)
#   }
#   
#   # Copy that column to "feature" for convenience
#   df_gene$feature <- df_gene[[expr_col]]
#   
#   # Get all cells for gray background (no filtering)
#   df_all_cells <- df_gene
#   
#   # Only keep cells in the top 50% expression for the target gene for density contours
#   expr_threshold <- quantile(df_gene$feature, probs = 0.9, na.rm = TRUE)
#   df_gene <- df_gene %>% filter(feature > expr_threshold)
#   
#   # 3) Parse timepoint and group from "split_label" (e.g., "Control_Pre")
#   df_all_cells$TimePoint <- sub("^.+_", "", df_all_cells$split_label)  # after underscore
#   df_all_cells$Group     <- sub("_.*$", "", df_all_cells$split_label)  # before underscore
#   df_gene$TimePoint      <- sub("^.+_", "", df_gene$split_label)
#   df_gene$Group          <- sub("_.*$", "", df_gene$split_label)
#   
#   # We'll make 3 subplots for Pre, C1, C2
#   timepoints   <- c("Pre", "C1", "C2")
#   group_colors <- c("Control" = "yellow", "ShortTerm" = "blue", "LongTerm" = "red")
#   
#   p_list <- list()
#   for (tp in timepoints) {
#     # Gray background using all cells
#     df_tp_all <- df_all_cells %>% filter(TimePoint == tp)
#     
#     # Density contours using top 50% expression cells
#     df_tp_top <- df_gene %>% filter(TimePoint == tp)
#     
#     # Base plot with gray background
#     p <- ggplot(df_tp_all, aes(x = UMAP_1, y = UMAP_2)) +
#       geom_point(color = "gray70", size = 0.4, alpha = 0.5) +
#       coord_fixed() + 
#       coord_flip() +
#       theme_minimal(base_size = 12) +
#       ggtitle(paste(target_gene, "-", tp))
#     
#     # Overlays: each group gets density, weighted by the expression
#     for (grp in c("Control", "ShortTerm", "LongTerm")) {
#       df_grp <- filter(df_tp_top, Group == grp)
#       if (nrow(df_grp) >= 2) {
#         p <- p + geom_density_2d(
#           data = df_grp,
#           aes(weight = feature),
#           color = group_colors[grp],
#           size  = 0.7
#         )
#       }
#     }
#     p_list[[tp]] <- p
#   }
#   return(plot_grid(plotlist = p_list, nrow = 1))
# }


################################################################################
# A function to create 3 subplots (Pre, C1, C2) from the FeaturePlot data
################################################################################
plot_2d_density <- function(df_all, target_gene) {
  # 1) Subset the rows for just this gene
  df_gene <- df_all %>% filter(gene == target_gene)
  
  # 2) Figure out which column holds the expression for this gene
  if (target_gene == "IDO1") {
    expr_col <- "IDO1"
  } else {
    expr_col <- paste0("rna_", target_gene)
  }
  if (! expr_col %in% colnames(df_gene)) {
    stop("No column named ", expr_col, " for gene=", target_gene)
  }
  
  # Copy that column to "feature" for convenience
  df_gene$feature <- df_gene[[expr_col]]
  
  # Get all cells for gray background (no filtering)
  df_all_cells <- df_gene
  
  # Only keep cells in the top 50% expression for the target gene for density contours
  expr_threshold <- quantile(df_gene$feature, probs = 0.9, na.rm = TRUE)
  df_gene <- df_gene %>% filter(feature > expr_threshold)
  
  # Parse timepoint and group from "split_label" (e.g., "Control_Pre")
  df_all_cells$TimePoint <- sub("^.+_", "", df_all_cells$split_label)  # after underscore
  df_all_cells$Group     <- sub("_.*$", "", df_all_cells$split_label)  # before underscore
  df_gene$TimePoint      <- sub("^.+_", "", df_gene$split_label)
  df_gene$Group          <- sub("_.*$", "", df_gene$split_label)
  
  # We'll make 3 subplots for Pre, C1, C2
  timepoints   <- c("Pre", "C1", "C2")
  group_colors <- c("Control" = "yellow", "ShortTerm" = "blue", "LongTerm" = "red")
  
  p_list <- list()
  for (tp in timepoints) {
    # Gray background using all cells
    df_tp_all <- df_all_cells %>% filter(TimePoint == tp)
    
    # Density contours using top 50% expression cells
    df_tp_top <- df_gene %>% filter(TimePoint == tp)
    
    # Base plot with gray background
    p <- ggplot(df_tp_all, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(color = "gray70", size = 0.4, alpha = 0.5) +
      coord_fixed() + 
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank()) +
      ggtitle(paste(target_gene, "-", tp))
    
    # Overlays: each group gets density, weighted by expression
    for (grp in c("Control", "ShortTerm", "LongTerm")) {
      df_grp <- filter(df_tp_top, Group == grp)
      if (nrow(df_grp) >= 2) {
        # Determine the number of contour lines based on the number of cells
        n_cells <- nrow(df_grp)
        n_bins <- max(3, min(10, floor(n_cells / 20)))  # At least 3 bins, at most 10
        
        p <- p + geom_density_2d(
          data = df_grp,
          aes(weight = feature),
          bins = n_bins,  # Set the number of contour lines
          color = group_colors[grp],
          size  = 0.7
        )
      }
    }
    p_list[[tp]] <- p
  }
  return(plot_grid(plotlist = p_list, nrow = 1))
}

# Loop over your genes and make a page per gene
pdf(pdf_output, width=15, height=5)
for (g in unique(df_all$gene)) {
  cat("Making 2D density for gene:", g, "\n")
  p_out <- plot_2d_density(df_all, target_gene=g)
  print(p_out)
}
dev.off()

cat("Saved 2D Density Plots to:", pdf_output, "\n")









####################################################################################
# graph for individual patients
####################################################################################
################################################################################
# Immune‑checkpoint expression across cohorts & time‑points
#   • Three cohorts: Control, Short‑Term, Long‑Term survivors
#   • Three time‑points: Pre, C1, C2
#   • For every gene, draw a panel of three box‑plots (one per cohort)
#     – each box shows the distribution of patient‑level mean expression
#     – points are patients, connected by lines across time‑points
################################################################################

# -----------------------------------------------------------------------------#
# 1.  Load libraries                                                            #
# -----------------------------------------------------------------------------#
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)     # still used in earlier steps
  library(patchwork)   # for arranging plots; install.packages("patchwork") if needed
  library(tibble) 
})

# -----------------------------------------------------------------------------#
# 2.  Define survivor groups                                                    #
# -----------------------------------------------------------------------------#
control_group            <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------------------------------------------------------------#
# 3.  Load & prepare survival metadata                                         #
# -----------------------------------------------------------------------------#
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"

survival_df <- read.csv(survival_data) |>
  mutate(patient_id = as.integer(patient_id)) |>
  filter(site == "UF") |>
  mutate(
    SurvivalGroup = case_when(
      patient_id %in% short_term_survivor_group ~ "ShortTerm",
      patient_id %in% long_term_survivor_group  ~ "LongTerm",
      patient_id %in% control_group             ~ "Control",
      TRUE                                      ~ NA_character_
    )
  ) |>
  filter(!is.na(SurvivalGroup)) |>
  arrange(OS.months.)

# sanity check
required_columns <- c("OS.months.", "patient_id")
missing_columns  <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0)
  stop("Missing columns in survival_df: ", paste(missing_columns, collapse = ", "))

# -----------------------------------------------------------------------------#
# 4.  Load Seurat object & keep only patients of interest                       #
# -----------------------------------------------------------------------------#
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# Remove cluster 13 if desired
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data))
  stop("'Patient' column not found in Seurat object metadata.")

# keep overlapping patients
patients_to_keep <- intersect(
  unique(seurat_object_t_cells@meta.data$Patient),
  survival_df$patient_id
)
if (length(patients_to_keep) == 0)
  stop("No overlapping patients between Seurat object and survival table.")

seurat_subset <- subset(
  seurat_object_t_cells,
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)

message("Patients after subsetting: ", length(unique(seurat_subset$Patient)))

# -----------------------------------------------------------------------------#
# 5.  Map SurvivalGroup onto Seurat metadata                                    #
# -----------------------------------------------------------------------------#
seurat_subset@meta.data <- seurat_subset@meta.data |>
  tibble::rownames_to_column("CellBarcode") |>
  left_join(
    survival_df |>
      select(patient_id, SurvivalGroup) |>
      rename(Patient = patient_id),
    by = "Patient"
  ) |>
  column_to_rownames("CellBarcode")

seurat_subset$SurvivalGroup <- factor(
  seurat_subset$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------------------------------------------------------------#
# 6.  Keep three key time‑points & down‑sample equal cells per group            #
# -----------------------------------------------------------------------------#
valid_timepoints <- c("Pre", "C1", "C2")
seurat_subset    <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)

set.seed(123)
final_cells <- c()
for (grp in levels(seurat_subset$SurvivalGroup)) {
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  sampled <- lapply(valid_timepoints, function(tp) {
    cells_tp <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == tp))
    cells_tp
  })
  names(sampled) <- valid_timepoints
  
  n_min <- min(lengths(sampled))
  if (n_min == 0) {
    warning("Skipping group ", grp, ": zero cells in at least one time‑point.")
    next
  }
  final_cells <- c(
    final_cells,
    unlist(lapply(sampled, sample, n_min))
  )
}
seurat_tp_ds <- subset(seurat_subset, cells = final_cells)

cat("Cell counts after down‑sampling:\n")
print(table(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint))

# -----------------------------------------------------------------------------#
# 7.  Create combined Group_TimePoint label (used for ordering only)            #
# -----------------------------------------------------------------------------#
seurat_tp_ds$Group_TimePoint <- factor(
  paste(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint, sep = "_"),
  levels = c(
    "Control_Pre",    "Control_C1",    "Control_C2",
    "ShortTerm_Pre",  "ShortTerm_C1",  "ShortTerm_C2",
    "LongTerm_Pre",   "LongTerm_C1",   "LongTerm_C2"
  )
)

# -----------------------------------------------------------------------------#
# 8.  Define genes of interest                                                 #
# -----------------------------------------------------------------------------#
immune_checkpoint_genes <- c(
  "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
  "BTLA","VSIR","CD276","VTCN1","IDO1","IDO2","CD47","CD28",
  "ICOS","TNFRSF4","TNFRSF9","TNFRSF18","NCR3LG1","HHLA2"
)

# -----------------------------------------------------------------------------#
# 9.  Build per‑gene box‑plots                                                 #
# -----------------------------------------------------------------------------#
boxplot_list <- list()

for (gene in immune_checkpoint_genes) {
  
  # pull expression + key metadata
  expr_df <- FetchData(
    seurat_tp_ds,
    vars = c(gene, "Patient", "TimePoint", "SurvivalGroup")
  )
  colnames(expr_df)[1] <- "expr"
  
  # average per patient‑time‑point (mean over that patient’s cells)
  agg_df <- expr_df |>
    group_by(Patient, TimePoint, SurvivalGroup) |>
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
    mutate(
      TimePoint     = factor(TimePoint, levels = valid_timepoints),
      SurvivalGroup = factor(SurvivalGroup,
                             levels = c("Control", "ShortTerm", "LongTerm"))
    )
  
  p <- ggplot(agg_df, aes(x = TimePoint, y = expr)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_point(aes(group = Patient), size = 2,
               position = position_dodge(width = 0.4)) +
    geom_line(aes(group = Patient), alpha = 0.7,
              position = position_dodge(width = 0.4)) +
    facet_wrap(~ SurvivalGroup, nrow = 1) +
    labs(
      title = gene,
      x     = NULL,
      y     = "Mean expression per patient"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      strip.text       = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  boxplot_list[[gene]] <- p
}

# -----------------------------------------------------------------------------#
# 10. Save all plots to a single PDF                                            #
# -----------------------------------------------------------------------------#
output_pdf <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_boxplots.pdf"

pdf(output_pdf, width = 12, height = 2.5 * length(boxplot_list))
wrap_plots(boxplot_list, ncol = 1)
dev.off()

cat("Saved box‑plots to: ", output_pdf, "\n")






##################################################################################
# April 21 - new version
####################################################################################
################################################################################
# Immune‑checkpoint expression across cohorts & time‑points
#   • Three cohorts: Control, Short‑Term, Long‑Term survivors
#   • Three time‑points: Pre, C1, C2
#   • For every gene, draw a panel of three box‑plots (one per cohort)
#     – each box shows the distribution of patient‑level mean expression
#     – points are patients, connected by lines across time‑points
#     – p‑values from paired t‑tests shown for Pre‑C1 and C1‑C2
#   • OPTIONAL: set use_nonzero_only = TRUE to average only over non‑zero cells
################################################################################

# -----------------------------------------------------------------------------#
# 0.  USER OPTION                                                               #
# -----------------------------------------------------------------------------#
use_nonzero_only <- TRUE   # <- set TRUE to ignore zero‑expressing cells

# -----------------------------------------------------------------------------#
# 1.  Load libraries                                                            #
# -----------------------------------------------------------------------------#
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(tibble)
  library(ggpubr)      # for stat_pvalue_manual
})

# -----------------------------------------------------------------------------#
# 2.  Define survivor groups                                                    #
# -----------------------------------------------------------------------------#
control_group            <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------------------------------------------------------------#
# 3.  Load & prepare survival metadata                                         #
# -----------------------------------------------------------------------------#
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"

survival_df <- read.csv(survival_data) |>
  mutate(patient_id = as.integer(patient_id)) |>
  filter(site == "UF") |>
  mutate(
    SurvivalGroup = case_when(
      patient_id %in% short_term_survivor_group ~ "ShortTerm",
      patient_id %in% long_term_survivor_group  ~ "LongTerm",
      patient_id %in% control_group             ~ "Control",
      TRUE                                      ~ NA_character_
    )
  ) |>
  filter(!is.na(SurvivalGroup)) |>
  arrange(OS.months.)

# sanity check
required_columns <- c("OS.months.", "patient_id")
missing_columns  <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0)
  stop("Missing columns in survival_df: ", paste(missing_columns, collapse = ", "))

# -----------------------------------------------------------------------------#
# 0·b  Cluster sets for every T‑cell sub‑population                             #
# -----------------------------------------------------------------------------#
celltype_cluster_map <- list(
  "all"                       = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
  "Activated_CD4"             = c(0),
  "Effector_CD8"              = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Exhausted_T"               = c(3),
  "Gamma_Delta_T"             = c(4),
  "Active_CD4"                = c(5),
  "Naive_CD4"                 = c(6, 9, 18),
  "Memory_CD4"                = c(7),
  "Stem_Like_CD8"             = c(8),
  "Effector_Memory_CD8"       = c(10),
  "Central_Memory_CD8"        = c(12),
  "Proliferating_Effector"    = c(14, 16, 17),
  "All_CD8"                   = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)


# -----------------------------------------------------------------------------#
# 4.  Load Seurat object & keep only patients of interest                       #
# -----------------------------------------------------------------------------#
seurat_t_cells_all <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# Remove cluster 13 if desired
seurat_t_cells_all <- subset(seurat_t_cells_all, subset = seurat_clusters != 13)

if (!"Patient" %in% colnames(seurat_t_cells_all@meta.data))
  stop("'Patient' column not found in Seurat object metadata.")

# keep overlapping patients
patients_to_keep <- intersect(
  unique(seurat_t_cells_all@meta.data$Patient),
  survival_df$patient_id
)

for (celltype in names(celltype_cluster_map)) {
  
  cluster_ids <- celltype_cluster_map[[celltype]]
  
  # --- 4·b·1  subset to the desired clusters AND patients --------------------#
  seurat_subset <- subset(
    seurat_t_cells_all,
    subset = (seurat_clusters %in% cluster_ids) &
      (Patient %in% patients_to_keep)
  )
  
  message(">>>  Working on: ", celltype,
          "   |  cells = ", ncol(seurat_subset))
  
  # -----------------------------------------------------------------------------#
  # 5.  Map SurvivalGroup onto Seurat metadata                                    #
  # -----------------------------------------------------------------------------#
  seurat_subset@meta.data <- seurat_subset@meta.data |>
    tibble::rownames_to_column("CellBarcode") |>
    left_join(
      survival_df |>
        select(patient_id, SurvivalGroup) |>
        rename(Patient = patient_id),
      by = "Patient"
    ) |>
    column_to_rownames("CellBarcode")
  
  seurat_subset$SurvivalGroup <- factor(
    seurat_subset$SurvivalGroup,
    levels = c("Control", "ShortTerm", "LongTerm")
  )
  
  # -----------------------------------------------------------------------------#
  # 6.  Keep three key time‑points & down‑sample equal cells per group            #
  # -----------------------------------------------------------------------------#
  valid_timepoints <- c("Pre", "C1", "C2")
  seurat_subset    <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)
  
  set.seed(123)
  final_cells <- c()
  for (grp in levels(seurat_subset$SurvivalGroup)) {
    grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
    
    sampled <- lapply(valid_timepoints, function(tp) {
      cells_tp <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == tp))
      cells_tp
    })
    names(sampled) <- valid_timepoints
    
    n_min <- min(lengths(sampled))
    if (n_min == 0) {
      warning("Skipping group ", grp, ": zero cells in at least one time‑point.")
      next
    }
    final_cells <- c(
      final_cells,
      unlist(lapply(sampled, sample, n_min))
    )
  }
  seurat_tp_ds <- subset(seurat_subset, cells = final_cells)
  
  cat("Cell counts after down‑sampling:\n")
  print(table(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint))
  
  # -----------------------------------------------------------------------------#
  # 7.  Create combined Group_TimePoint label (used for ordering only)            #
  # -----------------------------------------------------------------------------#
  seurat_tp_ds$Group_TimePoint <- factor(
    paste(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint, sep = "_"),
    levels = c(
      "Control_Pre",    "Control_C1",    "Control_C2",
      "ShortTerm_Pre",  "ShortTerm_C1",  "ShortTerm_C2",
      "LongTerm_Pre",   "LongTerm_C1",   "LongTerm_C2"
    )
  )
  
  # -----------------------------------------------------------------------------#
  # 8.  Define genes of interest                                                 #
  # -----------------------------------------------------------------------------#
  immune_checkpoint_genes <- c(
    "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
    "BTLA","VSIR","CD276","VTCN1","IDO1","IDO2","CD47","CD28",
    "ICOS","TNFRSF4","TNFRSF9","TNFRSF18","NCR3LG1","HHLA2"
  )
  
  # helper to compute mean with / without zeros
  mean_expr <- function(x) {
    if (use_nonzero_only) {
      nz <- x[x > 0]
      if (length(nz) == 0) return(0)
      return(mean(nz))
    } else {
      return(mean(x))
    }
  }
  
  # -----------------------------------------------------------------------------#
  # 9.  Build per‑gene box‑plots (fixed facet order + uncluttered p‑values)       #
  # -----------------------------------------------------------------------------#
  
  boxplot_list <- list()
  
  for (gene in immune_checkpoint_genes) {
    
    # ── 9·1  Patient‑level means ────────────────────────────────────────────────
    expr_df <- FetchData(
      seurat_tp_ds,
      vars = c(gene, "Patient", "TimePoint", "SurvivalGroup")
    )
    colnames(expr_df)[1] <- "expr"
    
    agg_df <- expr_df |>
      group_by(Patient, TimePoint, SurvivalGroup) |>
      summarise(expr = mean_expr(expr), .groups = "drop") |>
      mutate(
        TimePoint     = factor(TimePoint, levels = valid_timepoints),
        SurvivalGroup = factor(SurvivalGroup,
                               levels = c("Control", "ShortTerm", "LongTerm")) # <- explicit order
      )
    
    # ── 9·2  Paired t‑tests for Pre‑C1 and C1‑C2 in each cohort ────────────────
    stat_tests <- data.frame()
    global_max <- max(agg_df$expr, na.rm = TRUE)            # for nicer spacing
    
    for (grp in levels(agg_df$SurvivalGroup)) {
      
      wide_df <- agg_df |>
        filter(SurvivalGroup == grp) |>
        select(Patient, TimePoint, expr) |>
        pivot_wider(names_from = TimePoint, values_from = expr)
      
      # paired Pre vs C1
      p1 <- { if (all(c("Pre","C1") %in% names(wide_df))) {
        ok <- complete.cases(wide_df$Pre, wide_df$C1)
        if (sum(ok) >= 2) t.test(wide_df$Pre[ok], wide_df$C1[ok],
                                 paired = TRUE)$p.value else NA_real_
      } else NA_real_ }
      
      # paired C1 vs C2
      p2 <- { if (all(c("C1","C2") %in% names(wide_df))) {
        ok <- complete.cases(wide_df$C1, wide_df$C2)
        if (sum(ok) >= 2) t.test(wide_df$C1[ok], wide_df$C2[ok],
                                 paired = TRUE)$p.value else NA_real_
      } else NA_real_ }
      
      # y positions sufficiently above the highest point in *this* cohort,
      # plus extra head‑room so the strip label never overlaps.
      y_cohort_max <- max(agg_df$expr[agg_df$SurvivalGroup == grp], na.rm = TRUE)
      y_start      <- max(y_cohort_max, global_max * 0.05)  # avoid tiny offsets
      stat_tests   <- rbind(
        stat_tests,
        data.frame(SurvivalGroup = grp,
                   group1 = "Pre", group2 = "C1",
                   y.position = y_start + 0.10 * global_max,
                   p.value = p1),
        data.frame(SurvivalGroup = grp,
                   group1 = "C1", group2 = "C2",
                   y.position = y_start + 0.20 * global_max,
                   p.value = p2)
      )
    }
    
    stat_tests$label <- dplyr::case_when(
      is.na(stat_tests$p.value)          ~ "NA",
      stat_tests$p.value < 0.001         ~ "p < 0.001",
      TRUE                               ~ paste0("p = ", signif(stat_tests$p.value, 2))
    )
    
    # ── 9·3  Build the plot ─────────────────────────────────────────────────────
    p <- ggplot(agg_df, aes(x = TimePoint, y = expr)) +
      geom_boxplot(outlier.shape = NA, fill = "grey90") +
      geom_point(aes(group = Patient), size = 2,
                 position = position_dodge(width = 0.4)) +
      geom_line(aes(group = Patient), alpha = 0.7,
                position = position_dodge(width = 0.4)) +
      facet_wrap(
        ~ factor(SurvivalGroup,
                 levels = c("Control", "ShortTerm", "LongTerm")),
        nrow = 1, drop = FALSE
      ) +
      ggpubr::stat_pvalue_manual(
        stat_tests,
        label = "label",
        xmin  = "group1",
        xmax  = "group2",
        bracket.size = 0.4,
        tip.length   = 0,
        hide.ns      = FALSE
      ) +
      # extra head‑room to ensure brackets/labels never collide with the strip
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
      labs(
        title = gene,
        x     = NULL,
        y     = paste0(
          "Mean expression per patient",
          ifelse(use_nonzero_only, " (non‑zero cells only)", "")
        )
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        # allow annotation to extend into the extra top margin
        plot.margin        = margin(t = 10, r = 5, b = 5, l = 5),
        panel.spacing.x    = unit(1.1, "lines")    # little breathing‑room between facets
      ) +
      coord_cartesian(clip = "off")   # don’t clip the p‑value brackets
    
    boxplot_list[[gene]] <- p
  }
  
  # --- 10·b  Save one PDF per cell‑type ------------------------------------- #
  suffix     <- ifelse(use_nonzero_only, "_nonzero", "_allcells")
  output_pdf <- file.path(
    "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_boxplots",
    paste0(celltype, suffix, ".pdf")
  )
  
  ## ① ensure directory exists
  out_dir <- dirname(output_pdf)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ## ② skip if nothing to plot
  if (length(boxplot_list) == 0) {
    warning("No plots for ", celltype, "; skipping.")
    next
  }
  
  ## ③ open device, PRINT the plot, close
  pdf(output_pdf, width = 12, height = 2.8 * length(boxplot_list))
  print(wrap_plots(boxplot_list, ncol = 1))   # <- must be printed
  dev.off()
  
  cat("Saved box‑plots for **", celltype, "** to: ", output_pdf, "\n\n")
  
  
} # end for‑celltype loop





























# CD8 T Cell clusters
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters %in% c(1, 2, 3, 8, 10, 12, 14, 16, 17))

if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data))
  stop("'Patient' column not found in Seurat object metadata.")

# keep overlapping patients
patients_to_keep <- intersect(
  unique(seurat_object_t_cells@meta.data$Patient),
  survival_df$patient_id
)
if (length(patients_to_keep) == 0)
  stop("No overlapping patients between Seurat object and survival table.")

seurat_subset <- subset(
  seurat_object_t_cells,
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)

message("Patients after subsetting: ", length(unique(seurat_subset$Patient)))

# -----------------------------------------------------------------------------#
# 5.  Map SurvivalGroup onto Seurat metadata                                    #
# -----------------------------------------------------------------------------#
seurat_subset@meta.data <- seurat_subset@meta.data |>
  tibble::rownames_to_column("CellBarcode") |>
  left_join(
    survival_df |>
      select(patient_id, SurvivalGroup) |>
      rename(Patient = patient_id),
    by = "Patient"
  ) |>
  column_to_rownames("CellBarcode")

seurat_subset$SurvivalGroup <- factor(
  seurat_subset$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------------------------------------------------------------#
# 6.  Keep three key time‑points & down‑sample equal cells per group            #
# -----------------------------------------------------------------------------#
valid_timepoints <- c("Pre", "C1", "C2")
seurat_subset    <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)

set.seed(123)
final_cells <- c()
for (grp in levels(seurat_subset$SurvivalGroup)) {
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  sampled <- lapply(valid_timepoints, function(tp) {
    cells_tp <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == tp))
    cells_tp
  })
  names(sampled) <- valid_timepoints
  
  n_min <- min(lengths(sampled))
  if (n_min == 0) {
    warning("Skipping group ", grp, ": zero cells in at least one time‑point.")
    next
  }
  final_cells <- c(
    final_cells,
    unlist(lapply(sampled, sample, n_min))
  )
}
seurat_tp_ds <- subset(seurat_subset, cells = final_cells)

cat("Cell counts after down‑sampling:\n")
print(table(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint))

# -----------------------------------------------------------------------------#
# 7.  Create combined Group_TimePoint label (used for ordering only)            #
# -----------------------------------------------------------------------------#
seurat_tp_ds$Group_TimePoint <- factor(
  paste(seurat_tp_ds$SurvivalGroup, seurat_tp_ds$TimePoint, sep = "_"),
  levels = c(
    "Control_Pre",    "Control_C1",    "Control_C2",
    "ShortTerm_Pre",  "ShortTerm_C1",  "ShortTerm_C2",
    "LongTerm_Pre",   "LongTerm_C1",   "LongTerm_C2"
  )
)

# -----------------------------------------------------------------------------#
# 8.  Define genes of interest                                                 #
# -----------------------------------------------------------------------------#
immune_checkpoint_genes <- c(
  "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
  "BTLA","VSIR","CD276","VTCN1","IDO1","IDO2","CD47","CD28",
  "ICOS","TNFRSF4","TNFRSF9","TNFRSF18","NCR3LG1","HHLA2"
)

# helper to compute mean with / without zeros
mean_expr <- function(x) {
  if (use_nonzero_only) {
    nz <- x[x > 0]
    if (length(nz) == 0) return(0)
    return(mean(nz))
  } else {
    return(mean(x))
  }
}

# -----------------------------------------------------------------------------#
# 9.  Build per‑gene box‑plots (fixed facet order + uncluttered p‑values)       #
# -----------------------------------------------------------------------------#

boxplot_list <- list()

for (gene in immune_checkpoint_genes) {
  
  # ── 9·1  Patient‑level means ────────────────────────────────────────────────
  expr_df <- FetchData(
    seurat_tp_ds,
    vars = c(gene, "Patient", "TimePoint", "SurvivalGroup")
  )
  colnames(expr_df)[1] <- "expr"
  
  agg_df <- expr_df |>
    group_by(Patient, TimePoint, SurvivalGroup) |>
    summarise(expr = mean_expr(expr), .groups = "drop") |>
    mutate(
      TimePoint     = factor(TimePoint, levels = valid_timepoints),
      SurvivalGroup = factor(SurvivalGroup,
                             levels = c("Control", "ShortTerm", "LongTerm")) # <- explicit order
    )
  
  # ── 9·2  Paired t‑tests for Pre‑C1 and C1‑C2 in each cohort ────────────────
  stat_tests <- data.frame()
  global_max <- max(agg_df$expr, na.rm = TRUE)            # for nicer spacing
  
  for (grp in levels(agg_df$SurvivalGroup)) {
    
    wide_df <- agg_df |>
      filter(SurvivalGroup == grp) |>
      select(Patient, TimePoint, expr) |>
      pivot_wider(names_from = TimePoint, values_from = expr)
    
    # paired Pre vs C1
    p1 <- { if (all(c("Pre","C1") %in% names(wide_df))) {
      ok <- complete.cases(wide_df$Pre, wide_df$C1)
      if (sum(ok) >= 2) t.test(wide_df$Pre[ok], wide_df$C1[ok],
                               paired = TRUE)$p.value else NA_real_
    } else NA_real_ }
    
    # paired C1 vs C2
    p2 <- { if (all(c("C1","C2") %in% names(wide_df))) {
      ok <- complete.cases(wide_df$C1, wide_df$C2)
      if (sum(ok) >= 2) t.test(wide_df$C1[ok], wide_df$C2[ok],
                               paired = TRUE)$p.value else NA_real_
    } else NA_real_ }
    
    # y positions sufficiently above the highest point in *this* cohort,
    # plus extra head‑room so the strip label never overlaps.
    y_cohort_max <- max(agg_df$expr[agg_df$SurvivalGroup == grp], na.rm = TRUE)
    y_start      <- max(y_cohort_max, global_max * 0.05)  # avoid tiny offsets
    stat_tests   <- rbind(
      stat_tests,
      data.frame(SurvivalGroup = grp,
                 group1 = "Pre", group2 = "C1",
                 y.position = y_start + 0.10 * global_max,
                 p.value = p1),
      data.frame(SurvivalGroup = grp,
                 group1 = "C1", group2 = "C2",
                 y.position = y_start + 0.20 * global_max,
                 p.value = p2)
    )
  }
  
  stat_tests$label <- dplyr::case_when(
    is.na(stat_tests$p.value)          ~ "NA",
    stat_tests$p.value < 0.001         ~ "p < 0.001",
    TRUE                               ~ paste0("p = ", signif(stat_tests$p.value, 2))
  )
  
  # ── 9·3  Build the plot ─────────────────────────────────────────────────────
  p <- ggplot(agg_df, aes(x = TimePoint, y = expr)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_point(aes(group = Patient), size = 2,
               position = position_dodge(width = 0.4)) +
    geom_line(aes(group = Patient), alpha = 0.7,
              position = position_dodge(width = 0.4)) +
    facet_wrap(
      ~ factor(SurvivalGroup,
               levels = c("Control", "ShortTerm", "LongTerm")),
      nrow = 1, drop = FALSE
    ) +
    ggpubr::stat_pvalue_manual(
      stat_tests,
      label = "label",
      xmin  = "group1",
      xmax  = "group2",
      bracket.size = 0.4,
      tip.length   = 0,
      hide.ns      = FALSE
    ) +
    # extra head‑room to ensure brackets/labels never collide with the strip
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    labs(
      title = gene,
      x     = NULL,
      y     = paste0(
        "Mean expression per patient",
        ifelse(use_nonzero_only, " (non‑zero cells only)", "")
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      strip.text         = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      # allow annotation to extend into the extra top margin
      plot.margin        = margin(t = 10, r = 5, b = 5, l = 5),
      panel.spacing.x    = unit(1.1, "lines")    # little breathing‑room between facets
    ) +
    coord_cartesian(clip = "off")   # don’t clip the p‑value brackets
  
  boxplot_list[[gene]] <- p
}

# -----------------------------------------------------------------------------#
# 10. Save all plots to a single PDF                                            #
# -----------------------------------------------------------------------------#
# suffix <- ifelse(use_nonzero_only, "_nonzero", "_allcells")
suffix <- ifelse(use_nonzero_only, "_nonzero_CD8", "_allcells_CD8")
output_pdf <- paste0(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_boxplots",
  suffix, ".pdf"
)

pdf(output_pdf, width = 12, height = 2.8 * length(boxplot_list))
wrap_plots(boxplot_list, ncol = 1)
dev.off()

cat("Saved box‑plots to: ", output_pdf, "\n")

