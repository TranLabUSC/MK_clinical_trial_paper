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


