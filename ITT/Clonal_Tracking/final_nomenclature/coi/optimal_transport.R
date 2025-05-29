# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Define the mapping with unique subpopulation names
# mapping <- list(
#   "Activated_CD4" = c(0),
#   "Effector_CD8" = c(1),
#   "Effector_Memory_CD8" = c(2),
#   "Exhausted_T" = c(3),
#   "Gamma_Delta_T" = c(4),
#   "Active_CD4" = c(5),
#   "Naive_CD4" = c(6, 9, 18),
#   "Memory_CD4" = c(7),
#   "Memory_CD8" = c(8),
#   "Anergic_CD8" = c(10),
#   "Naive_CD8" = c(12),
#   # "Hyperactivated_CD8" = c(13),
#   "Proliferating_Effector" = c(14, 16, 17),
#   "CD8" = c(1, 2, 13, 14, 16, 17)
# )

mapping <- list(
  "Effector_CD8" = c(1),
  "Memory_Precursor_Effector_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "Proliferating_Effector" = c(14, 16, 17),
  "All_CD8" = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)

# Ordered timepoints
all_timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")
all_patients <- c(1, 4, 8, 9, 7, 10, 12, 14, 18, 2, 3, 5, 13, 19, 20, 21)

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF" & Patient %in% all_patients)
# remove the patient specific seurat cluster
full_seurat_obj <- subset(full_seurat_obj, subset = seurat_clusters != 13)
# keep only clusters of interest
full_seurat_obj <- subset(full_seurat_obj, subset = seurat_clusters %in% c(1, 2, 3, 8, 10, 12, 14, 16, 17))

# Define the cohorts
control_group <- c(1, 4, 8, 9)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Add SurvivorGroup metadata to the full_seurat_obj
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term",
         ifelse(full_seurat_obj$Patient %in% control_group, "control", NA))
)

# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/final_nomenclature/coi"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"

# # Function to extract data for optimal transport visualization with clone identification
# extractDataForOptimalTransport <- function(subpop_name, clusters, cohort_name, patient_ids, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder) {
#   
#   # Subset the Seurat object for the cohort
#   cohort_seurat_obj <- subset(full_seurat_obj, Patient %in% patient_ids)
#   
#   # Subset for the subpopulation
#   subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
#   
#   # Check if there are any cells in the subpopulation
#   if (ncol(subpop_seurat_obj) == 0) {
#     message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Read in clone_cell_barcode_file
#   clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
#   clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
#   
#   # Ensure 'barcode' and 'cdr3' columns exist
#   if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
#     stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
#   }
#   
#   # Read in sample_clone_file
#   sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
#   colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))
#   
#   # Ensure sample IDs are in columns
#   if (ncol(sample_clone_df) == 0) {
#     stop("sample_clone_file must contain sample IDs as columns.")
#   }
#   
#   # Get UMAP coordinates
#   umap_coords <- Embeddings(cohort_seurat_obj, "umap")
#   
#   # Prepare data frame
#   umap_df <- data.frame(
#     cell_id = rownames(cohort_seurat_obj@meta.data),
#     UMAP_1 = umap_coords[,1],
#     UMAP_2 = umap_coords[,2],
#     TimePoint = cohort_seurat_obj@meta.data$TimePoint,
#     SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup,
#     seurat_clusters = cohort_seurat_obj@meta.data$seurat_clusters,
#     origin = cohort_seurat_obj@meta.data$origin
#   )
#   
#   # Identify cell barcodes in the clusters (subpopulation) across all timepoints
#   cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
#   if (length(cells_in_subpop) == 0) {
#     # No cells in subpopulation
#     message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Get clones corresponding to these cell barcodes from clone_cell_barcode_df
#   clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
#   clones_from_cells <- unique(clones_from_cells)
#   clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
#   
#   if (length(clones_from_cells) == 0) {
#     # No clones found for the subpopulation
#     message("No clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Identify clones with proportion > 0 in samples corresponding to the subpopulation
#   sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
#   
#   # Get clones with proportion > 0 in these samples
#   clones_from_samples <- c()
#   for (sample_id in sample_ids) {
#     if (sample_id %in% colnames(sample_clone_df)) {
#       clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
#       clones_from_samples <- c(clones_from_samples, clones_in_sample)
#     }
#   }
#   clones_from_samples <- unique(clones_from_samples)
#   
#   # Take intersection of the two
#   intersected_clones <- intersect(clones_from_cells, clones_from_samples)
#   
#   if (length(intersected_clones) == 0) {
#     # No intersected clones
#     message("No intersected clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
#   all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
#   # Remove NAs and duplicates
#   all_cells_for_clones <- unique(all_cells_for_clones)
#   all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
#   
#   # These are the cells to highlight
#   cells_to_plot <- all_cells_for_clones
#   
#   cohort_timepoints <- intersect(all_timepoints, unique(subpop_seurat_obj@meta.data$TimePoint))
#   
#   # Loop over timepoints except the last one
#   for (i in 1:(length(cohort_timepoints)-1)) {
#     source_tp <- cohort_timepoints[i]
#     target_tp <- cohort_timepoints[i+1]
#     
#     # Get source cells (highlighted cells at source timepoint)
#     source_cells <- umap_df %>%
#       filter(TimePoint == source_tp & cell_id %in% cells_to_plot)
#     
#     # Get target cells (highlighted cells at target timepoint)
#     target_cells <- umap_df %>%
#       filter(TimePoint == target_tp & cell_id %in% cells_to_plot)
#     
#     # Check if we have both source and target cells
#     if (nrow(source_cells) == 0 | nrow(target_cells) == 0) {
#       message("No cells found for source timepoint: ", source_tp, " or target timepoint: ", target_tp, " for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
#       next
#     }
#     
#     # Get background gray cells at source timepoint
#     gray_cells <- umap_df %>%
#       filter(TimePoint == source_tp & !(cell_id %in% cells_to_plot))
#     
#     # Save data to CSV files
#     output_subfolder <- file.path(output_folder, paste0("Timepoint_", source_tp))
#     if (!dir.exists(output_subfolder)) {
#       dir.create(output_subfolder, recursive = TRUE)
#     }
#     
#     write.csv(source_cells, file = file.path(output_subfolder, "source_cells.csv"), row.names = FALSE)
#     write.csv(target_cells, file = file.path(output_subfolder, "target_cells.csv"), row.names = FALSE)
#     write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells.csv"), row.names = FALSE)
#   }
#   
#   # Handle last timepoint (only source cells and gray cells)
#   last_tp <- cohort_timepoints[length(cohort_timepoints)]
#   source_cells <- umap_df %>%
#     filter(TimePoint == last_tp & cell_id %in% cells_to_plot)
#   
#   if (nrow(source_cells) == 0) {
#     message("No cells found for timepoint: ", last_tp, " for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
#     return()
#   }
#   
#   gray_cells <- umap_df %>%
#     filter(TimePoint == last_tp & !(cell_id %in% cells_to_plot))
#   
#   output_subfolder <- file.path(output_folder, paste0("Timepoint_", last_tp))
#   if (!dir.exists(output_subfolder)) {
#     dir.create(output_subfolder, recursive = TRUE)
#   }
#   
#   write.csv(source_cells, file = file.path(output_subfolder, "source_cells.csv"), row.names = FALSE)
#   write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells.csv"), row.names = FALSE)
# }
# 
# # Loop over subpopulations and cohorts
# cohorts <- list(
#   "control" = control_group,
#   "short_term" = short_term_survivor_group,
#   "long_term" = long_term_survivor_group
# )
# 
# for (subpop_name in names(mapping)) {
#   clusters <- mapping[[subpop_name]]
#   
#   for (cohort_name in names(cohorts)) {
#     patient_ids <- cohorts[[cohort_name]]
#     
#     output_folder <- file.path(base_output_dir, subpop_name, cohort_name)
#     if (!dir.exists(output_folder)) {
#       dir.create(output_folder, recursive = TRUE)
#     }
#     
#     # Call the function
#     extractDataForOptimalTransport(subpop_name, clusters, cohort_name, patient_ids, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
#   }
# }






# Modified Function to extract data for optimal transport visualization with clone identification,
# retaining both UMAP and PCA embeddings (all PCs included)
extractDataForOptimalTransport <- function(subpop_name, clusters, cohort_name, patient_ids, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder) {
  
  # Subset the Seurat object for the cohort
  cohort_seurat_obj <- subset(full_seurat_obj, Patient %in% patient_ids)
  
  # Subset for the subpopulation
  subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in the subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Read in clone_cell_barcode_file
  clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
  clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
  
  # Ensure 'barcode' and 'cdr3' columns exist
  if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
    stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
  }
  
  # Read in sample_clone_file
  sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
  colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))
  
  # Ensure sample IDs are in columns
  if (ncol(sample_clone_df) == 0) {
    stop("sample_clone_file must contain sample IDs as columns.")
  }
  
  # Get UMAP coordinates
  umap_coords <- Embeddings(cohort_seurat_obj, "umap")
  
  # Get PCA coordinates (all PCs)
  pca_coords <- Embeddings(cohort_seurat_obj, "pca")
  
  # Assign column names to PCA coordinates if needed
  pca_col_names <- paste0("PC_", seq_len(ncol(pca_coords)))
  colnames(pca_coords) <- pca_col_names
  
  # Prepare data frame with both UMAP and PCA
  combined_df <- data.frame(
    cell_id = rownames(cohort_seurat_obj@meta.data),
    TimePoint = cohort_seurat_obj@meta.data$TimePoint,
    SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup,
    seurat_clusters = cohort_seurat_obj@meta.data$seurat_clusters,
    origin = cohort_seurat_obj@meta.data$origin,
    Patient = cohort_seurat_obj@meta.data$Patient,
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    pca_coords
  )
  
  # Identify cell barcodes in the clusters (subpopulation)
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    # No cells in subpopulation
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Get clones corresponding to these cell barcodes
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Identify clones with proportion > 0 in samples corresponding to the subpopulation
  sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
  
  # Get clones with proportion > 0 in these samples
  clones_from_samples <- c()
  for (sample_id in sample_ids) {
    if (sample_id %in% colnames(sample_clone_df)) {
      clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
      clones_from_samples <- c(clones_from_samples, clones_in_sample)
    }
  }
  clones_from_samples <- unique(clones_from_samples)
  
  # Take intersection of the two
  intersected_clones <- intersect(clones_from_cells, clones_from_samples)
  
  if (length(intersected_clones) == 0) {
    # No intersected clones
    message("No intersected clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Identify all cell barcodes corresponding to the intersected set of clones
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells to highlight
  cells_to_plot <- all_cells_for_clones
  
  # Assuming 'all_timepoints' is defined globally or earlier in your script
  cohort_timepoints <- intersect(all_timepoints, unique(subpop_seurat_obj@meta.data$TimePoint))
  
  # Loop over timepoints except the last one
  for (i in 1:(length(cohort_timepoints)-1)) {
    source_tp <- cohort_timepoints[i]
    target_tp <- cohort_timepoints[i+1]
    
    # Get source cells (highlighted cells at source timepoint)
    source_cells <- combined_df %>%
      dplyr::filter(TimePoint == source_tp & cell_id %in% cells_to_plot)
    
    # Get target cells (highlighted cells at target timepoint)
    target_cells <- combined_df %>%
      dplyr::filter(TimePoint == target_tp & cell_id %in% cells_to_plot)
    
    # Check if we have both source and target cells
    if (nrow(source_cells) == 0 | nrow(target_cells) == 0) {
      message("No cells found for source timepoint: ", source_tp, " or target timepoint: ", target_tp, " for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
      next
    }
    
    # Get background gray cells at source timepoint
    gray_cells <- combined_df %>%
      dplyr::filter(TimePoint == source_tp & !(cell_id %in% cells_to_plot))
    
    # Save data to CSV files
    output_subfolder <- file.path(output_folder, paste0("Timepoint_", source_tp))
    if (!dir.exists(output_subfolder)) {
      dir.create(output_subfolder, recursive = TRUE)
    }
    
    write.csv(source_cells, file = file.path(output_subfolder, "source_cells_new.csv"), row.names = FALSE)
    write.csv(target_cells, file = file.path(output_subfolder, "target_cells_new.csv"), row.names = FALSE)
    write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells_new.csv"), row.names = FALSE)
  }
  
  # Handle last timepoint (only source cells and gray cells)
  last_tp <- cohort_timepoints[length(cohort_timepoints)]
  source_cells <- combined_df %>%
    dplyr::filter(TimePoint == last_tp & cell_id %in% cells_to_plot)
  
  if (nrow(source_cells) == 0) {
    message("No cells found for timepoint: ", last_tp, " for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  gray_cells <- combined_df %>%
    dplyr::filter(TimePoint == last_tp & !(cell_id %in% cells_to_plot))
  
  output_subfolder <- file.path(output_folder, paste0("Timepoint_", last_tp))
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder, recursive = TRUE)
  }
  
  write.csv(source_cells, file = file.path(output_subfolder, "source_cells_new.csv"), row.names = FALSE)
  write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells_new.csv"), row.names = FALSE)
}


# Loop over subpopulations and cohorts remains unchanged
cohorts <- list(
  "control" = control_group,
  "short_term" = short_term_survivor_group,
  "long_term" = long_term_survivor_group
)

for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  for (cohort_name in names(cohorts)) {
    patient_ids <- cohorts[[cohort_name]]
    
    output_folder <- file.path(base_output_dir, subpop_name, cohort_name)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Call the modified function
    extractDataForOptimalTransport(subpop_name, clusters, cohort_name, patient_ids, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }
}
