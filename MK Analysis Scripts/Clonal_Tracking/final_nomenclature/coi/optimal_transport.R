################################################################################
# Optimal Transport Data Preparation Script
#
# PURPOSE:
# This script prepares source and target cell data for optimal transport analysis
# of T cell clones across timepoints. It extracts UMAP and PCA coordinates,
# identifies subpopulation-enriched TCR clones, and saves structured CSV files
# that are consumed by the Python notebook (t_cell_optimal_transport.ipynb) for
# optimal transport visualization and analysis.
#
# BIOLOGICAL RATIONALE:
# Optimal transport theory predicts cell fate transitions by finding the most
# efficient transformation from one cell distribution (source timepoint) to
# another (target timepoint). By tracking TCR clones (cells with identical
# CDR3β sequences), we can:
# - Visualize spatial migration patterns of clonally-related T cells
# - Identify transcriptional state transitions during treatment
# - Compare clone dynamics between patient survival cohorts
# - Understand how antigen-specific T cell responses evolve over time
#
# WORKFLOW:
# 1. Load T cell Seurat object and filter for UF site patients
# 2. Define CD8+ T cell subpopulations and patient survival cohorts
# 3. For each subpopulation × cohort combination:
#    a. Identify subpopulation-enriched TCR clones
#    b. Extract all cells bearing these clones (across all timepoints)
#    c. For consecutive timepoint pairs (Pre→C1, C1→C2, etc.):
#       - Save source cells (timepoint A) with UMAP/PCA coordinates
#       - Save target cells (timepoint B) with UMAP/PCA coordinates
#       - Save background cells (all other cells at timepoint A)
# 4. Output organized CSV files for Python optimal transport notebook
#
# OUTPUT STRUCTURE:
# base_output_dir/
# ├── Effector_CD8/
# │   ├── control/
# │   │   ├── Timepoint_Pre/
# │   │   │   ├── source_cells_new.csv (highlighted clones at Pre)
# │   │   │   ├── target_cells_new.csv (highlighted clones at C1)
# │   │   │   └── gray_cells_new.csv (background cells at Pre)
# │   │   ├── Timepoint_C1/
# │   │   │   ├── source_cells_new.csv
# │   │   │   ├── target_cells_new.csv
# │   │   │   └── gray_cells_new.csv
# │   │   └── ... (subsequent timepoints)
# │   ├── short_term/
# │   └── long_term/
# └── ... (other subpopulations)
#
# DOWNSTREAM USAGE:
# The CSV files generated here are loaded by t_cell_optimal_transport.ipynb,
# which performs optimal transport calculations, generates arrow plots showing
# cell fate transitions, and creates density shift visualizations.
################################################################################

# Load required libraries
library(Seurat)      # Single-cell analysis framework
library(dplyr)       # Data manipulation
library(ggplot2)     # Plotting (though not used in this data prep script)

# Define CD8+ T cell subpopulation cluster mappings
# Maps cell type annotations to Seurat cluster IDs
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
  
  # STEP 3: Load TCR barcode-to-clone mapping
  # Maps each cell barcode to its TCR CDR3β amino acid sequence (clone identifier)
  clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
  
  # Filter for TRB chain only (β chain, more diverse than α chain)
  # TCR clonality is primarily defined by CDR3β sequence
  clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
  
  # Validate required columns exist in TCR data
  if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
    stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
  }
  
  # STEP 4: Load clone proportion table
  # Contains quantified clone abundances (proportions) per patient sample
  sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
  
  # Clean column names: remove "_sc" suffix if present
  colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))
  
  # Validate that sample IDs exist as columns
  if (ncol(sample_clone_df) == 0) {
    stop("sample_clone_file must contain sample IDs as columns.")
  }
  
  # STEP 5: Extract UMAP coordinates for visualization
  # UMAP (Uniform Manifold Approximation and Projection) provides 2D representation
  # of high-dimensional single-cell transcriptomic data
  umap_coords <- Embeddings(cohort_seurat_obj, "umap")
  
  # STEP 6: Extract PCA coordinates for optimal transport distance calculations
  # PCA captures variance in gene expression, used for computing cell-cell distances
  # Extracting all principal components (PCs) for comprehensive distance matrix
  pca_coords <- Embeddings(cohort_seurat_obj, "pca")
  
  # Assign standardized column names to PCA coordinates
  pca_col_names <- paste0("PC_", seq_len(ncol(pca_coords)))
  colnames(pca_coords) <- pca_col_names
  
  # STEP 7: Create comprehensive data frame combining metadata and embeddings
  # This structure enables Python notebook to access all necessary information:
  # - Cell identifiers and metadata (patient, timepoint, cluster, survivor group)
  # - UMAP coordinates (for visualization)
  # - PCA coordinates (for optimal transport distance calculations)
  combined_df <- data.frame(
    cell_id = rownames(cohort_seurat_obj@meta.data),
    TimePoint = cohort_seurat_obj@meta.data$TimePoint,
    SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup,
    seurat_clusters = cohort_seurat_obj@meta.data$seurat_clusters,
    origin = cohort_seurat_obj@meta.data$origin,
    Patient = cohort_seurat_obj@meta.data$Patient,
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    pca_coords  # All PCA dimensions appended as additional columns
  )
  
  ################################################################################
  # CLONE IDENTIFICATION PIPELINE
  ################################################################################
  
  # CLONE STEP 1: Identify cell barcodes in this subpopulation
  # These are the cells we want to track (e.g., all Effector_CD8 cells)
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # CLONE STEP 2: Map cell barcodes → TCR clones (CDR3β sequences)
  # Looks up CDR3β amino acid sequence for each cell in the subpopulation
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)  # Remove duplicates
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]  # Remove NA values
  
  if (length(clones_from_cells) == 0) {
    message("No clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # CLONE STEP 3: Identify clones with detectable abundance in patient samples
  # This filters for clones that actually have non-zero proportion in quantified repertoire
  sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
  
  # Extract clones with proportion > 0 in any sample from this subpopulation
  clones_from_samples <- c()
  for (sample_id in sample_ids) {
    if (sample_id %in% colnames(sample_clone_df)) {
      # Get clones with non-zero proportion in this sample
      clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
      clones_from_samples <- c(clones_from_samples, clones_in_sample)
    }
  }
  clones_from_samples <- unique(clones_from_samples)
  
  # CLONE STEP 4: Intersect the two clone sets
  # Final clone set = clones found in subpopulation cells AND with detectable proportion
  # This ensures we track only high-confidence, subpopulation-enriched clones
  intersected_clones <- intersect(clones_from_cells, clones_from_samples)
  
  if (length(intersected_clones) == 0) {
    message("No intersected clones found for subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # CLONE STEP 5: Retrieve ALL cells bearing the intersected clones
  # IMPORTANT: This retrieves cells across ALL timepoints and ALL subpopulations
  # Not limited to the focal subpopulation - enables tracking clone migration
  # across different T cell states (e.g., Effector → Memory transitions)
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These cells will be highlighted (colored) in optimal transport visualizations
  cells_to_plot <- all_cells_for_clones
  
  ################################################################################
  # TIMEPOINT PAIR PROCESSING
  ################################################################################
  
  # Identify timepoints present in this cohort
  # Not all cohorts have all timepoints (e.g., some patients lack late timepoints)
  cohort_timepoints <- intersect(all_timepoints, unique(subpop_seurat_obj@meta.data$TimePoint))
  
  # Loop over consecutive timepoint pairs (Pre→C1, C1→C2, C2→C4, etc.)
  # Excludes the last timepoint since it has no subsequent target
  for (i in 1:(length(cohort_timepoints)-1)) {
    source_tp <- cohort_timepoints[i]      # Source timepoint (e.g., Pre)
    target_tp <- cohort_timepoints[i+1]    # Target timepoint (e.g., C1)
    
    # Extract source cells: Highlighted clones at source timepoint
    # These are the "origin" cells for optimal transport calculation
    source_cells <- combined_df %>%
      dplyr::filter(TimePoint == source_tp & cell_id %in% cells_to_plot)
    
    # Extract target cells: Highlighted clones at target timepoint
    # These are the "destination" cells for optimal transport calculation
    target_cells <- combined_df %>%
      dplyr::filter(TimePoint == target_tp & cell_id %in% cells_to_plot)
    
    # Validate that both source and target cells exist
    # Skip this timepoint pair if either is missing
    if (nrow(source_cells) == 0 | nrow(target_cells) == 0) {
      message("No cells found for source timepoint: ", source_tp, " or target timepoint: ", target_tp, " for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
      next
    }
    
    # Extract background cells: All other cells at source timepoint
    # Provides spatial context for visualization (shown in gray)
    gray_cells <- combined_df %>%
      dplyr::filter(TimePoint == source_tp & !(cell_id %in% cells_to_plot))
    
    ################################################################################
    # CSV FILE EXPORT
    ################################################################################
    
    # Create output subfolder organized by source timepoint
    # Structure: base_dir/subpopulation/cohort/Timepoint_X/
    output_subfolder <- file.path(output_folder, paste0("Timepoint_", source_tp))
    if (!dir.exists(output_subfolder)) {
      dir.create(output_subfolder, recursive = TRUE)
    }
    
    # Save three CSV files for this timepoint pair:
    # 1. source_cells_new.csv: Highlighted clones at source timepoint (with UMAP + PCA)
    # 2. target_cells_new.csv: Highlighted clones at target timepoint (with UMAP + PCA)
    # 3. gray_cells_new.csv: Background cells at source timepoint (with UMAP + PCA)
    write.csv(source_cells, file = file.path(output_subfolder, "source_cells_new.csv"), row.names = FALSE)
    write.csv(target_cells, file = file.path(output_subfolder, "target_cells_new.csv"), row.names = FALSE)
    write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells_new.csv"), row.names = FALSE)
  }
  
  ################################################################################
  # SPECIAL HANDLING: LAST TIMEPOINT
  ################################################################################
  
  # Last timepoint has no subsequent target, so only save source and background cells
  # This enables visualization of the final timepoint distribution
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
  
  # Save only source and gray cells (no target cells for last timepoint)
  write.csv(source_cells, file = file.path(output_subfolder, "source_cells_new.csv"), row.names = FALSE)
  write.csv(gray_cells, file = file.path(output_subfolder, "gray_cells_new.csv"), row.names = FALSE)
}

################################################################################
# MAIN EXECUTION: ITERATION OVER SUBPOPULATIONS AND COHORTS
################################################################################

# Define cohort-to-patient mapping for iteration
cohorts <- list(
  "control" = control_group,
  "short_term" = short_term_survivor_group,
  "long_term" = long_term_survivor_group
)

# Nested loop structure:
# Outer loop: Iterate over CD8+ T cell subpopulations
# Inner loop: Iterate over patient survival cohorts
#
# This generates data for all combinations:
# 8 subpopulations × 3 cohorts = 24 total analysis units
#
# Each analysis unit produces CSV files for all available timepoint pairs

for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]  # Get Seurat cluster IDs for this subpopulation
  
  for (cohort_name in names(cohorts)) {
    patient_ids <- cohorts[[cohort_name]]  # Get patient IDs for this cohort
    
    # Create output directory structure: base_dir/subpopulation/cohort/
    output_folder <- file.path(base_output_dir, subpop_name, cohort_name)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Execute data extraction and CSV export for this subpopulation-cohort combination
    extractDataForOptimalTransport(subpop_name, clusters, cohort_name, patient_ids, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }
}

################################################################################
# SCRIPT COMPLETION
################################################################################
#
# After running this script, the output directory structure contains organized
# CSV files ready for Python-based optimal transport analysis.
#
# NEXT STEP: Run t_cell_optimal_transport.ipynb
# The Jupyter notebook will:
# 1. Load source_cells_new.csv, target_cells_new.csv, gray_cells_new.csv
# 2. Compute optimal transport plan using Earth Mover's Distance
# 3. Generate arrow plots showing individual cell fate transitions
# 4. Create density shift visualizations colored by survival cohort
# 5. Export predecessor mapping for downstream clone tracking analysis
#
# The combination of this R script (data preparation) and the Python notebook
# (optimal transport computation) provides a complete pipeline for analyzing
# T cell clone dynamics and spatial reorganization during immunotherapy treatment.
################################################################################
