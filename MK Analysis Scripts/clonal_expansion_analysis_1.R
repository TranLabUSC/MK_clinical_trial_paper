################################################################################
# CLONAL EXPANSION ANALYSIS - STEP 1: DATA PREPARATION
################################################################################
#
# PURPOSE:
#   This is the FIRST step (1/4) in the clonal expansion analysis pipeline.
#   It prepares clonotyping data tables (absolute counts and proportions) for
#   all T cell subpopulations by processing TCR sequencing data.
#
# ANALYSIS PIPELINE:
#   Step 1 (THIS SCRIPT): Prepare clonotype count and proportion tables
#   Step 2: diversity_calculation_2.R - Calculate diversity metrics
#   Step 3: clonal_expansion_calculation_3.R - Calculate expansion metrics
#   Step 4: clonal_expansion_plots_4.R - Generate visualizations (Figure 6d)
#
# MANUSCRIPT FIGURE:
#   Figure 6d: Clonal expansion visualization (generated after completing all 4 steps)
#
# WORKFLOW:
#   1. Load T cell Seurat object and filter for UF site only
#   2. Load TCR contig annotations (CDR3 sequences)
#   3. Define T cell subpopulation mappings (cell types to cluster IDs)
#   4. For each T cell subpopulation:
#      a. Extract cells belonging to that subpopulation
#      b. Count clone (CDR3) frequency per sample
#      c. Create absolute count table (clones × samples)
#      d. Create proportional table (normalized within each sample)
#      e. Save both tables for downstream analysis
#
# KEY CONCEPTS:
#   - Clone: Group of T cells sharing the same CDR3 amino acid sequence
#   - Clonotype: Unique TCR sequence (we use CDR3β chain)
#   - Clonal Expansion: Presence of the same clone across multiple cells
#
# INPUT FILES:
#   - T cell Seurat object with cluster annotations
#   - TCR VDJ contig annotations (filtered_contig_annotations.csv)
#
# OUTPUT FILES (per cell type):
#   - clonotype_df_absolute_*.csv: Raw clone counts (rows=clones, cols=samples)
#   - clonotype_df_proportion_*.csv: Proportions (normalized within each sample)
#
################################################################################

# Load required libraries
library(Seurat)        # Single-cell analysis
library(dplyr)         # Data manipulation
library(data.table)    # Fast data handling
library(ggplot2)       # Plotting
library(survival)      # Survival analysis
library(broom)         # Tidy model outputs
library(ggfortify)     # For autoplot
library("survminer")   # Survival plots
library("Rcpp")        # C++ integration
library(cowplot)       # Plot composition
library(tidyr)         # Data tidying
library(rlang)         # Programming tools

################################################################################
# SECTION 1: LOAD AND PREPARE DATA
################################################################################

# Load T cell Seurat object (from annotation step)
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# Filter for UF site only (WUSTL TCR data is poor quality and excluded)
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata <- seurat_object_t_cells@meta.data

# Load TCR contig annotations from 10x VDJ data
# This file contains barcode-to-CDR3 sequence mappings
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)

# Filter for TRB chain only (β chain is more stable and reliable than α chain)
clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]

################################################################################
# SECTION 2: DEFINE T CELL SUBPOPULATION MAPPINGS
################################################################################
# Map cell types to cluster IDs based on marker gene expression from annotation
# Each cluster number corresponds to a specific T cell subpopulation
# Some cell types span multiple clusters (e.g., Naive_CD4 includes clusters 6, 9, 18)

mapping <- list(
  "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
  
  # CD4+ T cell subpopulations
  "Activated_CD4" = c(0),
  "Active_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  
  # CD8+ T cell subpopulations
  "Effector_CD8" = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "GZMK_Effector_Memory_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17),
  
  # Other T cell types
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  
  # Combined categories for analysis
  "Effector_and_Central_Memory_CD8" = c(1, 12),
  "Effector_Memory_Precursor_and_Central_Memory_CD8" = c(2, 10, 12),
  "Memory_Precursor_and_Central_Memory_CD8" = c(2, 12),
  "Effector_Memory_and_Central_Memory_CD8" = c(10, 12),
  "All_CD8" = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)

# Convert mapping list to a long-format dataframe
# Each row represents one cluster belonging to one cell type
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

################################################################################
# SECTION 3: PROCESS EACH T CELL SUBPOPULATION
################################################################################
# For each cell type, create clonotype frequency tables

for (celltype in unique(celltype_to_cluster$celltype)){
  print(celltype)
  
  # Get cluster IDs for this cell type
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Subset Seurat object for cells in these clusters
  seurat_object_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
  
  # Get cell barcodes for this cell type
  barcodes <- rownames(seurat_object_subset@meta.data)
  
  # Filter TCR data to include only cells from this cell type
  clone_cell_barcode_df_subset <- clone_cell_barcode_df[clone_cell_barcode_df$barcode %in% barcodes,]

  ################################################################################
  # SUBSECTION 3.1: BUILD CLONE FREQUENCY DATA STRUCTURE
  ################################################################################
  # Create nested list structure: sample -> CDR3 sequence -> count
  # Also track barcodes for each sample
  
  final_list <- list()
  
  for (i in 1:nrow(clone_cell_barcode_df_subset)) {
    current_row <- clone_cell_barcode_df_subset[i, ]
    barcode <- current_row$barcode              # Cell barcode
    aa_sequence <- current_row$cdr3             # CDR3 amino acid sequence (clone identifier)
    origin <- current_row$origin                # Sample ID (e.g., patient_timepoint)
    
    # Initialize sample in the list if not present
    if (origin %in% names(final_list)){
      print("sample exists")
    } else {
      final_list[[origin]] <- list()
      final_list[[origin]][["barcodes"]] <- c()
    }
    
    # Initialize CDR3 sequence counter if not present
    if (aa_sequence %in% names(final_list[[origin]])){
      print("AA sequence exists")
    } else {
      final_list[[origin]][[aa_sequence]] <- 0
    }

    # Update the data structure
    final_list[[origin]][["barcodes"]] <- c(final_list[[origin]][["barcodes"]], barcode)
    final_list[[origin]][[aa_sequence]] <- final_list[[origin]][[aa_sequence]] + 1
  }

  ################################################################################
  # SUBSECTION 3.2: CONVERT TO DATAFRAME FORMAT
  ################################################################################
  # Function to convert nested list to a dataframe
  # Output: rows = CDR3 sequences (clones), columns = samples
  convert_to_df <- function(final_list) {
    clone_df <- data.frame()
    
    for (sample_name in names(final_list)) {
      print(sample_name)
      sample_list <- final_list[[sample_name]]
      
      # Remove barcodes list (we only need counts)
      sample_list$barcodes <- NULL
      
      # Convert to dataframe: rows = CDR3, col = count
      temp_df <- as.data.frame(as.list(sample_list), stringsAsFactors = FALSE)
      temp_df <- as.data.frame(t(temp_df))
      colnames(temp_df) <- sample_name
      
      # Merge with main dataframe (outer join to keep all clones)
      clone_df <- merge(clone_df, temp_df, by="row.names", all = TRUE)
      rownames(clone_df) <- clone_df$Row.names
      clone_df$Row.names <- NULL
    }
    return(clone_df)
  }

  ################################################################################
  # SUBSECTION 3.3: CREATE AND SAVE ABSOLUTE COUNT TABLE
  ################################################################################
  # Convert to dataframe and fill missing values with 0
  clonotype_df_absolute <- convert_to_df(final_list)
  clonotype_df_absolute[is.na(clonotype_df_absolute)] <- 0
  clonotype_df_absolute$CDR3.aa <- rownames(clonotype_df_absolute)
  
  # Save absolute counts
  write.csv(clonotype_df_absolute, 
           paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/clonotype_df_absolute_", celltype, "_T_cells.csv"))

  ################################################################################
  # SUBSECTION 3.4: FILTER LOW-QUALITY SAMPLES
  ################################################################################
  # Remove samples with <1 cell (likely technical artifacts or contamination)
  
  numeric_cols <- sapply(clonotype_df_absolute, is.numeric)
  print(numeric_cols)
  print(colnames(clonotype_df_absolute)[numeric_cols])
  
  # Calculate sum for each sample column
  # Remove columns where total cell count is less than 1
  cols_to_remove <- sapply(clonotype_df_absolute[, numeric_cols, drop = FALSE], sum) < 1
  clonotype_df_absolute <- clonotype_df_absolute[, !cols_to_remove]

  ################################################################################
  # SUBSECTION 3.5: CREATE AND SAVE PROPORTIONAL TABLE
  ################################################################################
  # Convert absolute counts to proportions (normalize within each sample)
  # This allows fair comparison between samples with different total cell counts
  
  # Function to calculate proportions for all numeric columns
  calculate_proportions_all_numeric <- function(df) {
    numeric_cols <- sapply(df, is.numeric)

    # Copy the original dataframe to keep non-numeric data intact
    new_df <- df

    # Divide each sample column by its sum (normalize to proportions)
    # Each sample column will now sum to 1.0
    new_df[, numeric_cols] <- sapply(df[, numeric_cols, drop = FALSE], function(x) x / sum(x))

    # Verify that sums are now 1.0 (quality check)
    print(sapply(new_df[, numeric_cols, drop = FALSE], sum))
    return(new_df)
  }

  # Calculate proportions for all numeric columns
  proportion_df <- calculate_proportions_all_numeric(clonotype_df_absolute)
  
  # Save proportional table
  write.csv(proportion_df, 
           paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/clonotype_df_proportion_", celltype, "_T_cells.csv"))
}

################################################################################
# END OF SCRIPT
################################################################################
# 
# OUTPUT SUMMARY:
# For each T cell subpopulation, two CSV files are created:
#   1. clonotype_df_absolute_*.csv  - Raw clone counts per sample
#   2. clonotype_df_proportion_*.csv - Clone frequencies (proportions) per sample
#
# These tables serve as input for:
#   - Step 2: diversity_calculation_2.R (diversity metrics calculation)
#   - Step 3: clonal_expansion_calculation_3.R (expansion metrics)
#   - Step 4: clonal_expansion_plots_4.R (Figure 6d visualization)
#
################################################################################