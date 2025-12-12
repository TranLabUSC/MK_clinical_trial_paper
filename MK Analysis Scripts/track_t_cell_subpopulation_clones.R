
################################################################################
# T Cell Clone Tracking Across Timepoints - Figure 7a Generation
################################################################################
#
# PURPOSE:
# This script tracks T cell clones across timepoints to generate Figure 7a,
# which visualizes CD8+ T cell clone distribution using density contours.
# The analysis identifies subpopulation-specific clones (defined by TCR CDR3β
# sequences) and tracks their spatial distribution in UMAP space over time.
#
# FIGURE 7A DESCRIPTION:
# - Displays CD8+ T cell clone density across timepoints: C1, C2, C4, C6, C9, C17, C34
# - Uses density contour plots to visualize clone spatial distribution
# - Gray background shows all cells in the timepoint
# - Colored contours show tracked clones from specific subpopulations
#
# PATIENT COHORTS (Color-Coded):
# - Yellow: Control group (NLS+PEM = Non-Long Survivors + PEM)
# - Blue: Short-term survivors (LTT+PEM < mOS, median Overall Survival)
# - Red: Long-term survivors (LTT+PEM > mOS)
#
# CLONE TRACKING METHODOLOGY:
# 1. TCR-based clonotyping: Identifies clones by CDR3β amino acid sequence
# 2. Subpopulation specificity: Filters for clones present in target subpopulation
# 3. Temporal tracking: Follows clone spatial distribution across all timepoints
# 4. Cohort comparison: Color-codes clones by patient survival group
# 5. Spatial visualization: Uses UMAP embeddings to show clone localization
#
# BIOLOGICAL INTERPRETATION:
# - Clone expansion/contraction: Changes in density indicate clonal dynamics
# - Spatial migration: Shifts in contour location suggest differentiation
# - Cohort differences: Color patterns reveal survival-associated clone behaviors
# - Temporal evolution: Sequential timepoints show longitudinal clone trajectories
#
# OUTPUT:
# - Contour density plots (Figure 7a style): Density contours on gray background
# - Filled density plots: Enhanced visualization with polygon fills
# - Separate plots for each T cell subpopulation across all timepoints
#
################################################################################

##############################################################################################################################################################################################################################################################
# SETUP: Load Libraries and Define Subpopulations
##############################################################################################################################################################################################################################################################

# Load required libraries
library(Seurat)    # Single-cell RNA-seq analysis
library(ggplot2)   # Visualization
library(cowplot)   # Plot arrangement
library(dplyr)     # Data manipulation

##############################################################################################################################################################################################################################################################
# T Cell Subpopulation Definitions
##############################################################################################################################################################################################################################################################

# Define cluster-to-subpopulation mapping
# Each subpopulation is defined by one or more Seurat cluster IDs
# These mappings group functionally related T cell states based on marker expression
mapping <- list(
  "Activated_CD4" = c(0),
  "Effector_CD8" = c(1),
  "Effector_Memory_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  "Active_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  "Memory_CD8" = c(8),
  "Anergic_CD8" = c(10),
  "Naive_CD8" = c(12),
  "Hyperactivated_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17),
  "CD8" = c(1, 2, 13, 14, 16, 17)
)

# Define all possible timepoints
all_timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")
all_patients <- c(1, 4, 8, 9, 7, 10, 12, 14, 18, 2, 3, 5, 13, 19, 20, 21)

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# Filter for UF site samples and selected patients
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF" & Patient %in% all_patients)

# Remove cluster 13 (excluded from analysis)
full_seurat_obj <- subset(full_seurat_obj, subset = seurat_clusters != 13)

##############################################################################################################################################################################################################################################################
# Define Patient Cohorts by Survival Outcomes
##############################################################################################################################################################################################################################################################

# Clinical cohort definitions based on treatment response and survival
# Control group: Non-Long Survivors + PEM (NLS+PEM) - visualized in YELLOW
control_group <- c(1, 4, 8, 9)

# Short-term survivors: LTT+PEM < median OS - visualized in BLUE
short_term_survivor_group <- c(7, 10, 12, 14, 18)

# Long-term survivors: LTT+PEM > median OS - visualized in RED
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character for metadata matching
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Add SurvivorGroup metadata to Seurat object for color-coding
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term",
         ifelse(full_seurat_obj$Patient %in% control_group, "control", NA))
)

##############################################################################################################################################################################################################################################################
# FUNCTION 1: createCloneTrackingUMAPs - Contour Density Visualization (Figure 7a Style)
##############################################################################################################################################################################################################################################################
#
# PURPOSE:
# Tracks T cell clones specific to a subpopulation across timepoints and visualizes
# their spatial distribution using density contour plots (Figure 7a format).
#
# CLONE IDENTIFICATION WORKFLOW:
# Step 1: Extract all cells belonging to the target subpopulation (by cluster ID)
# Step 2: Map cell barcodes to TCR clone IDs (CDR3β amino acid sequences)
# Step 3: Identify clones present in subpopulation samples (proportion > 0)
# Step 4: Take intersection to find subpopulation-specific clones
# Step 5: Track ALL cells bearing these clones across all timepoints
# Step 6: Color-code cells by patient survival group (red/blue/yellow)
#
# VISUALIZATION STRATEGY:
# - Gray background: All T cells in each timepoint (spatial context)
# - Colored contours: Density of tracked clones from subpopulation
# - Temporal series: One panel per timepoint showing clone evolution
# - Cohort comparison: Color indicates survival group affiliation
#
# GENERATES FIGURE 7A:
# CD8+ T cell clone distribution across treatment timepoints with density contours
#
# PARAMETERS:
# @param subpop_name: Name of T cell subpopulation (e.g., "Effector_CD8")
# @param clusters: Vector of cluster IDs defining the subpopulation
# @param full_seurat_obj: Complete Seurat object with all T cells
# @param clone_cell_barcode_file: CSV mapping cell barcodes to TCR CDR3β sequences
# @param sample_clone_file: CSV with clone proportions per sample
# @param output_folder: Directory to save output plots
#
# RETURNS:
# None (saves PDF files):
# - <subpop>_clone_tracking_umap_*.pdf: Colored cell UMAP plots
# - <subpop>_density_plots_with_background_*.pdf: Density contour plots (Figure 7a)
#
################################################################################

createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder) {

  # Use complete dataset to include all patients and timepoints
  cohort_seurat_obj <- full_seurat_obj

  # Identify timepoints available in the dataset
  cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))

  # Extract cells belonging to target subpopulation (by cluster IDs)
  subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)

  # Validate subpopulation has cells
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name)
    return()
  }

  # CLONE BARCODE MAPPING: Load TCR data linking cell barcodes to CDR3β sequences
  clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
  # Filter for TRB chain (β chain defines clonotype)
  clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]

  # Validate required columns exist
  if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
    stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
  }

  # CLONE PROPORTION DATA: Load clone abundance information per sample
  sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
  colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))

  # Validate sample IDs are present
  if (ncol(sample_clone_df) == 0) {
    stop("sample_clone_file must contain sample IDs as columns.")
  }

  # UMAP COORDINATES: Extract spatial embeddings for visualization
  umap_coords <- Embeddings(cohort_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = cohort_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(cohort_seurat_obj@meta.data),
    SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup
  )

  # Fix axes across all timepoints for consistent spatial reference
  x_limits <- range(umap_df$UMAP_1)
  y_limits <- range(umap_df$UMAP_2)

  ##############################################################################
  # CLONE IDENTIFICATION PROCEDURE (Subpopulation-Specific Filtering)
  ##############################################################################
  
  # STEP 1: Extract cell barcodes from target subpopulation
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    message("No cells found in subpopulation: ", subpop_name)
    return()
  }

  # STEP 2: Map cell barcodes to TCR clone IDs (CDR3β sequences)
  # Extract CDR3β sequences for cells in subpopulation
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]

  if (length(clones_from_cells) == 0) {
    message("No clones found for subpopulation: ", subpop_name)
    return()
  }

  # STEP 3: Identify clones present in subpopulation samples (proportion > 0)
  # This ensures clones are actually detected in the quantified repertoire
  sample_ids <- unique(subpop_seurat_obj@meta.data$origin)

  clones_from_samples <- c()
  for (sample_id in sample_ids) {
    if (sample_id %in% colnames(sample_clone_df)) {
      clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
      clones_from_samples <- c(clones_from_samples, clones_in_sample)
    }
  }
  clones_from_samples <- unique(clones_from_samples)

  # STEP 4: Intersection = Subpopulation-specific clones
  # These clones are both (1) found in subpopulation cells and (2) detected in samples
  intersected_clones <- intersect(clones_from_cells, clones_from_samples)

  if (length(intersected_clones) == 0) {
    message("No intersected clones found for subpopulation: ", subpop_name)
    return()
  }

  # STEP 5: Track ALL cells bearing these clones across ALL timepoints
  # This reveals how subpopulation-specific clones are distributed temporally
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]

  # These are the cells to visualize (clone-bearing cells)
  cells_to_plot <- all_cells_for_clones

  ##############################################################################
  # VISUALIZATION GENERATION (Temporal Tracking Across Timepoints)
  ##############################################################################
  
  # Initialize plot storage
  umap_plot_list <- list()
  density_plot_list <- list()
  plot_counter <- 1

  # Generate plots for each timepoint
  for (col_tp in cohort_timepoints) {
    # Subset to current timepoint
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)

    # Prepare plot dataframe
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"

    # COLOR CODING: Assign colors based on survival group
    # Gray = background cells, Colored = tracked clones colored by patient cohort
    plot_df$color <- "gray"
    plot_df$color[plot_df$highlight == "highlighted"] <- plot_df$SurvivorGroup[plot_df$highlight == "highlighted"]
    plot_df$color[is.na(plot_df$color)] <- "unknown"

    # Map survival groups to visualization colors
    # Control (NLS+PEM) = Yellow, Short-term = Blue, Long-term = Red
    color_mapping <- c("gray" = "gray", "long_term" = "red", "short_term" = "blue", "control" = "yellow", "unknown" = "black")
    plot_df$color <- color_mapping[plot_df$color]

    # UMAP PLOT: Scatter plot with clone cells colored by cohort
    p_umap <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = color), size = 0.5, alpha = ifelse(plot_df$color == "gray", 0.5, 1)) +
      scale_color_identity() +
      coord_fixed() +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

    # Store the UMAP plot
    umap_plot_list[[plot_counter]] <- p_umap + coord_flip()

    # DENSITY CONTOUR PLOT (Figure 7a Format): Shows clone spatial distribution
    # Remove non-clone cells to calculate density only for tracked clones
    density_df <- plot_df[plot_df$color != "gray", ]

    # Layer 1: Gray background showing all cells (spatial context)
    p_density <- ggplot() +
      geom_point(data = plot_df[plot_df$color == "gray", ], aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5)

    # Layer 2: Add density contours for tracked clones (colored by cohort)
    # Contours reveal spatial concentration of subpopulation-specific clones
    if (nrow(density_df) > 0) {
      p_density <- p_density +
        geom_density_2d(data = density_df, aes(x = UMAP_1, y = UMAP_2, color = color), size = 0.7) +
        scale_color_identity()
    }

    p_density <- p_density +
      coord_fixed() +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

    # Store the density plot with flip for correct orientation
    density_plot_list[[plot_counter]] <- p_density + coord_flip()
    plot_counter <- plot_counter + 1
  }

  message("Total number of UMAPs and density plots: ", plot_counter - 1)

  # Validate plot generation
  if (length(umap_plot_list) == 0 || length(density_plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }

  ##############################################################################
  # SAVE OUTPUTS: Arrange plots in temporal grid and save
  ##############################################################################
  
  # Arrange plots: 1 row × N columns (one column per timepoint)
  n_rows <- 1
  n_cols <- length(umap_plot_list)

  # Create UMAP grid
  umap_grid_combined <- plot_grid(plotlist = umap_plot_list, ncol = n_cols)

  # Save UMAP plots (colored cells showing cohort affiliation)
  output_file_umap <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf"))

  ggsave(filename = output_file_umap, plot = umap_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)

  message("Saved UMAP grid for subpopulation: ", subpop_name)

  # Create density contour grid
  density_grid_combined <- plot_grid(plotlist = density_plot_list, ncol = n_cols)

  # Save density contour plots (FIGURE 7A FORMAT)
  # These plots show clone spatial distribution with contours on gray background
  output_file_density <- file.path(output_folder, paste0(subpop_name, "_density_plots_with_background_all_patients_all_cells_fixed_axes.pdf"))

  ggsave(filename = output_file_density, plot = density_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)

  message("Saved density plot grid for subpopulation: ", subpop_name)
}

##############################################################################################################################################################################################################################################################
# FUNCTION 2: createCloneTrackingUMAPs - Filled Density Visualization (Enhanced Version)
##############################################################################################################################################################################################################################################################
#
# PURPOSE:
# Alternative visualization of clone tracking using filled density polygons instead
# of contour lines. Provides more prominent visual representation of clone distributions.
#
# METHODOLOGY:
# Identical clone identification workflow as Function 1, but uses filled density
# polygons (stat_density_2d with geom="polygon") for enhanced visualization.
#
# DIFFERENCE FROM FUNCTION 1:
# - Function 1: Uses geom_density_2d for contour lines (Figure 7a style)
# - Function 2: Uses stat_density_2d with polygon fill for filled regions
#
# All other aspects (clone tracking, color coding, temporal analysis) are identical.
#
################################################################################

createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder) {
  
  # Use full_seurat_obj to include all patients
  cohort_seurat_obj <- full_seurat_obj
  
  # Define all timepoints present in the data
  cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))
  
  # Create subpop_seurat_obj
  subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name)
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
  
  # Prepare umap_df
  umap_coords <- Embeddings(cohort_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = cohort_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(cohort_seurat_obj@meta.data),
    SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup
  )
  
  # Get overall x and y limits
  x_limits <- range(umap_df$UMAP_1)
  y_limits <- range(umap_df$UMAP_2)
  
  # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    # No cells in subpopulation
    message("No cells found in subpopulation: ", subpop_name)
    return()
  }
  
  # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name)
    return()
  }
  
  # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
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
  
  # Step 5: Take intersection of the two
  intersected_clones <- intersect(clones_from_cells, clones_from_samples)
  
  if (length(intersected_clones) == 0) {
    # No intersected clones
    message("No intersected clones found for subpopulation: ", subpop_name)
    return()
  }
  
  # Step 6: Identify all cell barcodes corresponding to the intersected clones
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # Cells to plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot lists
  umap_plot_list <- list()
  density_plot_list <- list()
  plot_counter <- 1
  
  # Loop over timepoints
  for (col_tp in cohort_timepoints) {
    # Cells in current timepoint
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Prepare data frames
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"
    
    # Assign colors based on SurvivorGroup
    plot_df$color <- "gray"
    plot_df$color[plot_df$highlight == "highlighted"] <- plot_df$SurvivorGroup[plot_df$highlight == "highlighted"]
    plot_df$color[is.na(plot_df$color)] <- "unknown"
    
    # Map colors to actual colors
    color_mapping <- c("gray" = "gray", "long_term" = "red", "short_term" = "blue", "control" = "yellow", "unknown" = "black")
    plot_df$color <- color_mapping[plot_df$color]
    
    # UMAP PLOT: Same as Function 1
    p_umap <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = color), size = 0.5, alpha = ifelse(plot_df$color == "gray", 0.5, 1)) +
      scale_color_identity() +
      coord_fixed() +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))
    
    # Store the UMAP plot
    umap_plot_list[[plot_counter]] <- p_umap + coord_flip()
    
    # FILLED DENSITY PLOT: Enhanced visualization with polygon fills
    # Isolate clone cells for density calculation
    density_df <- plot_df[plot_df$color != "gray", ]
    
    # Layer 1: Gray background (all cells)
    p_density <- ggplot() +
      geom_point(data = plot_df[plot_df$color == "gray", ],
                 aes(x = UMAP_1, y = UMAP_2),
                 color = "gray", size = 0.5, alpha = 0.5)
    
    # Layer 2: Add filled density polygons for each survivor group
    # Creates regions of color showing clone concentration
    groups_present <- unique(density_df$SurvivorGroup)
    groups_present <- groups_present[!is.na(groups_present)]
    
    for (grp in groups_present) {
      grp_data <- density_df[density_df$SurvivorGroup == grp, ]
      grp_color <- color_mapping[[grp]]

      # stat_density_2d with polygon geom creates filled regions
      if (nrow(grp_data) > 2) {  # Minimum 3 points required for density estimation
        p_density <- p_density +
          stat_density_2d(data = grp_data,
                          aes(x = UMAP_1, y = UMAP_2, fill = ..level..),
                          geom = "polygon",
                          bins = 10,
                          alpha = 0.5,
                          color = NA,
                          fill = grp_color,
                          contour = TRUE,
                          inherit.aes = FALSE)
      }
    }
    
    p_density <- p_density +
      coord_fixed() +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")
      )
    
    # Store the filled density plot
    density_plot_list[[plot_counter]] <- p_density + coord_flip()
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs and density plots: ", plot_counter - 1)
  
  # Validate plot generation
  if (length(umap_plot_list) == 0 || length(density_plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }
  
  ##############################################################################
  # SAVE OUTPUTS: Arrange and save filled density plots
  ##############################################################################
  
  # Arrange plots in temporal grid (1 row × N columns)
  n_rows <- 1
  n_cols <- length(umap_plot_list)
  
  # Create UMAP grid
  umap_grid_combined <- plot_grid(plotlist = umap_plot_list, ncol = n_cols)
  
  # Save UMAP plots
  output_file_umap <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf"))
  
  ggsave(filename = output_file_umap, plot = umap_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name)
  
  # Create filled density grid
  density_grid_combined <- plot_grid(plotlist = density_plot_list, ncol = n_cols)
  
  # Save filled density plots (polygon-filled regions showing clone distribution)
  output_file_density <- file.path(output_folder, paste0(subpop_name, "_filled_density_plots_with_background_all_patients_all_cells_fixed_axes.pdf"))
  
  ggsave(filename = output_file_density, plot = density_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved density plot grid for subpopulation: ", subpop_name)
}

##############################################################################################################################################################################################################################################################
# MAIN EXECUTION: Process All T Cell Subpopulations
##############################################################################################################################################################################################################################################################

# Input data file paths
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/temp"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"


##############################################################################################################################################################################################################################################################
# Iterate Over All T Cell Subpopulations
##############################################################################################################################################################################################################################################################
#
# LOOP STRUCTURE:
# For each defined subpopulation in the mapping list:
# 1. Extract cluster IDs defining the subpopulation
# 2. Create dedicated output directory
# 3. Call createCloneTrackingUMAPs function (generates both visualization types)
# 4. Handle errors gracefully to continue processing other subpopulations
#
# BIOLOGICAL INSIGHT:
# Different subpopulations show distinct temporal clone dynamics:
# - Effector CD8: Rapid expansion at early timepoints, contraction later
# - Memory CD8: Persistent clones maintained across multiple timepoints
# - Exhausted T: Stable or increasing presence indicating chronic stimulation
# - Naive CD4/CD8: Minimal clonal expansion, diverse TCR repertoire
#
# COHORT COMPARISON:
# Clone tracking across cohorts reveals survival-associated patterns:
# - Long-term survivors (RED): May show sustained effector/memory clones
# - Short-term survivors (BLUE): Different clone dynamics or earlier exhaustion
# - Control group (YELLOW): Baseline clone behavior without long-term treatment
#
# OUTPUT FILES PER SUBPOPULATION:
# 1. *_clone_tracking_umap_*.pdf: Colored cell scatter plots
# 2. *_density_plots_with_background_*.pdf: Contour density plots (Figure 7a)
# 3. *_filled_density_plots_*.pdf: Filled polygon density plots
#
################################################################################

for (subpop_name in names(mapping)) {
  # Get cluster IDs for this subpopulation
  clusters <- mapping[[subpop_name]]
  
  # Create subpopulation-specific output directory
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Generate clone tracking visualizations
  # Both contour and filled density plots are created by the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }, error = function(e) {
    # Log errors but continue with remaining subpopulations
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}

################################################################################
# END OF SCRIPT
#
# EXPECTED OUTPUTS:
# - Figure 7a: CD8+ T cell clone distribution across timepoints with density contours
# - Additional subpopulation-specific tracking plots for all 14 T cell subpopulations
# - Temporal visualization of clone spatial dynamics across treatment course
# - Cohort comparison showing survival-associated clone behaviors
#
################################################################################