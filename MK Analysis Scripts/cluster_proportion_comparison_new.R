#########################################################################################################################
#########################################################################################################################
##
## CLUSTER PROPORTION COMPARISON: Treatment Cohort Analysis
##
## PURPOSE: Compare cluster proportions between MK-3475 Alone and MK-3475 + MLA treatment arms
##          across three timepoints (Pre-treatment, Cycle 1, Cycle 2) to identify differential
##          cell population dynamics in response to combination therapy.
##
## BIOLOGICAL CONTEXT:
## - Cluster proportion = (Number of cells in cluster) / (Total cells at timepoint for patient)
## - Changes in cluster proportions reveal shifts in immune cell composition over treatment
## - Higher proportions may indicate cell expansion, activation, or recruitment
## - Lower proportions may indicate differentiation, exhaustion, or migration
##
## OUTPUTS:
## - Supplemental Table S3: Comparison results for 5 major cell types (Classical Monocytes,
##                          Non-Classical Monocytes, cDC, pDC, NK cells)
## - Supplemental Table S4: Comparison results for 13 T cell subpopulations
##
## STATISTICAL APPROACH:
## - Wilcoxon rank-sum test: Compares proportion changes between treatment arms
## - Log2 fold change: Median change across patients (direction and magnitude)
## - Significance levels: * p<0.05, ** p<0.01, *** p<0.001
##
## STRUCTURE:
## - SECTION 1 (Lines 1-303): Analysis of 5 major cell types → Supplemental Table S3
## - SECTION 2 (Lines 311-659): Analysis of 13 T cell subpopulations → Supplemental Table S4
##
#########################################################################################################################
#########################################################################################################################

#########################################################################################################################
## SECTION 1: MAJOR CELL TYPES ANALYSIS (Generates Supplemental Table S3)
##
## Analyzes cluster proportion changes for 5 major immune cell populations:
## - Classical Monocytes (clusters 1, 3, 16, 12)
## - Non-Classical Monocytes (clusters 9, 33)
## - Conventional Dendritic Cells - cDC (cluster 31)
## - Plasmacytoid Dendritic Cells - pDC (cluster 36)
## - Natural Killer cells - NK (cluster 0)
##
## DATA SOURCE: MK_Cells_seurat_obj.RDS (all cells from UF site)
## COMPARISON: MK-3475 Alone vs MK-3475 + MLA across timepoint pairs (Pre-C1, Pre-C2, C1-C2)
#########################################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# =============================================================================================================
# FUNCTION: get_cluster_proportion
# =============================================================================================================
# PURPOSE: Calculate the proportion of cells in specified cluster(s) for each patient at a given timepoint
#
# BIOLOGICAL INTERPRETATION:
# - Cluster proportion represents the fraction of a patient's total cells that belong to specific cluster(s)
# - Formula: Proportion = (Cells in cluster) / (Total cells at timepoint)
# - Values range from 0 to 1 (or 0% to 100%)
# - Example: If patient has 1000 cells at Pre timepoint and 150 are Classical Monocytes,
#            then Classical Monocyte proportion = 150/1000 = 0.15 (15%)
#
# PARAMETERS:
# @param seurat_metadata: Seurat object metadata containing cell-level information
# @param patient_col: Column name containing patient identifiers
# @param timepoint_col: Column name containing timepoint information (Pre, C1, C2)
# @param cluster_col: Column name containing cluster assignments (seurat_clusters)
# @param clusters_of_interest: Vector of cluster numbers to calculate proportion for
# @param timepoint_of_interest: Specific timepoint to analyze (e.g., "Pre", "C1", "C2")
#
# RETURNS:
# @return data.frame with columns:
#         - patient_id: Patient identifier
#         - Cluster_Proportion_[timepoint]: Calculated proportion value for each patient
#
# EXAMPLE:
# If analyzing Classical Monocytes (clusters 1,3,16,12) at Pre timepoint:
#   result will have columns: patient_id, Cluster_Proportion_Pre
# =============================================================================================================
get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for cells in the cluster(s) of interest at the specified timepoint
  cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest &
                                              seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  
  # Filter metadata for ALL cells at the specified timepoint (denominator for proportion)
  timepoint_data <- seurat_metadata[seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  
  # Group filtered data by patient to calculate per-patient proportions
  cluster_timpoint_grouped_data <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
  timpoint_grouped_data <- split(timepoint_data, timepoint_data[[patient_col]])
  
  # Initialize vectors to store results for each patient
  subject_ids <- c()
  cluster_proportions <- c()
  
  # Calculate proportion for each patient
  for (key in names(cluster_timpoint_grouped_data)) {
    patient_cluster_timepoint_data <- cluster_timpoint_grouped_data[[key]]
    patient_timepoint_data <- timpoint_grouped_data[[key]]
    
    # Count cells: numerator (cluster cells) and denominator (total cells)
    patient_cluster_timepoint_cells <- nrow(patient_cluster_timepoint_data)
    patient_timepoint_cells <- nrow(patient_timepoint_data)
    
    # Only calculate proportion if patient has cells in the cluster
    if (patient_cluster_timepoint_cells > 0) {
      # Calculate cluster proportion: (cluster cells) / (total cells)
      cluster_proportion <- patient_cluster_timepoint_cells / patient_timepoint_cells
      
      # Store patient ID and their calculated proportion
      subject_ids <- c(subject_ids, key)
      cluster_proportions <- c(cluster_proportions, cluster_proportion)
    }
  }
  
  # Create dynamic column name based on timepoint (e.g., "Cluster_Proportion_Pre")
  custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
  
  # Build result data frame with patient IDs and proportions
  result_df <- data.frame(patient_id = subject_ids)
  result_df[[custom_column_name]] <- cluster_proportions
  
  return(result_df)
}

# =============================================================================================================
# FUNCTION: get_significance
# =============================================================================================================
# PURPOSE: Convert p-values to significance symbols for visualization
#
# @param p: P-value from statistical test (numeric)
# @return: Significance symbol (character)
#          "***" if p < 0.001 (highly significant)
#          "**"  if p < 0.01  (very significant)
#          "*"   if p < 0.05  (significant)
#          "ns"  if p >= 0.05 (not significant) or NA
# =============================================================================================================
get_significance <- function(p) {
  if (is.na(p)) {
    return("ns")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# =============================================================================================================
# Initialize Results Data Frame for Supplemental Table S3
# =============================================================================================================
# This data frame stores statistical comparison results for each cell type and timepoint pair
# Columns:
#   - celltype: Name of the cell population (e.g., Classical_Monocytes)
#   - timepoint_A: Earlier timepoint in comparison (e.g., Pre)
#   - timepoint_B: Later timepoint in comparison (e.g., C1)
#   - comparison: Type of comparison (control vs experiment = treatment arm comparison)
#   - p_value: Wilcoxon rank-sum test p-value comparing proportion changes between arms
#   - significance: Symbol representation of p-value (*, **, ***, ns)
#   - fold_change: Median log2 fold change across all patients (direction and magnitude)
comparison_results <- data.frame(
  celltype = character(),
  timepoint_A = character(),
  timepoint_B = character(),
  comparison = character(),
  p_value = numeric(),
  significance = character(),
  fold_change = numeric(),
  stringsAsFactors = FALSE
)

# =============================================================================================================
# Load Clinical Data and Configure Analysis Parameters
# =============================================================================================================

# Load patient survival and treatment arm data (filtered for UF site only)
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]  # UF site patients only

# Define metadata column names
patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# Define proportion change calculation method:
#   "difference" = Proportion_B - Proportion_A (absolute change)
#   "ratio" = Proportion_B / Proportion_A (fold change)
method = "difference"

# Define timepoints and create all pairwise combinations
# Timepoints: Pre (Pre-treatment), C1 (Cycle 1), C2 (Cycle 2)
# Pairs: (Pre, C1), (Pre, C2), (C1, C2)
timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

# Load Seurat object containing ALL cell types (myeloid + NK cells)
seurat_object <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
seurat_metadata <- seurat_object@meta.data

# =============================================================================================================
# Define Cell Type to Cluster Mapping
# =============================================================================================================
# Maps biological cell types to their corresponding cluster numbers from Seurat clustering
# These mappings are based on marker gene expression and functional annotation
#
# Cell Type Definitions:
#   - Classical Monocytes: CD14+ monocytes, major circulating monocyte population
#   - Non-Classical Monocytes: CD16+ monocytes, patrolling and inflammatory functions
#   - cDC (Conventional Dendritic Cells): Professional antigen-presenting cells
#   - pDC (Plasmacytoid Dendritic Cells): Type I interferon producers
#   - NK (Natural Killer): Cytotoxic innate lymphocytes
mapping <- list(
  "Classical_Monocytes" = c(1, 3, 16, 12),      # Multiple clusters represent Classical Monocytes
  "Non_Classical_Monocytes" = c(9, 33),         # CD16+ monocytes
  "cDC" = c(31),                                # Conventional dendritic cells
  "pDC" = c(36),                                # Plasmacytoid dendritic cells
  "NK" = c(0)                                   # Natural killer cells
)

# Convert mapping list to data frame for easier processing
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# =============================================================================================================
# Main Analysis Loop: Iterate Over Each Cell Type
# =============================================================================================================
# For each cell type:
#   1. Calculate cluster proportions at each timepoint for all patients
#   2. Compare proportion changes between treatment arms (MK-3475 Alone vs MK-3475 + MLA)
#   3. Perform statistical testing (Wilcoxon rank-sum test)
#   4. Calculate fold changes and generate visualizations
for (celltype in unique(celltype_to_cluster$celltype)){
  # Get cluster numbers for current cell type
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  seurat_metadata <- seurat_metadata
  
  # Initialize lists to store plots for current cell type
  plot_list_1 <- list()  # Line plots showing individual patient trajectories
  plot_list_2 <- list()  # Boxplots comparing treatment arms
  
  # ---------------------------------------------------------------------------------
  # Inner Loop: Compare Each Pair of Timepoints
  # ---------------------------------------------------------------------------------
  # For each timepoint pair (e.g., Pre vs C1, Pre vs C2, C1 vs C2):
  #   1. Calculate cluster proportions at both timepoints
  #   2. Merge with treatment arm information
  #   3. Calculate proportion changes
  #   4. Perform statistical test
  #   5. Generate visualizations
  # ---------------------------------------------------------------------------------
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]  # Earlier timepoint
    timepoint_b <- timepoint_pair[2]  # Later timepoint
    
    # -----------------------------------------
    # STEP 1: Calculate Cluster Proportions
    # -----------------------------------------
    # Calculate proportion of cells in this cell type for each patient at timepoint A
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata,
                                                                patient_col,
                                                                timepoint_col,
                                                                cluster_col,
                                                                cluster_list,
                                                                timepoint_a)
    
    # Calculate proportion of cells in this cell type for each patient at timepoint B
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata,
                                                                patient_col,
                                                                timepoint_col,
                                                                cluster_col,
                                                                cluster_list,
                                                                timepoint_b)
    
    # -----------------------------------------
    # STEP 2: Merge Data and Add Treatment Arm
    # -----------------------------------------
    # Merge proportions from both timepoints by patient ID
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df,
                                  timepoint_b_cluster_proportion_df,
                                  by = "patient_id")
    
    # Add treatment arm information (MK-3475 Alone vs MK-3475 + MLA)
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")],
                                  by = "patient_id")
    
    # Set factor levels for consistent ordering in plots (Control first, then Experimental)
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    
    # Define custom colors for treatment arms (consistent across all plots)
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    # -----------------------------------------
    # STEP 3: Calculate Proportion Change
    # -----------------------------------------
    # Calculate change in cluster proportion between timepoints for each patient
    # Method can be "difference" (absolute change) or "ratio" (fold change)
    #
    # BIOLOGICAL INTERPRETATION:
    #   - Positive difference: Cell population expanded between timepoints
    #   - Negative difference: Cell population contracted between timepoints
    #   - Ratio > 1: Cell population increased
    #   - Ratio < 1: Cell population decreased
    if (method == "difference") {
      merged_proportion_df$proportion_change <-
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] -
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change <-
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] /
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    # -----------------------------------------
    # STEP 4: Calculate Log2 Fold Change
    # -----------------------------------------
    # Calculate median log2 fold change across all patients as summary statistic
    # This represents the overall magnitude and direction of change
    #
    # INTERPRETATION:
    #   - log2(FC) = 0: No change
    #   - log2(FC) = 1: 2-fold increase (doubling)
    #   - log2(FC) = -1: 2-fold decrease (halving)
    #   - log2(FC) = 2: 4-fold increase
    #   - Positive values: Population expanded
    #   - Negative values: Population contracted
    proportion_a_vals <- merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    proportion_b_vals <- merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")]
    
    # Add small epsilon to avoid division by zero and log(0)
    eps <- 1e-9
    log2_ratios <- log2((proportion_b_vals + eps) / (proportion_a_vals + eps))
    
    # Calculate median log2 fold change across all patients
    median_fold_change <- median(log2_ratios, na.rm = TRUE)
    if (is.na(median_fold_change)) {
      median_fold_change <- 0  # Fallback if all values are NA
    }
    
    # -----------------------------------------
    # STEP 5: Statistical Testing
    # -----------------------------------------
    # Wilcoxon Rank-Sum Test (Mann-Whitney U test):
    #   - Non-parametric test comparing two independent groups
    #   - Tests whether proportion changes differ between treatment arms
    #   - Null hypothesis: No difference in proportion changes between MK-3475 Alone and MK-3475 + MLA
    #   - Alternative hypothesis: Proportion changes differ between treatment arms
    #
    # Requirements:
    #   - Both treatment arms present in data
    #   - At least 2 samples per treatment arm (minimum for statistical comparison)
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test comparing proportion_change between treatment arms
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value <- test_result$p.value
    } else {
      # Not enough samples for statistical test
      p_value <- NA_real_
    }
    
    # Convert p-value to significance symbol
    significance <- get_significance(p_value)
    
    # -----------------------------------------
    # STEP 6: Store Results
    # -----------------------------------------
    # Append statistical results to comparison_results data frame
    # This will be exported as Supplemental Table S3
    comparison_results <- rbind(
      comparison_results,
      data.frame(
        celltype = celltype,
        timepoint_A = timepoint_a,
        timepoint_B = timepoint_b,
        comparison = "control vs experiment",
        p_value = p_value,
        significance = significance,
        fold_change = median_fold_change,  # Median log2 fold change
        stringsAsFactors = FALSE
      )
    )
    
    # -----------------------------------------
    # STEP 7: Visualization Preparation
    # -----------------------------------------
    # Remove 'Cluster_Proportion_' prefix from column names for cleaner plotting
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert wide format to long format for ggplot2 visualization
    # This creates one row per patient per timepoint
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b),
                   names_to = "Time_Point",
                   values_to = "Cluster_Proportion")
    
    # Ensure timepoints are in correct chronological order
    data_long$Time_Point <- factor(data_long$Time_Point,
                                   levels = c(timepoint_a, timepoint_b))
    
    # -----------------------------------------
    # VISUALIZATION 1: Line Plot (Patient Trajectories)
    # -----------------------------------------
    # Shows individual patient trajectories across timepoints
    # Each line represents one patient's cluster proportion change
    # Color indicates treatment arm (blue = MK-3475 Alone, red = MK-3475 + MLA)
    #
    # INTERPRETATION:
    #   - Rising lines: Patient's cluster proportion increased over time
    #   - Falling lines: Patient's cluster proportion decreased over time
    #   - Parallel patterns by color: Similar responses within treatment arm
    #   - Diverging patterns by color: Differential responses between treatment arms
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion,
                               group = patient_id, color = Arm)) +
      geom_line(size = 1) +                              # Connect timepoints for each patient
      geom_point(size = 3) +                             # Show individual timepoint values
      scale_color_manual(values = custom_colors) +       # Apply custom treatment arm colors
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line()
      )
    plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p
    
    # -----------------------------------------
    # VISUALIZATION 2: Boxplot (Treatment Arm Comparison)
    # -----------------------------------------
    # Compares distribution of proportion changes between treatment arms
    # Boxplot shows median, quartiles, and outliers
    # Jittered points show individual patient values
    #
    # INTERPRETATION:
    #   - Boxplot center (median): Typical proportion change for treatment arm
    #   - Box height (IQR): Variability in proportion changes
    #   - Individual points: Each patient's proportion change
    #   - P-value annotation: Statistical significance of difference between arms
    #   - Significant p-value (<0.05): Treatment arms have different proportion changes
    p_box <- ggplot(merged_proportion_df,
                    aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +
      geom_jitter(position = position_jitterdodge(jitter.width = 1,
                                                  dodge.width = 1.2),
                  size = 4, alpha = 0.8) +               # Show individual patient values
      scale_fill_manual(values = custom_colors) +        # Apply custom treatment arm colors
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = NULL,
           y = "Proportion Change") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(),
        legend.position = "none"
      ) +
      annotate("text", x = 1.5,
               y = max(merged_proportion_df$proportion_change, na.rm = TRUE),
               label = paste("p-value =", ifelse(is.na(p_value),
                                                 "NA",
                                                 round(p_value, 4))),
               size = 4, vjust = -0.5)                    # Display p-value above plot
    
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
  }
  
  # -----------------------------------------
  # Combine and Save Visualizations
  # -----------------------------------------
  # Combine all plots into a grid: one row of line plots, one row of boxplots
  # Each column represents a different timepoint comparison
  plot_list <- c(plot_list_1, plot_list_2)
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list),
                          nrow = 2, align = 'v')
  
  # Define output directory based on calculation method
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Save combined plot as PDF
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_All_Cells_Comparison.pdf")
  ggsave(filename = file_name, plot = final_plot, device = "pdf",
         width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}

#########################################################################################################################
## EXPORT SUPPLEMENTAL TABLE S3
#########################################################################################################################
# Save comparison results to CSV file - this becomes Supplemental Table S3
# Contains statistical comparison results for all 5 major cell types across all timepoint pairs
#
# TABLE COLUMNS:
#   - celltype: Cell population name
#   - timepoint_A: Earlier timepoint
#   - timepoint_B: Later timepoint
#   - comparison: Type of comparison (control vs experiment)
#   - p_value: Wilcoxon test p-value
#   - significance: Significance symbol (*, **, ***, ns)
#   - fold_change: Median log2 fold change across patients
#
# BIOLOGICAL INTERPRETATION:
#   - Rows with significant p-values: Cell populations with differential dynamics between treatment arms
#   - Positive fold_change: Population expanded more (or contracted less) over treatment
#   - Negative fold_change: Population contracted more (or expanded less) over treatment
output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/comparison_results.csv"
write.csv(comparison_results, file = output_csv_path, row.names = FALSE)

cat("Comparison results have been saved to:", output_csv_path, "\n")
cat("This file serves as Supplemental Table S3\n")








#########################################################################################################################
#########################################################################################################################
##
## SECTION 2: T CELL SUBPOPULATIONS ANALYSIS (Generates Supplemental Table S4)
##
## PURPOSE: Perform detailed analysis of T cell subpopulation dynamics comparing treatment arms
##
## BIOLOGICAL CONTEXT:
## T cells are adaptive immune cells critical for anti-tumor immunity. This section analyzes 13
## functionally distinct T cell subpopulations to understand how combination therapy (MK-3475 + MLA)
## affects T cell composition compared to monotherapy (MK-3475 Alone).
##
## Key T cell populations analyzed:
## - Activated CD4: Helper T cells with activation markers
## - Effector CD8: Cytotoxic T cells with immediate effector function
## - Effector Memory CD8: Long-lived cytotoxic cells with memory potential
## - Exhausted T: Dysfunctional T cells with reduced effector capacity
## - Stem-Like CD8: Self-renewing memory precursors
## - Central Memory CD8: Long-term memory cells
##
## DATA SOURCE: MK_T_Cells_seurat_obj.RDS (T cells only, subset from total dataset)
## COMPARISON: Same methodology as Section 1, applied to T cell subpopulations
## OUTPUT: Supplemental Table S4 (comparison_results_T_Cells.csv)
##
## RATIONALE:
## T cell dynamics are particularly important in checkpoint blockade therapy:
## - Activated/Effector cells: Immediate anti-tumor response
## - Memory cells: Long-term immune surveillance
## - Exhausted cells: Marker of immune dysfunction to be reversed by therapy
## - Stem-like cells: Source of durable anti-tumor immunity
##
#########################################################################################################################
#########################################################################################################################

# Load Required Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# =============================================================================================================
# Helper Function: Calculate Cluster Proportion (T Cell Version)
# =============================================================================================================
# NOTE: This is a duplicate of the function from Section 1, using tidyverse syntax
#       Functionality is identical - calculates proportion of cells in specified clusters
#       See Section 1 for detailed documentation
# =============================================================================================================
get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for the clusters of interest and the specified timepoint
  cluster_timepoint_data <- seurat_metadata %>%
    filter((!!sym(cluster_col)) %in% clusters_of_interest & (!!sym(timepoint_col)) == timepoint_of_interest)
  
  timepoint_data <- seurat_metadata %>%
    filter((!!sym(timepoint_col)) == timepoint_of_interest)
  
  # Group the filtered data by patient
  cluster_timepoint_grouped <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
  timepoint_grouped <- split(timepoint_data, timepoint_data[[patient_col]])
  
  # Initialize vectors to store results
  subject_ids <- c()
  cluster_proportions <- c()
  
  # Loop through each patient to calculate proportions
  for (patient_id in names(cluster_timepoint_grouped)) {
    patient_cluster_data <- cluster_timepoint_grouped[[patient_id]]
    patient_timepoint_data <- timepoint_grouped[[patient_id]]
    
    patient_cluster_cells <- nrow(patient_cluster_data)
    patient_total_cells <- nrow(patient_timepoint_data)
    
    if (patient_cluster_cells > 0) {
      cluster_proportion <- patient_cluster_cells / patient_total_cells
      subject_ids <- c(subject_ids, patient_id)
      cluster_proportions <- c(cluster_proportions, cluster_proportion)
    }
  }
  
  # Create a data frame with the results
  custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
  result_df <- data.frame(patient_id = subject_ids, stringsAsFactors = FALSE)
  result_df[[custom_column_name]] <- cluster_proportions
  
  return(result_df)
}

# =============================================================================================================
# Helper Function: Determine Significance
get_significance <- function(p) {
  if (is.na(p)) {
    return("ns")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# =============================================================================================================
# Initialize Results Data Frame for Supplemental Table S4
# =============================================================================================================
# This data frame stores statistical comparison results for T cell subpopulations
# Structure identical to Section 1, but applied to T cell clusters
comparison_results <- data.frame(
  celltype     = character(),
  timepoint_A  = character(),
  timepoint_B  = character(),
  comparison   = character(),
  p_value      = numeric(),
  significance = character(),
  fold_change  = numeric(),
  stringsAsFactors = FALSE
)

# =============================================================================================================
# Load Clinical Data and Configure Analysis Parameters
# =============================================================================================================

# Load patient survival and treatment arm data (UF site only)
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df %>% filter(site == "UF")

# Define metadata column names
patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# Define proportion change calculation method
# NOTE: Section 2 uses "ratio" method instead of "difference" used in Section 1
#       "ratio" may be more appropriate for T cells due to their dynamic expansion/contraction
method <- "ratio"

# Define timepoints and create pairwise combinations (same as Section 1)
timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- split(pair_combinations, col(pair_combinations))

# =============================================================================================================
# Load Seurat Object for T Cells
# =============================================================================================================
# This object contains ONLY T cells (subset from the full dataset used in Section 1)
# Clustering was performed specifically on T cell data to resolve finer T cell subtypes
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

# =============================================================================================================
# Define T Cell Subpopulation to Cluster Mapping
# =============================================================================================================
# Maps 13 functionally distinct T cell subpopulations to their cluster numbers
# Clusters determined by marker gene expression and functional annotation
#
# T CELL SUBPOPULATION DEFINITIONS:
#
# CD4+ T Helper Cells (regulate immune responses):
#   - Activated_CD4 (cluster 0): CD4+ T cells with activation markers (CD69, HLA-DR)
#   - Active_CD4 (cluster 5): CD4+ T cells in active state
#   - Naive_CD4 (clusters 6, 9, 18): Unprimed CD4+ T cells (CCR7+, CD45RA+)
#   - Memory_CD4 (cluster 7): Memory CD4+ T cells with recall capacity
#
# CD8+ Cytotoxic T Cells (kill tumor cells):
#   - Effector_CD8 (cluster 1): Terminally differentiated cytotoxic cells (GZMB+, PRF1+)
#   - Effector_Memory_Precursor_CD8 (cluster 2): Cells transitioning to memory state
#   - Stem_Like_CD8 (cluster 8): Self-renewing memory precursors (TCF7+, LEF1+)
#   - Effector_Memory_CD8 (cluster 10): Cytotoxic cells with memory markers
#   - Central_Memory_CD8 (cluster 12): Long-lived memory cells (CCR7+, CD62L+)
#   - GZMK_Effector_Memory_CD8 (cluster 13): GZMK+ effector memory subset
#
# Specialized T Cell Subsets:
#   - Exhausted_T (cluster 3): Dysfunctional T cells (PD-1+, TIM-3+, LAG-3+)
#   - Gamma_Delta_T (cluster 4): Innate-like T cells (TCRγδ+)
#   - Proliferating_Effector (clusters 14, 16, 17): Actively dividing effector T cells (MKI67+)
#
# BIOLOGICAL SIGNIFICANCE:
#   - Effector populations: Immediate anti-tumor activity
#   - Memory populations: Long-term immune surveillance
#   - Exhausted cells: Targets for checkpoint blockade reversal
#   - Stem-like cells: Potential for sustained response
mapping <- list(
  "Activated_CD4"                = c(0),          # Activated CD4+ helper cells
  "Effector_CD8"                 = c(1),          # Terminal effector CD8+ cells
  "Effector_Memory_Precursor_CD8"= c(2),          # Transitional effector-memory cells
  "Exhausted_T"                  = c(3),          # Exhausted/dysfunctional T cells
  "Gamma_Delta_T"                = c(4),          # γδ T cells
  "Active_CD4"                   = c(5),          # Active CD4+ T cells
  "Naive_CD4"                    = c(6, 9, 18),   # Naive CD4+ T cells (multiple clusters)
  "Memory_CD4"                   = c(7),          # Memory CD4+ T cells
  "Stem_Like_CD8"                = c(8),          # Stem-like CD8+ memory precursors
  "Effector_Memory_CD8"          = c(10),         # Effector memory CD8+ cells
  "Central_Memory_CD8"           = c(12),         # Central memory CD8+ cells
  "GZMK_Effector_Memory_CD8"     = c(13),         # GZMK+ effector memory CD8+ cells
  "Proliferating_Effector"       = c(14, 16, 17)  # Proliferating effector cells (multiple clusters)
)

# Convert mapping list to data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype),
             cluster  = mapping[[celltype]],
             stringsAsFactors = FALSE)
}))

# =============================================================================================================
# Main Analysis Loop: Iterate Over Each T Cell Subpopulation
# =============================================================================================================
# Same workflow as Section 1, but applied to T cell subpopulations:
#   1. Calculate cluster proportions at each timepoint
#   2. Compare changes between treatment arms
#   3. Perform Wilcoxon testing
#   4. Calculate fold changes
#   5. Generate visualizations
for (celltype in unique(celltype_to_cluster$celltype)) {
  # Get the list of clusters for the current celltype
  cluster_list <- celltype_to_cluster %>%
    filter(celltype == !!celltype) %>%
    pull(cluster)
  
  # Assign the T cell metadata
  seurat_metadata <- seurat_metadata_t_cells
  
  # Initialize lists to store plots
  plot_list_1 <- list()  # Line Plots (Dotplots)
  plot_list_2 <- list()  # Boxplots
  
  # Iterate Over Each Pair of Timepoints
  for (i in 1:ncol(pair_combinations)) {
    timepoint_pair <- pair_combinations[, i]
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    # Calculate Cluster Proportions for Both Timepoints
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(
      seurat_metadata, patient_col, timepoint_col, cluster_col,
      clusters_of_interest = cluster_list, timepoint_of_interest = timepoint_a
    )
    
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(
      seurat_metadata, patient_col, timepoint_col, cluster_col,
      clusters_of_interest = cluster_list, timepoint_of_interest = timepoint_b
    )
    
    # Merge Proportion Data Frames by Patient ID
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, 
                                  timepoint_b_cluster_proportion_df, 
                                  by = "patient_id")
    
    # Merge with Survival Data to Get Arm Information
    merged_proportion_df <- merge(merged_proportion_df, 
                                  survival_df %>% select(patient_id, Arm), 
                                  by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, 
                                       levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    
    # Define Custom Colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    # Calculate Proportion Change Based on Method
    if (method == "difference") {
      merged_proportion_df$proportion_change <- 
        merged_proportion_df[[paste("Cluster_Proportion", timepoint_b, sep = "_")]] -
        merged_proportion_df[[paste("Cluster_Proportion", timepoint_a, sep = "_")]]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change <- 
        merged_proportion_df[[paste("Cluster_Proportion", timepoint_b, sep = "_")]] /
        merged_proportion_df[[paste("Cluster_Proportion", timepoint_a, sep = "_")]]
    }
    
    # ------------------------------------------------------------------------------
    # NEW STEP: Calculate a "fold_change" (median log2 ratio) across all patients
    #           to indicate both magnitude and direction of change.
    proportion_a_vals <- merged_proportion_df[[paste("Cluster_Proportion", timepoint_a, sep = "_")]]
    proportion_b_vals <- merged_proportion_df[[paste("Cluster_Proportion", timepoint_b, sep = "_")]]
    
    eps <- 1e-9
    log2_ratios <- log2((proportion_b_vals + eps) / (proportion_a_vals + eps))
    median_fold_change <- median(log2_ratios, na.rm = TRUE)
    
    if (is.na(median_fold_change)) {
      # Fallback if everything is NA
      median_fold_change <- 0
    }
    # ------------------------------------------------------------------------------
    
    # Perform Wilcoxon Test if Conditions are Met
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value <- test_result$p.value
    } else {
      p_value <- NA_real_
    }
    
    # Determine Significance
    significance <- get_significance(p_value)
    
    # Append Results to comparison_results Data Frame, Including fold_change
    comparison_results <- rbind(
      comparison_results,
      data.frame(
        celltype     = celltype,
        timepoint_A  = timepoint_a,
        timepoint_B  = timepoint_b,
        comparison   = "control vs experiment",
        p_value      = p_value,
        significance = significance,
        fold_change  = median_fold_change,  # <- NEW COLUMN
        stringsAsFactors = FALSE
      )
    )
    
    # Remove 'Cluster_Proportion_' Prefix from Column Names for Plotting
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert Data to Long Format for ggplot (Line/Dot Plots)
    data_long <- merged_proportion_df %>%
      pivot_longer(
        cols = c(timepoint_a, timepoint_b),
        names_to = "Time_Point",
        values_to = "Cluster_Proportion"
      )
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the Line Plot (Dotplot) without x-axis label
    p_line <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, 
                                    group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(
        title = "Cluster Proportion Change by Arm",
        x = NULL,  # Remove x-axis label
        y = "Cluster Proportion"
      ) +
      theme_minimal() +
      theme(
        axis.line         = element_line(color = "black"),
        axis.title.x      = element_blank(),
        axis.title.y      = element_text(size = 20),
        plot.title        = element_text(size = 20),
        axis.text.y       = element_text(size = 15),
        axis.text.x       = element_text(size = 12),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.ticks.y      = element_line()
      )
    
    # Add the Line Plot to the List
    plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p_line
    
    # Create the Boxplot
    p_box <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2), outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), 
                  size = 4, alpha = 0.8) +
      scale_fill_manual(values = custom_colors) +
      labs(
        title = paste("Change between", timepoint_a, "and", timepoint_b),
        x = NULL,
        y = "Proportion Change"
      ) +
      theme_minimal() +
      theme(
        axis.line         = element_line(color = "black"),
        axis.title.x      = element_blank(),
        axis.title.y      = element_text(size = 20),
        plot.title        = element_text(size = 20),
        axis.text.y       = element_text(size = 15),
        axis.text.x       = element_blank(),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.ticks.y      = element_line(),
        legend.position   = "none"
      ) +
      annotate(
        "text",
        x = 1.5,
        y = max(merged_proportion_df$proportion_change, na.rm = TRUE),
        label = ifelse(is.na(p_value), "p-value = NA", paste("p-value =", round(p_value, 4))),
        size = 4,
        vjust = -0.5
      )
    
    # Add the Boxplot to the List
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
  }
  
  # -----------------------------------------
  # Combine and Save T Cell Visualizations
  # -----------------------------------------
  # Combine all plots into grid format (same as Section 1)
  plot_list <- c(plot_list_1, plot_list_2)
  final_plot <- plot_grid(
    plotlist = plot_list,
    ncol = length(pairs_list),    # One column per timepoint comparison
    nrow = 2,                      # Two rows: line plots and boxplots
    align = 'v'
  )
  
  # Define output directory based on method (ratio vs difference)
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  
  # Create output directory if needed
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Save combined plot for this T cell subpopulation
  # File naming: [CellType]_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  ggsave(
    filename = file_name,
    plot     = final_plot,
    device   = "pdf",
    width    = 6 * length(pairs_list),
    height   = 6 * 2,
    limitsize= FALSE
  )
}

#########################################################################################################################
## EXPORT SUPPLEMENTAL TABLE S4
#########################################################################################################################
# Save T cell comparison results to CSV file - this becomes Supplemental Table S4
# Contains statistical comparison results for all 13 T cell subpopulations across all timepoint pairs
#
# TABLE STRUCTURE: Identical to Supplemental Table S3, but for T cell subpopulations
# COLUMNS:
#   - celltype: T cell subpopulation name (e.g., Effector_CD8, Exhausted_T)
#   - timepoint_A: Earlier timepoint
#   - timepoint_B: Later timepoint
#   - comparison: Type of comparison (control vs experiment)
#   - p_value: Wilcoxon test p-value comparing treatment arms
#   - significance: Significance symbol (*, **, ***, ns)
#   - fold_change: Median log2 fold change across patients
#
# CLINICAL INTERPRETATION:
# Significant changes in T cell subpopulations can indicate:
#   - Increased Effector/Memory cells: Enhanced anti-tumor immunity
#   - Decreased Exhausted cells: Successful checkpoint blockade reversal
#   - Increased Stem-like cells: Potential for durable response
#   - Differential patterns between arms: Combination therapy benefit
#
# NOTE: File saved to "ratio" folder (method = "ratio" in Section 2)
output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/comparison_results_T_Cells.csv"
write.csv(comparison_results, file = output_csv_path, row.names = FALSE)

cat("T cell comparison results have been saved to:", output_csv_path, "\n")
cat("This file serves as Supplemental Table S4\n")

#########################################################################################################################
## END OF SCRIPT
#########################################################################################################################
