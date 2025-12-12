# ====================================================================================================
# PATHWAY ACTIVITY COMPARISON ANALYSIS - MANUSCRIPT FIGURES 6B AND 6C
# ====================================================================================================
#
# PURPOSE:
# This script calculates pathway activity scores across T cell subtypes and compares treatment cohorts
# to generate Figures 6b and 6c for the manuscript. Pathway activity is defined as the mean expression
# of all genes within a specific pathway, representing the overall functional state of that pathway
# in different cell populations.
#
# BIOLOGICAL RATIONALE:
# Pathway activity scores provide a robust measure of cellular functional state by aggregating
# expression across multiple genes in a biological pathway. This approach:
# 1. Reduces noise compared to single-gene analysis
# 2. Captures coordinated biological processes (e.g., immune activation, T cell response)
# 3. Enables comparison of functional states between treatment cohorts (MK-3475 Alone vs MK-3475 + MLA)
# 4. Reveals temporal dynamics of immune pathway activation across treatment timepoints
#
# KEY METHODOLOGY:
# - Pathway Activity: Mean expression of all genes in pathway (from GMT file)
# - Cluster Proportion Weighting: Adjusts pathway signal by the frequency of cell subtype
# - Signal Change Metrics:
#   * Ratio method: (Activity at timepoint B) / (Activity at timepoint A)
#   * Difference method: (Activity at timepoint B) - (Activity at timepoint A)
# - Statistical Testing: Wilcoxon rank-sum test comparing cohorts at each timepoint pair
#
# MANUSCRIPT FIGURES GENERATED:
# - Figure 6b: Pathway activity changes across timepoints for different T cell subtypes
# - Figure 6c: Statistical comparison of pathway activity between treatment cohorts
#
# SCRIPT STRUCTURE:
# - SECTION 1 (Lines 1-426): Original implementation with basic functions and analysis
# - SECTION 2 (Lines 432-862): Refactored implementation with enhanced plotting and significance symbols
#
# ====================================================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)


# ====================================================================================================
# SECTION 1: ORIGINAL IMPLEMENTATION - PATHWAY ACTIVITY CALCULATION AND COHORT COMPARISON
# ====================================================================================================
# This section contains the initial implementation analyzing 10 immune pathways across 13 T cell
# subtypes, comparing MK-3475 Alone vs MK-3475 + MLA cohorts across multiple timepoints.
# ====================================================================================================

# ----------------------------------------------------------------------------------------------------
# FUNCTION: create_survival_data
# ----------------------------------------------------------------------------------------------------
# PURPOSE:
# Calculate pathway activity scores for each patient at multiple timepoints by computing the mean
# expression of all genes within a specified pathway.
#
# PARAMETERS:
# @param gmt_file: Path to GMT file containing pathway gene sets
# @param pathway_name: Name of the specific pathway to analyze (e.g., "GO_T_CELL_ACTIVATION")
# @param seurat_obj: Seurat object containing normalized single-cell RNA-seq data
# @param survival_data: Path to CSV file with patient clinical and survival data
#
# RETURNS:
# Data frame with patient IDs and mean pathway expression at each timepoint (Mean_Expr_Pre, etc.)
#
# METHODOLOGY:
# 1. Extract pathway gene list from GMT file
# 2. For each patient and timepoint:
#    - Subset cells belonging to that patient
#    - Calculate mean expression of pathway genes across those cells (pathway activity score)
# 3. Pathway activity = average normalized expression of all genes in the pathway
#
# BIOLOGICAL INTERPRETATION:
# Higher pathway activity indicates increased expression of genes in that pathway, suggesting
# enhanced functional activity (e.g., higher T cell activation pathway activity = more activated cells)
# ----------------------------------------------------------------------------------------------------
create_survival_data <- function(gmt_file, pathway_name, seurat_obj, survival_data) {
  # Read GMT file and extract genes for the specified pathway
  gmt <- readLines(gmt_file)
  pathway_genes <- NULL
  for (line in gmt) {
    split_line <- strsplit(line, "\t")[[1]]
    if (split_line[1] == pathway_name) {
      pathway_genes <- split_line[-c(1,2)]  # Assuming the first two columns are pathway name and description
      break
    }
  }
  
  # Check if pathway was found
  if (is.null(pathway_genes)) {
    stop("Pathway not found in GMT file.")
  }
  
  # Read survival data and sort patients by survival
  survival_df <- read.csv(survival_data)
  sorted_patients <- survival_df[order(-survival_df$OS.months.), 'patient_id']
  print(sorted_patients)
  
  # Create a grid of feature plots
  timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")
  # Initialize columns for mean expression at each timepoint in survival_df
  for(tp in timepoints) {
    survival_df[[paste0("Mean_Expr_", tp)]] <- NA_real_
  }
  
  for (patient in sorted_patients) {
    # print(patient)
    patient_specific_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$donor == patient,])
    if (length(patient_specific_cells) == 0) {
      next
    }
    patient_seurat_obj <- subset(seurat_obj, cells = patient_specific_cells)
    # Get the genes present in the Seurat object
    seurat_genes <- rownames(patient_seurat_obj@assays$RNA@counts)
    
    # Find intersection of pathway genes and Seurat object genes
    common_genes <- intersect(pathway_genes, seurat_genes)
    
    # Check if there are any common genes
    if (length(common_genes) == 0) {
      stop("None of the pathway genes are found in the Seurat object.")
    }
    
    # Calculate average expression of pathway genes per cell
    # Note: This assumes that the data is already normalized
    pathway_avg_expression <- tryCatch({
      apply(GetAssayData(patient_seurat_obj, assay = "RNA", slot = "data")[common_genes, ], 2, mean, na.rm = TRUE)
    }, error = function(e) {
      for (tp in timepoints) {
        # Update survival_df with mean expression for this patient and timepoint
        survival_df[survival_df$patient_id == patient, paste0("Mean_Expr_", tp)] <- NA
      }
      return(NULL)
    })
    
    # Skip the rest of the loop if an error occurred
    if (is.null(pathway_avg_expression)) {
      next
    }
    
    # # Check for cells at timepoint N
    # baseline_timepoint <- intersect(timepoints, unique(patient_seurat_obj@meta.data$TimePoint))[1]
    # timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == baseline_timepoint, ])
    # 
    # if (length(timepoint_N_barcodes) > 0) {
    #   # Calculate mean and standard deviation for cells at timepoint N
    #   mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
    #   sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])
    #   print(sd_N)
    #   if (sd_N == 0 | is.na(sd_N)){
    #     scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    #   } else {
    #     # Scale based on timepoint N
    #     scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
    #   }
    # } else {
    #   # Use scale function if no cells at timepoint N
    #   scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    # }
    # 
    # pathway_avg_expression <- scaled_pathway_avg_expression
    names(pathway_avg_expression) <- colnames(patient_seurat_obj)
    
    # Add this as a metadata column
    patient_seurat_obj[["pathway_avg_expression"]] <- pathway_avg_expression
    
    for (tp in timepoints) {
      patient_tp_cells <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == tp, ])
      if (length(patient_tp_cells) > 0) {
        # Calculate mean expression for this timepoint
        mean_expression <- mean(pathway_avg_expression[patient_tp_cells], na.rm = TRUE)
      } else {
        mean_expression <- NA
      }
      # Update survival_df with mean expression for this patient and timepoint
      survival_df[survival_df$patient_id == patient, paste0("Mean_Expr_", tp)] <- mean_expression
    }
  }
  # Return the updated survival dataframe
  return(survival_df)
}




# ----------------------------------------------------------------------------------------------------
# FUNCTION: get_cluster_proportion
# ----------------------------------------------------------------------------------------------------
# PURPOSE:
# Calculate the proportion of cells belonging to specific clusters for each patient at a given timepoint.
# This proportion is used to weight pathway activity scores, accounting for the relative abundance
# of specific T cell subtypes in each patient's sample.
#
# PARAMETERS:
# @param seurat_metadata: Metadata from Seurat object containing cell-level information
# @param patient_col: Column name containing patient/donor identifiers
# @param timepoint_col: Column name containing timepoint information
# @param cluster_col: Column name containing cluster assignments
# @param clusters_of_interest: Vector of cluster IDs to include in proportion calculation
# @param timepoint_of_interest: Specific timepoint to analyze (e.g., "Pre", "C1", "C2")
#
# RETURNS:
# Data frame with patient_id and Cluster_Proportion_X (where X is the timepoint)
#
# METHODOLOGY:
# Cluster proportion = (# cells in specified clusters) / (# total cells at timepoint)
# This weighting accounts for inter-patient variability in cell subset frequencies
#
# BIOLOGICAL RATIONALE:
# Different patients may have varying proportions of specific T cell subtypes. By weighting
# pathway activity by cluster proportion, we account for this heterogeneity and calculate
# a more accurate representation of overall pathway signal contribution.
# ----------------------------------------------------------------------------------------------------
get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for clusters of interest and specified timepoint
  cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest & seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  timepoint_data <- seurat_metadata[seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  # Group the filtered data by patient
  cluster_timpoint_grouped_data <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
  timpoint_grouped_data <- split(timepoint_data, timepoint_data[[patient_col]])
  
  # Initialize lists to store patient_id and Cluster_Proportion_X
  subject_ids <- c()
  cluster_proportions <- c()
  
  # Loop through each patient
  for (key in names(cluster_timpoint_grouped_data)) {
    patient_cluster_timepoint_data <- cluster_timpoint_grouped_data[[key]]
    patient_timepoint_data <- timpoint_grouped_data[[key]]
    
    # Calculate the total number of cells for the patient
    patient_cluster_timepoint_cells <- nrow(patient_cluster_timepoint_data)
    patient_timepoint_cells <- nrow(patient_timepoint_data)
    
    # Add the total number of cells to the results if there are any cells
    if (patient_cluster_timepoint_cells > 0) {
      # Calculate the cluster proportion
      cluster_proportion <- patient_cluster_timepoint_cells / patient_timepoint_cells
      
      # Append the patient_id and Cluster_Proportion_X to the lists
      subject_ids <- c(subject_ids, key)
      cluster_proportions <- c(cluster_proportions, cluster_proportion)
    }
  }
  
  custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
  # Create a data frame with patient_id and custom Cluster_Proportion_X column
  result_df <- data.frame(patient_id = subject_ids)
  result_df[[custom_column_name]] <- cluster_proportions
  
  return(result_df)
}



# ====================================================================================================
# MAIN ANALYSIS LOOP: PATHWAY ACTIVITY COMPARISON ACROSS T CELL SUBTYPES
# ====================================================================================================
# This section analyzes 10 immune pathways across 13 different T cell subtypes, comparing pathway
# activity between MK-3475 Alone and MK-3475 + MLA treatment cohorts at multiple timepoints.
# The analysis generates visualizations showing temporal pathway dynamics and statistical comparisons.
# ====================================================================================================

# Load T cell Seurat object (contains pre-clustered T cell populations)

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data


survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# ANALYSIS PARAMETERS
include_cluster_proportion <- TRUE  # Weight pathway activity by cluster proportion
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

# Define 10 immune-related pathways to analyze
pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

# Signal change calculation method
# "ratio": Calculate fold-change (Activity_B / Activity_A) - preferred for multiplicative effects
# "difference": Calculate absolute change (Activity_B - Activity_A) - preferred for additive effects
method = "ratio"

# ----------------------------------------------------------------------------------------------------
# T CELL SUBTYPE DEFINITIONS
# ----------------------------------------------------------------------------------------------------
# Map 13 biological T cell subtypes to their corresponding Seurat cluster IDs
# This mapping enables analysis of pathway activity in functionally distinct T cell populations
# ----------------------------------------------------------------------------------------------------
mapping <- list(
  "Activated_CD4" = c(0),
  "Effector_CD8" = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  "Active_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "GZMK_Effector_Memory_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))


timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()

# ----------------------------------------------------------------------------------------------------
# NESTED LOOP: ITERATE THROUGH PATHWAYS AND CELL TYPES
# ----------------------------------------------------------------------------------------------------
# For each pathway-celltype combination, calculate pathway activity at multiple timepoints,
# compare cohorts, perform statistical testing, and generate visualizations
# ----------------------------------------------------------------------------------------------------
for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    
    # Subset Seurat object to specific T cell subtype clusters
    seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
    
    # Normalize data using Relative Counts (RC) method scaled to 1 million
    # This ensures comparable expression values across cells
    seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    tcell_cluster_list <- cluster_list
    seurat_metadata <- seurat_metadata_t_cells
    
    # Calculate pathway activity (mean expression of pathway genes) for each patient and timepoint
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
    
    # Exclude IDH-positive patients (these have different biology and are analyzed separately)
    updated_survival_df <- updated_survival_df[updated_survival_df$IDH != "POS",]
    print(updated_survival_df)
    
    # Initialize a list to store the plots
    plot_list_1 <- list()
    plot_list_2 <- list()
    
    for (timepoint_pair in pairs_list){
      temp_survival_df <- updated_survival_df
      timepoint_a <- timepoint_pair[1]
      timepoint_b <- timepoint_pair[2]
      
      if (include_cluster_proportion) {
        # Calculate cluster proportions for both timepoints
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
        
        # Merge cluster proportion data with survival data
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        
        # CLUSTER PROPORTION WEIGHTING:
        # Weight pathway activity by the proportion of cells in this subtype
        # Pathway Signal = (Mean pathway expression) × (Proportion of cells in this cluster)
        # This accounts for cell subset frequency variations across patients
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
        # SIGNAL CHANGE CALCULATION:
        # Quantify how pathway activity changes between timepoint A and timepoint B
        if (method == "difference") {
          # Difference method: Absolute change in pathway signal
          # Positive = increased activity; Negative = decreased activity
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          # Ratio method: Fold-change in pathway signal
          # >1 = increased activity; <1 = decreased activity; 1 = no change
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      } else {
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      }
      print(temp_survival_df[, "signal_change"])
      
      merged_proportion_df <- temp_survival_df
      # COHORT COMPARISON SETUP:
      # Compare two treatment arms: MK-3475 alone (monotherapy) vs MK-3475 + MLA (combination therapy)
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      
      # Define custom colors for visualization
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      # STATISTICAL TESTING:
      # Perform Wilcoxon rank-sum test only if both cohorts have sufficient samples (n >= 2)
      if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
          sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
          sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
        
        # Calculate the median signal change for each group
        medians <- merged_proportion_df %>%
          group_by(Arm) %>%
          summarize(median_signal_change = median(signal_change, na.rm = TRUE))
        
        # Compute the log fold change (log2)
        control_median <- medians$median_signal_change[medians$Arm == "MK-3475 Alone"]
        treatment_median <- medians$median_signal_change[medians$Arm == "MK-3475 + MLA"]
        if (is.na(treatment_median) | is.na(control_median)){
          fold_change <- NA
          log_fold_change <- NA
        } else {
          if (method == "difference") {
            fold_change <- treatment_median - control_median
            log_fold_change <- NA
          } else {
            # Calculate absolute fold change
            if (control_median * treatment_median < 0) {
              # If medians have different signs, treat fold change as undefined
              fold_change <- NA
              log_fold_change <- NA
            } else {
              # Calculate fold change
              fold_change <- abs(treatment_median / control_median)
              
              # Calculate log fold change
              log_fold_change <- log2(abs(treatment_median / control_median))
              
              # Apply the sign to log fold change
              if (treatment_median < 0 & control_median < 0) {
                log_fold_change <- -log_fold_change
              } else if (treatment_median > 0 & control_median < 0) {
                log_fold_change <- -log_fold_change
              }
            }
          }
        }
        
        
        # WILCOXON RANK-SUM TEST:
        # Non-parametric test comparing signal_change distributions between treatment cohorts
        # Null hypothesis: No difference in pathway activity change between cohorts
        # Tests whether one cohort shows significantly different pathway dynamics
        p_value <- tryCatch({
          test_result <- wilcox.test(signal_change ~ Arm, data = merged_proportion_df)
          test_result$p.value
        }, error = function(e) {
          NA_real_
        })
      } else {
        p_value = NA_real_  # Use NA_real_ for numerical context missing values
        fold_change <- NA_real_
        log_fold_change <- NA_real_
      }
      
      pathway_list <- c(pathway_list, pathway)
      celltype_list <- c(celltype_list, celltype)
      comparison_list <- c(comparison_list, paste(timepoint_a, timepoint_b, sep = "_vs_"))
      fc_list <- c(fc_list, fold_change)
      logfc_list <- c(logfc_list, log_fold_change)
      pvalue_list <- c(pvalue_list, p_value)
      
      # Remove the 'pre_' prefix from column names
      colnames(merged_proportion_df) <- sub("^Pathway_Signal_", "", colnames(merged_proportion_df))
      
      # Convert data to long format for ggplot
      data_long <- merged_proportion_df %>%
        pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Pathway_Signal")
      
      data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
      
      # Create the plot
      p <- ggplot(data_long, aes(x = Time_Point, y = Pathway_Signal, group = patient_id, color = Arm)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        scale_color_manual(values = custom_colors) +
        labs(title = "Pathway Signal Change by Arm",
             x = "Time Point",
             y = "Pathway Signal") +
        theme_minimal() +
        theme(axis.line = element_line(color = "black"),
              axis.title.x = element_text(size = 20),  # Adjust x axis label size
              axis.title.y = element_text(size = 20),  # Adjust y axis label size
              plot.title = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),  # Add axis lines
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              axis.ticks.y = element_line())
      
      # Add the plot to the list
      plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p
      
      # Create boxplot
      p <- ggplot(merged_proportion_df, aes(x = Arm, y = signal_change, fill = Arm)) +
        geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
        geom_jitter(position = position_dodge(width = 1.2)) +
        scale_fill_manual(values = custom_colors) +
        labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
             x = "Arm",
             y = "Signal Change") +
        theme_minimal() +
        theme(axis.line = element_line(color = "black"),
              axis.title.x = element_text(size = 20),  # Adjust x axis label size
              axis.title.y = element_text(size = 20),  # Adjust y axis label size
              plot.title = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),  # Add axis lines
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              axis.ticks.y = element_line()) +
        annotate("text", x = 1.5, y = max(merged_proportion_df$signal_change, na.rm = TRUE), 
                 label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
      
      # Add the plot to the list
      plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
    }
    
    # *** MANUSCRIPT FIGURE GENERATION (FIGURE 6B AND 6C) ***
    # Combine line plots (temporal dynamics) and boxplots (cohort comparisons) into grid
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    
    if(include_cluster_proportion){
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/ratio/"
      }
    } else {
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/ratio/"
      }
    }
    file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison.pdf")
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, FC = fc_list, logFC = logfc_list, pvalue = pvalue_list)
if(include_cluster_proportion){
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/difference/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/ratio/"
  }
} else {
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/difference/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/ratio/"
  }
}
file_name <- paste0(dir, "signal_comparison_in_T_cells.csv")
write.csv(df, file_name)






# ====================================================================================================
# SECTION 2: REFACTORED IMPLEMENTATION WITH ENHANCED VISUALIZATION
# ====================================================================================================
# This section contains an improved version of the analysis with:
# - Better code organization using dplyr/tidyverse syntax
# - Addition of significance symbols (*, **, ***) based on p-values
# - Improved plot aesthetics (removed x-axis labels, legends for clarity)
# - More robust error handling
#
# The core methodology remains identical to Section 1, calculating pathway activity scores,
# weighting by cluster proportions, and comparing treatment cohorts using Wilcoxon tests.
# ====================================================================================================

# Load Required Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# ----------------------------------------------------------------------------------------------------
# FUNCTION: get_significance
# ----------------------------------------------------------------------------------------------------
# PURPOSE:
# Convert p-values to standard significance symbols for publication-ready figures
#
# PARAMETERS:
# @param p: Numeric p-value from statistical test
#
# RETURNS:
# Character string with significance symbol:
#   "***" for p < 0.001 (highly significant)
#   "**"  for p < 0.01  (very significant)
#   "*"   for p < 0.05  (significant)
#   "ns"  for p >= 0.05 (not significant)
#
# USAGE:
# This function standardizes significance reporting across all comparisons and ensures
# consistent interpretation of statistical results in figures
# ----------------------------------------------------------------------------------------------------
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

# ==================================================================================================================================================================
# Helper Function to Calculate Cluster Proportions
# (Assuming this function is defined elsewhere in your script. If not, include its definition here.)

# get_cluster_proportion <- function(...) {
#   # Your existing implementation
# }

# ==================================================================================================================================================================
# Load Seurat Object for T Cells
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

# ==================================================================================================================================================================
# Load Survival Data
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df %>% filter(site == "UF")

# Define Columns
patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# Additional Parameters
include_cluster_proportion <- FALSE
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

# Define Pathways
pathways <- c(
  "GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE",
  "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
  "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY",
  "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION",
  "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"
)

# Define Method
method <- "ratio"  # Options: "difference" or "ratio"

# Define the Mapping from Celltypes to Clusters
mapping <- list(
  "Activated_CD4" = c(0),
  "Effector_CD8" = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  "Active_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "GZMK_Effector_Memory_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17)
)

# Convert the Mapping to a Data Frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]], stringsAsFactors = FALSE)
}))

# Define Timepoints and Pair Combinations
timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

# Initialize Vectors to Store Comparison Results
pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()
significance_list <- c()

# ==================================================================================================================================================================
# Iterate Over Each Pathway and Celltype
for (pathway in pathways) {
  for (celltype in unique(celltype_to_cluster$celltype)) {
    print(celltype)
    
    # Get the list of clusters for the current celltype
    cluster_list <- celltype_to_cluster %>%
      filter(celltype == !!celltype) %>%
      pull(cluster)
    
    # Subset the Seurat Object for T Cells of Interest
    seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_t_cell_subset <- NormalizeData(
      seurat_object_t_cell_subset,
      normalization.method = "RC",
      scale.factor = 1e6,
      verbose = TRUE,
      assay = "RNA"
    )
    
    tcell_cluster_list <- cluster_list
    seurat_metadata <- seurat_metadata_t_cells
    
    # Update Survival Data with Pathway Information
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
    updated_survival_df <- updated_survival_df %>% filter(IDH != "POS")
    print(updated_survival_df)
    
    # Initialize Lists to Store Plots
    plot_list_1 <- list()  # Line Plots (Dotplots)
    plot_list_2 <- list()  # Boxplots
    
    # Iterate Over Each Pair of Timepoints
    for (timepoint_pair in pairs_list) {
      temp_survival_df <- updated_survival_df
      timepoint_a <- timepoint_pair[1]
      timepoint_b <- timepoint_pair[2]
      
      if (include_cluster_proportion) {
        # Calculate Cluster Proportions for Both Timepoints
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(
          seurat_metadata, patient_col, timepoint_col, cluster_col,
          clusters_of_interest = tcell_cluster_list, timepoint_of_interest = timepoint_a
        )
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(
          seurat_metadata, patient_col, timepoint_col, cluster_col,
          clusters_of_interest = tcell_cluster_list, timepoint_of_interest = timepoint_b
        )
        
        # Merge Cluster Proportion Data
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        
        # Calculate Pathway Signal
        temp_survival_df <- temp_survival_df %>%
          mutate(
            !!paste0("Pathway_Signal_", timepoint_a) := !!sym(paste0("Mean_Expr_", timepoint_a)) * !!sym(paste0("Cluster_Proportion_", timepoint_a)),
            !!paste0("Pathway_Signal_", timepoint_b) := !!sym(paste0("Mean_Expr_", timepoint_b)) * !!sym(paste0("Cluster_Proportion_", timepoint_b))
          )
        
        # Calculate Signal Change
        if (method == "difference") {
          temp_survival_df <- temp_survival_df %>%
            mutate(signal_change = !!sym(paste0("Pathway_Signal_", timepoint_b)) - !!sym(paste0("Pathway_Signal_", timepoint_a)))
        } else if (method == "ratio") {
          temp_survival_df <- temp_survival_df %>%
            mutate(signal_change = !!sym(paste0("Pathway_Signal_", timepoint_b)) / !!sym(paste0("Pathway_Signal_", timepoint_a)))
        }
      } else {
        # Calculate Pathway Signal Without Cluster Proportion
        temp_survival_df <- temp_survival_df %>%
          mutate(
            !!paste0("Pathway_Signal_", timepoint_a) := !!sym(paste0("Mean_Expr_", timepoint_a)),
            !!paste0("Pathway_Signal_", timepoint_b) := !!sym(paste0("Mean_Expr_", timepoint_b))
          )
        
        # Calculate Signal Change
        if (method == "difference") {
          temp_survival_df <- temp_survival_df %>%
            mutate(signal_change = !!sym(paste0("Pathway_Signal_", timepoint_b)) - !!sym(paste0("Pathway_Signal_", timepoint_a)))
        } else if (method == "ratio") {
          temp_survival_df <- temp_survival_df %>%
            mutate(signal_change = !!sym(paste0("Pathway_Signal_", timepoint_b)) / !!sym(paste0("Pathway_Signal_", timepoint_a)))
        }
      }
      
      print(temp_survival_df[, "signal_change"])
      
      # Assign to merged_proportion_df for Consistency
      merged_proportion_df <- temp_survival_df
      
      # Convert Arm to Factor with Specified Levels
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      
      # Define Custom Colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      # Perform Wilcoxon Test if Conditions are Met
      if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
          sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
          sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
        
        # Calculate Medians (Optional, Not Used Here)
        medians <- merged_proportion_df %>%
          group_by(Arm) %>%
          summarize(median_signal_change = median(signal_change, na.rm = TRUE))
        
        # Compute Fold Change (Optional, Not Used Here)
        control_median <- medians$median_signal_change[medians$Arm == "MK-3475 Alone"]
        treatment_median <- medians$median_signal_change[medians$Arm == "MK-3475 + MLA"]
        
        if (is.na(treatment_median) | is.na(control_median)) {
          fold_change <- NA
          log_fold_change <- NA
        } else {
          if (method == "difference") {
            fold_change <- treatment_median - control_median
            log_fold_change <- NA
          } else {
            # Calculate absolute fold change
            if (control_median * treatment_median < 0) {
              # If medians have different signs, treat fold change as undefined
              fold_change <- NA
              log_fold_change <- NA
            } else {
              # Calculate fold change
              fold_change <- abs(treatment_median / control_median)
              
              # Calculate log fold change
              log_fold_change <- log2(abs(treatment_median / control_median))
              
              # Apply the sign to log fold change
              if (treatment_median < 0 & control_median < 0) {
                log_fold_change <- -log_fold_change
              } else if (treatment_median > 0 & control_median < 0) {
                log_fold_change <- -log_fold_change
              }
            }
          }
        }
        
        # Perform Wilcoxon Test
        p_value <- tryCatch({
          test_result <- wilcox.test(signal_change ~ Arm, data = merged_proportion_df)
          test_result$p.value
        }, error = function(e) {
          NA_real_
        })
      } else {
        p_value <- NA_real_
        fold_change <- NA_real_
        log_fold_change <- NA_real_
      }
      
      # Determine Significance
      significance <- get_significance(p_value)
      
      # Append Results to Comparison Vectors
      pathway_list <- c(pathway_list, pathway)
      celltype_list <- c(celltype_list, celltype)
      comparison_list <- c(comparison_list, paste(timepoint_a, timepoint_b, sep = "_vs_"))
      fc_list <- c(fc_list, fold_change)
      logfc_list <- c(logfc_list, log_fold_change)
      pvalue_list <- c(pvalue_list, p_value)
      significance_list <- c(significance_list, significance)
      
      # Remove 'Pathway_Signal_' Prefix from Column Names for Plotting
      colnames(merged_proportion_df) <- sub("^Pathway_Signal_", "", colnames(merged_proportion_df))
      
      # Convert Data to Long Format for ggplot (Line Plots)
      data_long <- merged_proportion_df %>%
        pivot_longer(
          cols = c(timepoint_a, timepoint_b),
          names_to = "Time_Point",
          values_to = "Pathway_Signal"
        )
      
      data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
      
      # Create the Line Plot (Dotplot) without x-axis label
      p_line <- ggplot(data_long, aes(x = Time_Point, y = Pathway_Signal, group = patient_id, color = Arm)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        scale_color_manual(values = custom_colors) +
        labs(
          title = "Pathway Signal Change by Arm",
          x = NULL,  # Remove x-axis label
          y = "Pathway Signal"
        ) +
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black"),
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_line()
        )
      
      # Add the Line Plot to the List
      plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p_line
      
      # Create the Boxplot without x-axis label, x-axis text, and legend
      p_box <- ggplot(merged_proportion_df, aes(x = Arm, y = signal_change, fill = Arm)) +
        geom_boxplot(width = 0.5, position = position_dodge(width = 1.2), outlier.shape = NA) +
        # geom_jitter(position = position_dodge(width = 1.2)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), size = 4, alpha = 0.8) +
        scale_fill_manual(values = custom_colors) +
        labs(
          title = paste("Change between", timepoint_a, "and", timepoint_b, 
                        ifelse(is.na(p_value), "(p-value = NA)", paste("(p-value =", round(p_value, 4), ")"))),
          x = NULL,  # Remove x-axis label
          y = "Signal Change"
        ) +
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black"),
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_blank(),  # Remove x-axis text (tick labels)
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_line(),
          legend.position = "none"  # Remove legend
        )
        # annotate(
        #   "text",
        #   x = 1.5,
        #   y = max(merged_proportion_df$signal_change, na.rm = TRUE),
        #   label = ifelse(is.na(p_value), "p-value = NA", paste("p-value =", round(p_value, 4))),
        #   size = 4,
        #   vjust = -0.5
        # )
      
      # Add the Boxplot to the List
      plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
    }
    
    # -----------------------------------------------------------------------------------------------
    # *** MANUSCRIPT FIGURE GENERATION (FIGURE 6B AND 6C) ***
    # Combine all plots (line plots + boxplots) into publication-ready grid layout
    # -----------------------------------------------------------------------------------------------
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(
      plotlist = plot_list,
      ncol = length(pairs_list),
      nrow = 2,
      align = 'v'
    )
    
    # Define Plot Directory Based on Method and Cluster Proportion Inclusion
    if (include_cluster_proportion) {
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/ratio/"
      }
    } else {
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/ratio/"
      }
    }
    
    # Create the Plot Directory if It Doesn't Exist
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    
    # Define the File Name and Save the Plot
    file_name <- paste0(plot_dir, celltype, "_", pathway, "_Signal_Change_Comparison.pdf")
    ggsave(
      filename = file_name,
      plot = final_plot,
      device = "pdf",
      width = 7 * length(pairs_list),
      height = 6 * 2,
      limitsize = FALSE
    )
  }
}

# ====================================================================================================
# EXPORT STATISTICAL RESULTS
# ====================================================================================================
# Save comprehensive statistical comparison results including pathway, celltype, timepoint comparison,
# fold-change, log fold-change, p-values, and significance symbols for downstream analysis
# ====================================================================================================

# Combine All Comparison Vectors into a Data Frame
df <- data.frame(
  pathway = pathway_list,
  celltype = celltype_list,
  comparison = comparison_list,
  FC = fc_list,
  logFC = logfc_list,
  pvalue = pvalue_list,
  significance = significance_list,  # Add the significance column
  stringsAsFactors = FALSE
)

# Define Output Directory Based on Method and Cluster Proportion Inclusion
if (include_cluster_proportion) {
  if (method == "difference") {
    output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/difference/"
  } else {
    output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/with_cluster_proportion/ratio/"
  }
} else {
  if (method == "difference") {
    output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/difference/"
  } else {
    output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Pathway_Activity_Comparison_Results/without_cluster_proportion/ratio/"
  }
}

# Create the Output Directory if It Doesn't Exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the File Name for the CSV
csv_file_name <- paste0(output_dir, "signal_comparison_in_T_cells.csv")

# Save the Data Frame to CSV
write.csv(df, csv_file_name, row.names = FALSE)

# Optional: Print a Completion Message
cat("Comparison results have been saved to:", csv_file_name, "\n")
