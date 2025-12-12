################################################################################
# CLONAL EXPANSION ANALYSIS - STEP 4: VISUALIZATION (FIGURE 6d)
################################################################################
#
# PURPOSE:
#   This is the FINAL step (4/4) in the clonal expansion analysis pipeline.
#   It generates comprehensive visualizations of clonal expansion across
#   T cell subpopulations, timepoints, and patient cohorts.
#
# ANALYSIS PIPELINE:
#   Step 1: clonal_expansion_analysis_1.R - Prepare clonotype tables
#   Step 2: diversity_calculation_2.R - Calculate diversity indices
#   Step 3: clonal_expansion_calculation_3.R - Calculate expansion metrics
#   Step 4 (THIS SCRIPT): Generate visualizations → Figure 6d
#
# MANUSCRIPT FIGURE:
#   Figure 6d: Clonal expansion visualization showing diversity changes across
#              timepoints in different patient cohorts and T cell subpopulations
#
# VISUALIZATION APPROACHES:
#   1. Expansion Ratio Plots (Lines 15-183):
#      - Violin/boxplots of diversity ratios (C1/C2 vs 1.0)
#      - Separate plots per treatment arm
#      - T-test comparing ratio to 1.0 (null hypothesis: no change)
#
#   2. Paired Diversity Boxplots (Lines 187-577):
#      - Side-by-side boxplots showing actual diversity values
#      - Paired lines connecting the same patient across timepoints
#      - Faceted by patient cohort (control, short_term, long_term)
#      - Paired t-test for statistical significance
#      - Multiple timepoint comparisons (Pre vs C1, Pre vs C2, C1 vs C2)
#
# PATIENT COHORTS:
#   - Control: Patients 1, 4, 8, 9, 11 (no treatment/standard care)
#   - Short-term survivors: Patients 7, 10, 12, 14, 18
#   - Long-term survivors: Patients 2, 3, 5, 13, 19, 20, 21
#
# TIMEPOINTS:
#   - Pre (S): Baseline/pre-treatment
#   - C1: Cycle 1 (early treatment)
#   - C2: Cycle 2 (later treatment)
#
# STATISTICAL TESTS:
#   - T-test comparing expansion ratio to 1.0 (no change)
#   - Paired t-test for same patients across timepoints
#   - Significance threshold: p ≤ 0.05
#
# INPUT FILES:
#   - C1_vs_C2.txt: Expansion ratios (from Step 3)
#   - shannon/simpson_clonal_diversity_*.txt: Raw diversity values (from Step 2)
#   - patient_survival_data.csv: Clinical metadata
#
# OUTPUT FILES:
#   - {diversity_index}_Clonal_Diversity_{arm}_clonal_expansion.pdf
#   - {diversity_index}_Clonal_Diversity_{celltype}_allArms_clonal_expansion.pdf
#   - Paired boxplot PDFs organized by timepoint comparison
#
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(broom)
library(ggfortify)
library("survminer")
library("Rcpp")
library(cowplot)
library(tidyr)
library(rlang)

################################################################################
# VISUALIZATION APPROACH 1: EXPANSION RATIO PLOTS
################################################################################

################################################################################
# FUNCTION: make_comparison_value_plot
################################################################################
# Create violin/boxplot of clonal expansion ratios compared to 1.0
#
# PARAMETERS:
#   df - Dataframe containing expansion values
#   variable_name - Column name with expansion ratios
#   diversity_index - "shannon" or "simpson"
#   alternative - T-test alternative hypothesis ("l" = less than 1)
#   plot_dir - Output directory for plots
#   arm_type - Treatment arm name
#
# VISUALIZATION:
#   - Violin plot showing distribution
#   - Boxplot for quartiles
#   - Individual points for each patient
#   - Red dashed horizontal line at y=1 (no change)
#   - P-value annotation (red if significant)
#
# INTERPRETATION:
#   Values < 1: Expansion occurred (diversity decreased)
#   Values > 1: Diversification occurred (diversity increased)
#
################################################################################
make_comparison_value_plot = function(df,
                                      variable_name,
                                      diversity_index = "simpson",
                                      alternative = "l",
                                      plot_dir,
                                      arm_type) {
  # alternative: T-test hypothesis
  #   "l" = less than (testing if ratio < 1, i.e., expansion)
  #   "g" = greater than
  #   "two.sided" = different from
  
  print(df)
  
  # Filter out infinite values before statistical test
  finite_values <- df[, variable_name][is.finite(df[, variable_name])]
  
  # T-test comparing ratio to 1.0 (null hypothesis: no expansion)
  p_val = t.test(x = finite_values,
                 mu = 1,
                 alternative = alternative)[["p.value"]]
  
  # Format p-value label
  p_label = paste("p =", round(p_val, 3))
  
  # Color significant results in red
  if (p_val <= 0.05) {
    col = "red"
    p_label = paste(p_label, "*")
  } else {
    col = "black"
  }
  
  # Create plot labels
  plot_label = paste("Relative", variable_name)
  x_label = paste("relative", variable_name)
  
  # Filter dataframe for plotting
  df_plot <- df[is.finite(df[, variable_name]), ]
  
  # Create visualization
  # Shannon and Simpson have different y-axis positioning for p-value label
  if (diversity_index == "shannon") {
    p = ggplot(df_plot, aes(x = x_label , y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.2) +
      geom_boxplot(width = 0.1) +
      geom_point() +
      annotate(
        "text",
        label = p_label,
        x = 1,
        y = max(df_plot[, variable_name], na.rm = TRUE) * 1.5,
        size = 10,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 2
      )
  } else {
    p = ggplot(df_plot, aes(x = x_label , y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.2) +
      geom_boxplot(width = 0.1) +
      geom_point() +
      annotate(
        "text",
        label = p_label,
        x = 1,
        y = max(df_plot[, variable_name], na.rm = TRUE) * 1.1,
        size = 10,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 2
      )
  }
  
  # Apply clean theme
  p = p + theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_blank()
    ) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 10
    ))
  
  print(p)
  
  # Save plot
  file_name <- paste0(diversity_index, "_Clonal_Diversity_", arm_type, "_", variable_name, ".pdf")
  ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 4, height = 6, dpi = 300, limitsize = FALSE)
}

################################################################################
# SECTION 1: LOAD PATIENT METADATA AND DEFINE COHORTS
################################################################################

# Load clinical metadata
patient_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
patient_metadata_df <- read.csv(patient_metadata_file)
patient_metadata_df <- patient_metadata_df[patient_metadata_df$site == "UF",]

# Define patient cohorts based on treatment response and survival
control_group <- c(1, 4, 8, 9, 11)                    # No treatment/standard care
short_term_survivor_group <- c(7, 10, 12, 14, 18)     # Short survival after treatment
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21) # Long survival after treatment

# Assign group labels to each patient
patient_metadata_df$group <- ifelse(
  patient_metadata_df$patient_id %in% control_group, "control",
  ifelse(
    patient_metadata_df$patient_id %in% short_term_survivor_group, "short_term",
    ifelse(
      patient_metadata_df$patient_id %in% long_term_survivor_group, "long_term",
      NA_character_
    )
  )
)

################################################################################
# SECTION 2: DEFINE T CELL SUBPOPULATION MAPPINGS
################################################################################

mapping <- list(
  "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
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
  "Proliferating_Effector" = c(14, 16, 17),
  "All_CD8" = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)

celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

celltypes <- unique(celltype_to_cluster$celltype)

################################################################################
# SECTION 3: GENERATE EXPANSION RATIO PLOTS (Pre vs C1 Comparison)
################################################################################
# Plots showing diversity ratio distribution for each treatment arm

diversity_indices = c("shannon", "simpson")
dir = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/Pre_vs_C1"

for (celltype in celltypes) {
  print(celltype)
  
  for (diversity_index in diversity_indices){
    print(diversity_index)
    
    # Set plot directory for this cell type and diversity index
    plot_dir = paste0(dir, "/clonal_expansion_", diversity_index, "_", celltype, "_T_cells_division/")
    
    # Load expansion ratio data
    df <- read.delim(paste0(plot_dir, "S_vs_C1.txt"))
    
    # Merge with clinical metadata to get treatment arm information
    df_merged <- merge(df, 
                      patient_metadata_df[, c("patient_id", "Arm")], 
                      by.x = "row.names", 
                      by.y = "patient_id", 
                      all.x = TRUE, 
                      all.y = FALSE)
    
    # Generate separate plots for each treatment arm
    for (arm_type in unique(df_merged$Arm)){
      df <- df_merged[df_merged$Arm == arm_type, ]
      df <- df[!is.na(df$clonal_expansion), ]
      
      # Only plot if there are at least 2 data points
      if (nrow(df) > 1){
        make_comparison_value_plot(df, "clonal_expansion", 
                                  diversity_index = diversity_index, 
                                  alternative = "l",  # Test if ratio < 1 (expansion)
                                  plot_dir, 
                                  arm_type)
      } else {
        print("=======================================================================================================================")
        print(celltype)
        print(df)
      }
    }
  }
}

################################################################################
# VISUALIZATION APPROACH 2: PAIRED DIVERSITY BOXPLOTS
################################################################################
# Updated April 22: More sophisticated visualization showing:
#   - Actual diversity values (not ratios)
#   - Paired observations with connecting lines
#   - Multiple timepoint comparisons
#   - Cohort-stratified analysis
################################################################################

# Re-load libraries (redundant but ensures availability)
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(broom)
library(ggfortify)
library(survminer)
library(Rcpp)
library(cowplot)
library(tidyr)
library(rlang)

################################################################################
# SECTION 4: SETUP FOR PAIRED BOXPLOT VISUALIZATION
################################################################################

# Load patient metadata
patient_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
patient_metadata_df  <- fread(patient_metadata_file) %>%
  filter(site == "UF") %>%
  mutate(patient_id = as.integer(patient_id))

# Define patient cohorts
control_group             <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Assign cohort labels
patient_metadata_df <- patient_metadata_df %>%
  mutate(group = case_when(
    patient_id %in% control_group             ~ "control",
    patient_id %in% short_term_survivor_group ~ "short_term",
    patient_id %in% long_term_survivor_group  ~ "long_term",
    TRUE                                      ~ NA_character_
  ))

# Define cohort-specific color palettes
# Each cohort has light and dark shades for the two timepoints
cohort_shades <- list(
  control    = c("#FFF8B0", "#FFD700"),  # Light yellow / gold
  short_term = c("#AED7FF", "#1E90FF"),  # Light blue / dodger blue
  long_term  = c("#FFAD99", "#FF4500")   # Light salmon / orange red
)

cohort_order <- c("control", "short_term", "long_term")

################################################################################
# SECTION 5: CELL TYPE MAPPINGS
################################################################################

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
  "GZMK_Effector_Memory_CD8"  = c(13),
  "Proliferating_Effector"    = c(14, 16, 17),
  "All_CD8"                   = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)

celltypes <- names(celltype_cluster_map)

# Define diversity indices and timepoint comparisons to analyze
diversity_indices  <- c("shannon", "simpson")
timepoint_pairs    <- list(
  Pre_vs_C1 = c("S", "C1"),   # Baseline vs Cycle 1
  Pre_vs_C2 = c("S", "C2"),   # Baseline vs Cycle 2
  C1_vs_C2  = c("C1", "C2")   # Cycle 1 vs Cycle 2
)

# Mapping of timepoint codes to display labels
code2label <- c(S = "Pre", C1 = "C1", C2 = "C2")

# Set up directories
base_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature"
output_root <- file.path(base_dir, "diversity_boxplots")
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

################################################################################
# SECTION 6: HELPER FUNCTIONS FOR PAIRED BOXPLOT VISUALIZATION
################################################################################

################################################################################
# FUNCTION: parse_diversity_file
################################################################################
# Parse diversity index files and extract patient/timepoint information
#
# PARAMETERS:
#   fp - File path to diversity index file
#
# RETURNS:
#   Dataframe with columns: sample_id, diversity, patient_id, tp_code
#
# FORMAT:
#   Sample IDs are formatted as "{patient_id}{timepoint_code}"
#   Example: "5C1" = Patient 5, Cycle 1
#
################################################################################
parse_diversity_file <- function(fp) {
  fread(fp, header = FALSE, sep = "\t", col.names = c("sample_id", "diversity")) %>%
    filter(grepl("^[0-9]", sample_id)) %>%  # Keep only rows starting with numbers
    mutate(patient_id = as.integer(gsub("^([0-9]+).*", "\\1", sample_id)),  # Extract patient ID
           tp_code    = gsub("^[0-9]+", "", sample_id))                      # Extract timepoint code
}

################################################################################
# FUNCTION: plot_pair_for_cohort
################################################################################
# Create paired boxplot for a single cohort comparing two timepoints
#
# PARAMETERS:
#   df_pair - Dataframe with diversity values for two timepoints
#   cohort_name - Name of the cohort ("control", "short_term", "long_term")
#   tp_codes - Vector of two timepoint codes to compare
#   y_lim - Y-axis limits (consistent across cohorts)
#
# RETURNS:
#   ggplot object or NULL if insufficient data
#
# VISUALIZATION:
#   - Boxplots for each timepoint
#   - Lines connecting paired observations (same patient)
#   - Jittered points for individual samples
#   - Cohort-specific color scheme
#   - Paired t-test p-value in title
#
################################################################################
plot_pair_for_cohort <- function(df_pair, cohort_name, tp_codes, y_lim) {
  # Filter for this cohort and ensure paired data
  df_arm <- df_pair %>%
    filter(group == cohort_name) %>%
    group_by(patient_id) %>% 
    filter(n_distinct(tp_code) == 2) %>%  # Keep only patients with both timepoints
    ungroup()
  
  if (nrow(df_arm) == 0) return(NULL)
  
  # Paired t-test: compare diversity at two timepoints within same patients
  wide   <- df_arm %>% 
    select(patient_id, tp_code, diversity) %>%
    pivot_wider(names_from = tp_code, values_from = diversity) %>% 
    na.omit()
  
  p_val  <- t.test(wide[[tp_codes[1]]], wide[[tp_codes[2]]], paired = TRUE)$p.value
  p_lab  <- sprintf("paired t p = %.3g%s", p_val, ifelse(p_val <= 0.05, " *", ""))
  
  # Create timepoint labels
  df_arm <- df_arm %>% 
    mutate(tp_label = factor(code2label[tp_code], levels = code2label[tp_codes]))
  
  # Apply cohort-specific colors
  fills  <- setNames(cohort_shades[[cohort_name]], code2label[tp_codes])
  
  # Create paired boxplot
  ggplot(df_arm, aes(x = tp_label, y = diversity, fill = tp_label)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, colour = "black", position = position_dodge(width = 1.2)) +
    geom_line(aes(group = patient_id), colour = "black") +  # Connect paired observations
    geom_jitter(position = position_dodge(width = 1.2), colour = "black") +
    scale_fill_manual(values = fills) +
    scale_y_continuous(limits = y_lim, name = "Clonal diversity") +
    labs(x = "Timepoint", y = "Clonal diversity", title = sprintf("%s (%s)", cohort_name, p_lab)) +
    theme_minimal(base_size = 14) +
    theme(axis.line       = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(hjust = 0.5, size = 15),
          axis.title.x     = element_text(size = 16),
          axis.title.y     = element_text(size = 16),
          axis.text.x      = element_text(size = 14),
          axis.text.y      = element_text(size = 14),
          legend.position  = "none")
}

################################################################################
# SECTION 7: GENERATE EXPANSION RATIO PLOTS FOR ALL ARMS COMBINED
################################################################################

# Function to create combined visualization across all treatment arms
plot_comparison_all_arms <- function(df_merged,
                                     variable_name = "clonal_expansion",
                                     diversity_index = "shannon",
                                     alternative = "l",
                                     plot_dir,
                                     celltype,
                                     control_arm = "control") {
  
  # Filter out rows with missing clonal expansion values
  df_merged <- df_merged[!is.na(df_merged[[variable_name]]), ]
  
  # Get control arm data for comparison
  control_df <- df_merged[df_merged$group == control_arm, ]
  
  # Determine y-axis limits that accommodate all arms
  y_min <- min(df_merged[[variable_name]], na.rm = TRUE)
  y_max <- max(df_merged[[variable_name]], na.rm = TRUE)
  
  # Store individual plots for each arm
  plot_list <- list()
  unique_arms <- unique(df_merged$group)
  
  for (arm in unique_arms) {
    # Subset data for current arm
    df_arm <- df_merged[df_merged$group == arm, ]
    finite_values <- df_arm[[variable_name]][is.finite(df_arm[[variable_name]])]
    
    # Calculate p-value comparing to control arm
    p_val <- NA
    p_label <- ""
    col <- "black"
    
    if (arm != control_arm && nrow(control_df) > 1 && length(finite_values) > 1) {
      control_finite <- control_df[[variable_name]][is.finite(control_df[[variable_name]])]
      if (length(control_finite) > 1) {
        # Two-sample t-test vs control
        p_val <- t.test(x = finite_values, 
                        y = control_finite,
                        alternative = alternative)[["p.value"]]
        
        p_label <- paste("p =", round(p_val, 3))
        if (p_val <= 0.05) {
          p_label <- paste0(p_label, " *")
          col <- "red"
        }
      }
    }
    
    # Create plot title
    if (arm == control_arm) {
      plot_title <- paste0(arm, " (No p-value)")
    } else {
      plot_title <- paste0(arm, if (p_label != "") paste0(" (", p_label, ")"))
    }
    
    # Build plot for this arm
    p <- ggplot(df_arm, aes(x = arm, y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.3) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +   # Remove outlier dots to avoid overlap
      geom_point(position = position_jitter(width = 0.1)) +
      scale_y_continuous(name = variable_name, limits = c(y_min, y_max)) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_title) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 1) +  # Reference line at 1.0
      theme_bw() +
      theme(
        legend.position = "none",
        text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, color = col)
      )
    
    plot_list[[arm]] <- p
  }
  
  # Combine all arm plots into a single figure
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
  
  # Save combined plot
  file_name <- paste0(diversity_index, "_Clonal_Diversity_", celltype, "_allArms_", variable_name, ".pdf")
  pdf(file = file.path(plot_dir, file_name), width = 12, height = 5)
  print(combined_plot)
  dev.off()
}

################################################################################
# SECTION 8: GENERATE COMBINED PLOTS (C1 vs C2)
################################################################################

dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2"

for (celltype in celltypes) {
  cat("Processing celltype:", celltype, "\n")
  
  for (diversity_index in diversity_indices) {
    cat("  Diversity index:", diversity_index, "\n")
    
    # Build output directory
    plot_dir <- paste0(dir, "/clonal_expansion_", diversity_index, "_", celltype, "_T_cells_division/")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    
    # Load expansion data
    data_file <- file.path(plot_dir, "C1_vs_C2.txt")
    if (!file.exists(data_file)) {
      cat("  File not found:", data_file, "\n")
      next
    }
    df <- read.delim(data_file)
    
    # Merge with clinical metadata
    df_merged <- merge(df, 
                       patient_metadata_df[, c("patient_id", "Arm", "group")], 
                       by.x = "row.names", 
                       by.y = "patient_id", 
                       all.x = TRUE, 
                       all.y = FALSE)
    colnames(df_merged)[1] <- "patient_id"
    
    # Skip if insufficient data
    if (nrow(df_merged) < 2) {
      cat("   Not enough data. Skipping.\n")
      next
    }
    
    # Generate combined figure for all arms
    plot_comparison_all_arms(df_merged,
                             variable_name = "clonal_expansion",
                             diversity_index = diversity_index,
                             alternative = "l",
                             plot_dir,
                             celltype,
                             control_arm = "control")
  }
}

################################################################################
# SECTION 9: GENERATE PAIRED BOXPLOTS ACROSS ALL TIMEPOINT COMPARISONS
################################################################################
# Main loop to create comprehensive paired diversity boxplots

for (celltype in celltypes) {
  message("\nProcessing cell-type: ", celltype)
  
  for (div_idx in diversity_indices) {
    message("  Diversity index: ", div_idx)
    
    # Locate diversity file from Step 2
    fp <- file.path(base_dir, sprintf("%s_clonal_diversity_%s_T_cells.txt", div_idx, celltype))
    if (!file.exists(fp)) {
      message("     * File not found → skipping *")
      next
    }
    
    # Parse diversity data and merge with metadata
    raw_df   <- parse_diversity_file(fp)
    df_merge <- raw_df %>% 
      inner_join(patient_metadata_df %>% select(patient_id, Arm, group), by = "patient_id")
    
    if (nrow(df_merge) < 4) next  # Need minimum data for meaningful plots
    
    # Generate plots for each timepoint comparison
    for (pair_name in names(timepoint_pairs)) {
      tp_codes <- timepoint_pairs[[pair_name]]
      df_pair  <- df_merge %>% filter(tp_code %in% tp_codes)
      
      if (nrow(df_pair) < 4) next
      
      # Determine y-axis limits for consistent scaling across cohorts
      y_lim <- range(df_pair$diversity, na.rm = TRUE)
      
      # Generate plot for each cohort
      plot_list <- list()
      for (cohort in cohort_order) {
        if (cohort %in% df_pair$group) {
          gg <- plot_pair_for_cohort(df_pair, cohort, tp_codes, y_lim)
          if (!is.null(gg)) plot_list[[cohort]] <- gg
        }
      }
      
      if (length(plot_list) == 0) next
      
      # Combine cohort plots into single figure
      combined <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
      
      # Save to organized directory structure
      out_dir  <- file.path(output_root, pair_name, paste0(div_idx, "_", celltype))
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      
      pdf(file = file.path(out_dir, sprintf("%s_%s_%s_boxplot.pdf", pair_name, div_idx, celltype)), 
          width = 12, height = 5)
      print(combined)
      dev.off()
    }
  }
}

################################################################################
# END OF SCRIPT
################################################################################
# 
# OUTPUT SUMMARY:
# This script generates two types of visualizations:
#
# 1. EXPANSION RATIO PLOTS:
#    - Individual PDFs per treatment arm
#    - Combined PDFs with all arms side-by-side
#    - T-test p-values comparing to no-change (ratio = 1.0) or to control
#
# 2. PAIRED DIVERSITY BOXPLOTS (FIGURE 6d):
#    - Organized by timepoint comparison (Pre vs C1, Pre vs C2, C1 vs C2)
#    - Faceted by patient cohort (control, short_term, long_term)
#    - Paired lines showing individual patient trajectories
#    - Paired t-test p-values within each cohort
#    - Consistent y-axis scaling for fair comparison
#
# These visualizations reveal:
#   - Which T cell subpopulations undergo clonal expansion
#   - When expansion occurs (early vs late treatment)
#   - Differences between responder cohorts
#   - Individual patient heterogeneity
#
################################################################################
