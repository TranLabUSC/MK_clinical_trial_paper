library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# ==================================================================================================================================================================
# cluster proportion change comparison
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
    
    # Check for cells at timepoint N
    baseline_timepoint <- intersect(timepoints, unique(patient_seurat_obj@meta.data$TimePoint))[1]
    timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == baseline_timepoint, ])
    
    if (length(timepoint_N_barcodes) > 0) {
      # Calculate mean and standard deviation for cells at timepoint N
      mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
      sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])
      print(sd_N)
      if (sd_N == 0 | is.na(sd_N)){
        scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
      } else {
        # Scale based on timepoint N
        scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
      }
    } else {
      # Use scale function if no cells at timepoint N
      scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    }
    
    pathway_avg_expression <- scaled_pathway_avg_expression
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




get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for the cluster of interest and the specified timepoint
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







# =============================================================================================================================================================================================
# For T Cells only (cluster proportion in all cells)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_seurat_object_GBM.rds")
# seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- FALSE
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

method = "ratio"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Activated_CD4", "Effector_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Activated_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8"),
  cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8)
)

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
pvalue_list <- c()

for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    #######################################################################
    tcell_barcodes <- rownames(seurat_object_t_cell_subset@meta.data)
    temp <- seurat_metadata[tcell_barcodes,]
    tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
    #######################################################################
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
    updated_survival_df <- updated_survival_df[updated_survival_df$IDH != "POS",]
    print(updated_survival_df)
    
    # Initialize a list to store the plots
    plot_list <- list()
    
    for (timepoint_pair in pairs_list){
      temp_survival_df <- updated_survival_df
      timepoint_a <- timepoint_pair[1]
      timepoint_b <- timepoint_pair[2]
      
      if (include_cluster_proportion) {
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]) / (temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)])
        }
      } else {
        temp_survival_df[, "signal_change"] <- (temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] - temp_survival_df[, paste0("Mean_Expr_", timepoint_a)])
      }
      
      merged_proportion_df <- temp_survival_df
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      # Define custom colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      
      if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
          sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
          sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
        # Perform Wilcoxon test
        test_result <- wilcox.test(signal_change ~ Arm, data = merged_proportion_df)
        p_value = test_result$p.value
      } else {
        p_value = NA_real_  # Use NA_real_ for numerical context missing values
      }
      
      pathway_list <- c(pathway_list, pathway)
      celltype_list <- c(celltype_list, celltype)
      comparison_list <- c(comparison_list, paste(timepoint_a, timepoint_b, sep = "_vs_"))
      pvalue_list <- c(pvalue_list, p_value)
      
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
      plot_list[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
    }
    
    # Combine all patient plots into a grid
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    if (method == "difference") {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/difference/"
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/ratio/"
    }
    if(include_cluster_proportion){
      file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison_With_Cluster_Proportion.pdf")
    } else {
      file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison.pdf")
    }
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, pvalue = pvalue_list)
write.csv(df, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/ratio/signal_comparison_without_cluster_proportion.csv")


# =============================================================================================================================================================================================
# For subpopulations of T Cells only (cluster proportion in T cells)

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- FALSE
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

method = "ratio"

# # Example of a cell type to cluster mapping
# celltype_to_cluster <- data.frame(
#   celltype = c("Activated_CD4", "Effector_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Activated_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8", "Effector+Activated_CD8", "Effector+Activated_CD8", "Proliferating+Activated_CD8", "Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Memory+Activated_CD8", "Memory+Activated_CD8"),
#   cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8, 1, 8, 7, 8, 1, 7, 8, 4, 8)
# )

# celltype_to_cluster <- data.frame(
#   celltype = c("Activated_CD4", "Transitional_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Effector_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8", "Transitional+Effector_CD8", "Transitional+Effector_CD8", "Proliferating+Effector_CD8", "Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Memory+Effector_CD8", "Memory+Effector_CD8"),
#   cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8, 1, 8, 7, 8, 1, 7, 8, 4, 8)
# )

# Define the mapping
mapping <- list(
  "Activated_CD4" = c(0, 1, 11, 14, 22),
  "Effector_CD8" = c(2),
  "Memory_CD4" = c(3, 4, 7),
  "Naive_CD4" = c(5),
  "Transitional_CD8" = c(6, 16),
  "Memory_CD8" = c(8),
  "Naive_CD8" = c(9),
  "Energic_CD8" = c(10),
  "Transitional_CD4" = c(12),
  "Exhausted_T" = c(15),
  "Proliferative_effector_CD8" = c(17, 18),
  "Proliferative_memory_CD8" = c(21),
  "Effector+Proliferative_effector_CD8" = c(2, 17, 18),
  "Effector+Energic_CD8" = c(2, 10)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))


timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()

for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    tcell_cluster_list <- cluster_list
    seurat_metadata <- seurat_metadata_t_cells
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
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
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      } else {
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
        temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
      }
      
      merged_proportion_df <- temp_survival_df
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      # Define custom colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      
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
            log_fold_change <- log2(treatment_median - control_median)
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
        
        
        # Perform Wilcoxon test
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
    
    # Combine all patient plots into a grid
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    if (method == "difference") {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/difference/new_nomenclature/"
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/ratio/new_nomenclature/"
    }
    if(include_cluster_proportion){
      file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison_With_Cluster_Proportion_in_T_Cells.pdf")
    } else {
      file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison.pdf")
    }
    
    if(include_cluster_proportion){
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/ratio/"
      }
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/without_cluster_proportion/"
    }
    file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison.pdf")
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, FC = fc_list, logFC = logfc_list, pvalue = pvalue_list)
if(include_cluster_proportion){
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/difference/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/ratio/"
  }
} else {
  dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/without_cluster_proportion/"
}
file_name <- paste0(dir, "signal_comparison_in_T_cells.csv")
write.csv(df, file_name)



# =============================================================================================================================================================================================
# For subpopulations of T Cells only (cluster proportion in T cells) (UF Patients ONLY)

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- TRUE
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

method = "ratio"

# # Example of a cell type to cluster mapping
# celltype_to_cluster <- data.frame(
#   celltype = c("Activated_CD4", "Effector_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Activated_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8", "Effector+Activated_CD8", "Effector+Activated_CD8", "Proliferating+Activated_CD8", "Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Effector+Proliferating+Activated_CD8", "Memory+Activated_CD8", "Memory+Activated_CD8"),
#   cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8, 1, 8, 7, 8, 1, 7, 8, 4, 8)
# )

# celltype_to_cluster <- data.frame(
#   celltype = c("Activated_CD4", "Transitional_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Effector_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8", "Transitional+Effector_CD8", "Transitional+Effector_CD8", "Proliferating+Effector_CD8", "Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Transitional+Proliferating+Effector_CD8", "Memory+Effector_CD8", "Memory+Effector_CD8"),
#   cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8, 1, 8, 7, 8, 1, 7, 8, 4, 8)
# )

# Define the mapping
mapping <- list(
  "Activated_CD4" = c(0, 1, 11, 14, 22),
  "Effector_CD8" = c(2),
  "Memory_CD4" = c(3, 4, 7),
  "Naive_CD4" = c(5),
  "Transitional_CD8" = c(6, 16),
  "Memory_CD8" = c(8),
  "Naive_CD8" = c(9),
  "Energic_CD8" = c(10),
  "Transitional_CD4" = c(12),
  "Exhausted_T" = c(15),
  "Proliferative_effector_CD8" = c(17, 18),
  "Proliferative_memory_CD8" = c(21),
  "Effector+Proliferative_effector_CD8" = c(2, 17, 18),
  "Effector+Energic_CD8" = c(2, 10)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))


timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()

for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    tcell_cluster_list <- cluster_list
    seurat_metadata <- seurat_metadata_t_cells
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
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
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      } else {
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
        temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
      }
      
      merged_proportion_df <- temp_survival_df
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      # Define custom colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      
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
            log_fold_change <- log2(treatment_median - control_median)
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
        
        
        # Perform Wilcoxon test
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
    
    # Combine all patient plots into a grid
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    
    if(include_cluster_proportion){
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/difference/UF_only/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/ratio/UF_only/"
      }
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/without_cluster_proportion/UF_only/"
    }
    file_name <- paste0(plot_dir, celltype, "_T_Cells_", pathway, "_Signal_Change_Comparison.pdf")
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, FC = fc_list, logFC = logfc_list, pvalue = pvalue_list)
if(include_cluster_proportion){
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/difference/UF_only/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/with_cluster_proportion/ratio/UF_only/"
  }
} else {
  dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/t_cells/celltype/without_cluster_proportion/UF_only/"
}
file_name <- paste0(dir, "signal_comparison_in_T_cells.csv")
write.csv(df, file_name)





# =============================================================================================================================================================================================
# For subpopulations of Monocyte Cells only (cluster proportion in Monocyte cells)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_object_monocyte_cells <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(0, 5, 12, 6))
seurat_metadata <- seurat_object_monocyte_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- TRUE
gmt_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/monocyte_pathways.gmt"

pathways <- c("GO_CCR2_CHEMOKINE_RECEPTOR_BINDING", "GO_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN", "GO_POSITIVE_REGULATION_OF_PHAGOCYTOSIS_ENGULFMENT", "GO_TOLL_SIGNALING_PATHWAY", "GO_RESPIRATORY_BURST_AFTER_PHAGOCYTOSIS", "GO_REGULATION_OF_TOLL-LIKE_RECEPTOR_SIGNALING_PATHWAY", "GO_REGULATION_OF_NITRIC_OXIDE_MEDIATED_SIGNAL_TRANSDUCTION", "GO_ANTIBODY-DEPENDENT_CELLULAR_CYTOTOXICITY", "GO_MONOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GO_PHOSPHATIDYLINOSITOL_3-KINASE_REGULATOR_ACTIVITY", "GO_TOLL-LIKE_RECEPTOR_BINDING", "GO_CD80_BIOSYNTHETIC_PROCESS", "GO_ARGINASE_ACTIVITY", "GO_ICAM-3_RECEPTOR_ACTIVITY", "GO_MONOCYTE_DIFFERENTIATION", "GO_REGULATION_OF_MONOCYTE_CHEMOTAXIS", "GO_POSITIVE_REGULATION_OF_IMMUNE_COMPLEX_CLEARANCE_BY_MONOCYTES_AND_MACROPHAGES", "GO_FC-GAMMA_RECEPTOR_SIGNALING_PATHWAY_INVOLVED_IN_PHAGOCYTOSIS", "GO_MYELOID_LEUKOCYTE_MIGRATION", "GO_MONOCYTE_AGGREGATION", "GO_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_ACTIVITY")

method = "ratio"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Classical_Monocytes", "Intermediate_Monocytes", "Intermediate_Monocytes", "Non_Classical_Monocytes"),
  cluster = c(0, 5, 12, 6)
)

timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()

for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    seurat_object_subset <- subset(seurat_object_all_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_subset <- NormalizeData(seurat_object_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_subset, survival_data)
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
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      } else {
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
        temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
      }
      
      merged_proportion_df <- temp_survival_df
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      # Define custom colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      
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
            log_fold_change <- log2(treatment_median - control_median)
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
        
        
        # Perform Wilcoxon test
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
    
    # Combine all patient plots into a grid
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    if(include_cluster_proportion){
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/difference/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/ratio/"
      }
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/without_cluster_proportion/"
    }
    file_name <- paste0(plot_dir, celltype, "_", pathway, "_Signal_Change_Comparison.pdf")
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, FC = fc_list, logFC = logfc_list, pvalue = pvalue_list)
if(include_cluster_proportion){
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/difference/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/ratio/"
  }
} else {
  dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/without_cluster_proportion/"
}
file_name <- paste0(dir, "signal_comparison_in_monocyte_cells.csv")
write.csv(df, file_name)





# =============================================================================================================================================================================================
# For subpopulations of Monocyte Cells only (cluster proportion in Monocyte cells) (UF Patients ONLY)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_object_UF_cells <- subset(seurat_object_all_cells, subset = Site == "UF")
seurat_object_monocyte_cells <- subset(seurat_object_UF_cells, subset = seurat_clusters %in% c(0, 5, 12, 6))
seurat_metadata <- seurat_object_monocyte_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- TRUE
gmt_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/monocyte_pathways.gmt"

pathways <- c("GO_CCR2_CHEMOKINE_RECEPTOR_BINDING", "GO_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN", "GO_POSITIVE_REGULATION_OF_PHAGOCYTOSIS_ENGULFMENT", "GO_TOLL_SIGNALING_PATHWAY", "GO_RESPIRATORY_BURST_AFTER_PHAGOCYTOSIS", "GO_REGULATION_OF_TOLL-LIKE_RECEPTOR_SIGNALING_PATHWAY", "GO_REGULATION_OF_NITRIC_OXIDE_MEDIATED_SIGNAL_TRANSDUCTION", "GO_ANTIBODY-DEPENDENT_CELLULAR_CYTOTOXICITY", "GO_MONOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GO_PHOSPHATIDYLINOSITOL_3-KINASE_REGULATOR_ACTIVITY", "GO_TOLL-LIKE_RECEPTOR_BINDING", "GO_CD80_BIOSYNTHETIC_PROCESS", "GO_ARGINASE_ACTIVITY", "GO_ICAM-3_RECEPTOR_ACTIVITY", "GO_MONOCYTE_DIFFERENTIATION", "GO_REGULATION_OF_MONOCYTE_CHEMOTAXIS", "GO_POSITIVE_REGULATION_OF_IMMUNE_COMPLEX_CLEARANCE_BY_MONOCYTES_AND_MACROPHAGES", "GO_FC-GAMMA_RECEPTOR_SIGNALING_PATHWAY_INVOLVED_IN_PHAGOCYTOSIS", "GO_MYELOID_LEUKOCYTE_MIGRATION", "GO_MONOCYTE_AGGREGATION", "GO_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_ACTIVITY")

method = "ratio"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Classical_Monocytes", "Intermediate_Monocytes", "Intermediate_Monocytes", "Non_Classical_Monocytes"),
  cluster = c(0, 5, 12, 6)
)

timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

pathway_list <- c()
celltype_list <- c()
comparison_list <- c()
fc_list <- c()
logfc_list <- c()
pvalue_list <- c()

for (pathway in pathways){
  for (celltype in unique(celltype_to_cluster$celltype)){
    print(celltype)
    cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
    seurat_object_subset <- subset(seurat_object_all_cells, subset = seurat_clusters %in% cluster_list)
    seurat_object_subset <- NormalizeData(seurat_object_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_subset, survival_data)
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
        timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
        timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
        temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
        temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
        if (method == "difference") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        } else if (method == "ratio") {
          temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
        }
      } else {
        temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
        temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
        temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
      }
      
      merged_proportion_df <- temp_survival_df
      merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
      # Define custom colors
      custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
      
      
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
            log_fold_change <- log2(treatment_median - control_median)
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
        
        
        # Perform Wilcoxon test
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
    
    # Combine all patient plots into a grid
    plot_list <- c(plot_list_1, plot_list_2)
    final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
    if(include_cluster_proportion){
      if (method == "difference") {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/difference/UF_only/"
      } else {
        plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/ratio/UF_only/"
      }
    } else {
      plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/without_cluster_proportion/UF_only/"
    }
    file_name <- paste0(plot_dir, celltype, "_", pathway, "_Signal_Change_Comparison.pdf")
    # Save the grid plot as a PDF
    ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
  }
}

df = data.frame(pathway = pathway_list, celltype = celltype_list, comparison = comparison_list, FC = fc_list, logFC = logfc_list, pvalue = pvalue_list)
if(include_cluster_proportion){
  if (method == "difference") {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/difference/UF_only/"
  } else {
    dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/with_cluster_proportion/ratio/UF_only/"
  }
} else {
  dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/pathway_activity_plots/monocytes/without_cluster_proportion/UF_only/"
}
file_name <- paste0(dir, "signal_comparison_in_monocyte_cells.csv")
write.csv(df, file_name)