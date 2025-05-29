library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# ==================================================================================================================================================================
# cluster proportion change comparison

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


seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("NK", "B cells", "B cells", "B cells", "B cells", "cDC", "pDC", "Classical Monocytes", "Intermediate Monocytes", "Intermediate Monocytes", "Non-Classical Monocytes"),
  cluster = c(2, 3, 16, 20, 21, 18, 23, 0, 5, 12, 6)
)

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Create boxplot
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
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
  file_name <- paste0(plot_dir, celltype, "_All_Cells_Cluster_Proportion_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6, limitsize = FALSE)
}



# ==================================================================================================================================================================
# Monocyte Subpopulations cluster proportion change comparison in All Monocytes

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_object_monocyte_cells <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(0, 5, 12, 6))
seurat_metadata <- seurat_object_monocyte_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Classical_Monocytes", "Intermediate_Monocytes", "Intermediate_Monocytes", "Non_Classical_Monocytes"),
  cluster = c(0, 5, 12, 6)
)

timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/monocytes/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/monocytes/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cells_Cluster_Proportion_in_Monocyte_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}





# ==================================================================================================================================================================
# Monocyte Subpopulations cluster proportion change comparison in All Monocytes (UF Only patients)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_object_UF_cells <- subset(seurat_object_all_cells, subset = Site == "UF")
seurat_object_monocyte_cells <- subset(seurat_object_UF_cells, subset = seurat_clusters %in% c(0, 5, 12, 6))
seurat_metadata <- seurat_object_monocyte_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "ratio"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Classical_Monocytes", "Intermediate_Monocytes", "Intermediate_Monocytes", "Non_Classical_Monocytes"),
  cluster = c(0, 5, 12, 6)
)

timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/monocytes/difference/UF_only/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/monocytes/ratio/UF_only/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cells_Cluster_Proportion_in_Monocyte_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}





# =============================================================================================================================================================================================
# For T Cells only (proportion in All Cells)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_seurat_object_GBM.rds")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Activated_CD4", "Effector_CD8", "Memory_CD4", "Naïve_CD4", "Memory_CD8", "Foxp3_CTLA4_CD4", "Proliferating_memory", "Proliferating_memory", "Proliferating_effectors", "Activated_CD8", "CD4", "CD4", "CD4", "CD4", "CD4", "CD8", "CD8", "CD8", "CD8"),
  cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 0, 2, 3, 5, 10, 1, 4, 7, 8)
)

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  #######################################################################
  tcell_barcodes <- rownames(seurat_metadata_t_cells[seurat_metadata_t_cells[[cluster_col]] %in% cluster_list, ])
  temp <- seurat_metadata[tcell_barcodes,]
  tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
  #######################################################################
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Create boxplot
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
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
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6, limitsize = FALSE)
}


# =============================================================================================================================================================================================
# For sub populations of T Cells only (proportion in T Cells) 

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

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

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  tcell_cluster_list <- cluster_list
  seurat_metadata <- seurat_metadata_t_cells
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/celltype/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/celltype/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}



# =============================================================================================================================================================================================
# For sub populations of T Cells only (proportion in T Cells) (UF Patients ONLY)

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

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

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  tcell_cluster_list <- cluster_list
  seurat_metadata <- seurat_metadata_t_cells
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/celltype/difference/UF_only/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/celltype/ratio/UF_only/"
  }
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}




# =============================================================================================================================================================================================
# For all clusters in T Cells only (proportion in T Cells) 

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "ratio"


timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (cluster in unique(seurat_metadata_t_cells$seurat_clusters)){
  cluster_list <- as.vector(cluster)
  celltype <- paste0("Cluster_", cluster)
  tcell_cluster_list <- cluster_list
  seurat_metadata <- seurat_metadata_t_cells
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/cluster_number/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/cluster_number/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}





# =============================================================================================================================================================================================
# For all clusters in T Cells only (proportion in T Cells) ((UF Patients ONLY))

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"


timepoints <- c("C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

for (cluster in unique(seurat_metadata_t_cells$seurat_clusters)){
  cluster_list <- as.vector(cluster)
  celltype <- paste0("Cluster_", cluster)
  tcell_cluster_list <- cluster_list
  seurat_metadata <- seurat_metadata_t_cells
  
  # Initialize a list to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Remove the 'pre_' prefix from column names
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
    
    # Create the plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Arm",
           x = "Time Point",
           y = "Cluster Proportion") +
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
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = "Arm",
           y = "Proportion Change") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p
  }
  
  # ordered_list <- list(plot_list[[1]], plot_list[[3]], plot_list[[5]], plot_list[[2]], plot_list[[4]], plot_list[[6]])
  plot_list <- c(plot_list_1, plot_list_2)
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/cluster_number/difference/UF_only/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/t_cells/cluster_number/ratio/UF_only/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}


# ==============================================================================================================================================================================================
# cluster proportion comparison (not change)

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


seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("NK", "B cells", "B cells", "B cells", "B cells", "cDC", "pDC", "Classical Monocytes", "Intermediate Monocytes", "Intermediate Monocytes", "Non-Classical Monocytes"),
  cluster = c(2, 3, 16, 20, 21, 18, 23, 0, 5, 12, 6)
)

timepoints <- c("Pre", "C1", "C2")

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  for (timepoint in timepoints){
    timepoint_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint)
    colnames(timepoint_cluster_proportion_df) <- c("patient_id", "cluster_proportion")
    merged_proportion_df <- merge(timepoint_cluster_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(cluster_proportion ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Create boxplot
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = cluster_proportion, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Cluster Propoprtion at ", timepoint),
           x = "Arm",
           y = "Proportion") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$cluster_proportion, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list[[timepoint]] <- p
  }
  
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
  plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/"
  file_name <- paste0(plot_dir, celltype, "_All_Cells_Cluster_Proportion_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6, limitsize = FALSE)
}

# =============================================================================================================================================================================================
# For T Cells only

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_seurat_object_GBM.rds")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"


# Example of a cell type to cluster mapping
celltype_to_cluster <- data.frame(
  celltype = c("Activated CD4", "Effector CD8", "Memory CD4", "Naïve CD4", "Memory CD8", "Foxp3 CTLA4 CD4", "Proliferating memory", "Proliferating memory", "Proliferating effectors", "Activated CD8", "CD4"),
  cluster = c(0, 1, 2, 3, 4, 5, 6, 12, 7, 8, 10)
)

timepoints <- c("Pre", "C1", "C2")

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  #######################################################################
  tcell_barcodes <- rownames(seurat_metadata_t_cells[seurat_metadata_t_cells[[cluster_col]] %in% cluster_list, ])
  temp <- seurat_metadata[tcell_barcodes,]
  tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
  #######################################################################
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  for (timepoint in timepoints){
    timepoint_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint)
    colnames(timepoint_cluster_proportion_df) <- c("patient_id", "cluster_proportion")
    merged_proportion_df <- merge(timepoint_cluster_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      # Perform Wilcoxon test
      test_result <- wilcox.test(cluster_proportion ~ Arm, data = merged_proportion_df)
      p_value = test_result$p.value
    } else {
      p_value = NA_real_  # Use NA_real_ for numerical context missing values
    }
    
    # Create boxplot
    p <- ggplot(merged_proportion_df, aes(x = Arm, y = cluster_proportion, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Cluster Propoprtion at ", timepoint),
           x = "Arm",
           y = "Proportion") +
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
      annotate("text", x = 1.5, y = max(merged_proportion_df$cluster_proportion, na.rm = TRUE), 
               label = paste("p-value =", round(p_value, 4)), size = 4, vjust = -0.5)
    
    # Add the plot to the list
    plot_list[[timepoint]] <- p
  }
  
  # Combine all patient plots into a grid
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), align = 'v')
  plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/cluster_proportion_plots/"
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6, limitsize = FALSE)
}




seurat_object_all_cells@meta.data$TimePoint <- factor(seurat_object_all_cells@meta.data$TimePoint, levels = c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36"))
seurat_metadata <- seurat_object_all_cells@meta.data
control_cells <- rownames(seurat_metadata[seurat_metadata$Patient %in% c("8", "11", "9", "4"), ])
exp_cells <- rownames(seurat_metadata[seurat_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7"), ])

temp1 <- subset(seurat_object_all_cells, cells = control_cells)
temp2 <- subset(seurat_object_all_cells, cells = exp_cells)

DimPlot(temp1, split.by = "TimePoint", raster = FALSE)
DimPlot(temp2, split.by = "TimePoint", raster = FALSE)
