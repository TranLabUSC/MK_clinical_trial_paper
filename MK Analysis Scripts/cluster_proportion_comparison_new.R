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



survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "ratio"

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

seurat_object <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
seurat_metadata <- seurat_object@meta.data


# Define the mapping
mapping <- list(
  "Classical_Monocytes" = c(1, 3, 16, 12),
  "Non_Classical_Monocytes" = c(9,33),
  "cDC" = c(31),
  "pDC" = c(36),
  "NK" = c(0)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  seurat_metadata <- seurat_metadata
  
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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_All_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}




# ==================================================================================================================================================================
# Monocyte Subpopulations cluster proportion change comparison in All Monocytes

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
seurat_object_monocyte_cells <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(1, 3, 6, 12, 9, 33))
seurat_metadata <- seurat_object_monocyte_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

# Define the mapping
mapping <- list(
  "Classical_Monocytes" = c(1, 3, 16, 12),
  "Non_Classical_Monocytes" = c(9,33)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_Monocyte_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}





######################################################################################################################################################################################################################
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data


# Define the mapping
mapping <- list(
  "Activated_CD4" = c(0),
  "Effector_CD8" = c(1),
  "Effector_Memory_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  "Activated_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  "Memory_CD8" = c(8),
  "Anergic_CD8" = c(10),
  "Naive_CD8" = c(12),
  "Hyperactivated_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17),
  "CD8" = c(1, 2, 13, 14, 16, 17)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}


#############################################################################################################################################################################################################################
# Load Necessary Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# ============================================================================================================
# Cluster Proportion Change Comparison

# Enhanced Function to Calculate Cluster Proportions with Error Handling
get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Check if required columns exist
  required_cols <- c(patient_col, timepoint_col, cluster_col)
  missing_cols <- setdiff(required_cols, colnames(seurat_metadata))
  if(length(missing_cols) > 0){
    stop(paste("The following required columns are missing in seurat_metadata:", paste(missing_cols, collapse = ", ")))
  }
  
  # Filter metadata for the cluster of interest and the specified timepoint
  cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest & seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  timepoint_data <- seurat_metadata[seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  
  # Check if there is any data after filtering
  if(nrow(timepoint_data) == 0){
    warning(paste("No data found for timepoint:", timepoint_of_interest))
    return(data.frame())
  }
  
  # Group the filtered data by patient
  cluster_timpoint_grouped_data <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
  timpoint_grouped_data <- split(timepoint_data, timepoint_data[[patient_col]])
  
  # Initialize lists to store patient_id and Cluster_Proportion_X
  subject_ids <- c()
  cluster_proportions <- c()
  
  # Loop through each patient
  for (key in names(timpoint_grouped_data)) {
    patient_timepoint_data <- timpoint_grouped_data[[key]]
    patient_cluster_timepoint_data <- cluster_timpoint_grouped_data[[key]]
    
    # Calculate the total number of cells for the patient
    patient_timepoint_cells <- nrow(patient_timepoint_data)
    patient_cluster_timepoint_cells <- ifelse(is.null(patient_cluster_timepoint_data), 0, nrow(patient_cluster_timepoint_data))
    
    # Calculate the cluster proportion if total cells > 0
    if (patient_timepoint_cells > 0) {
      cluster_proportion <- patient_cluster_timepoint_cells / patient_timepoint_cells
      
      # Append the patient_id and Cluster_Proportion_X to the lists
      subject_ids <- c(subject_ids, key)
      cluster_proportions <- c(cluster_proportions, cluster_proportion)
    }
  }
  
  # If no subjects have cluster proportions, return empty dataframe
  if(length(subject_ids) == 0){
    warning(paste("No cluster proportions calculated for timepoint:", timepoint_of_interest))
    return(data.frame())
  }
  
  custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
  
  # Create a data frame with patient_id and custom Cluster_Proportion_X column
  result_df <- data.frame(patient_id = subject_ids)
  result_df[[custom_column_name]] <- cluster_proportions
  
  return(result_df)
}

# ============================================================================================================
# Define Patient Groups

control_group <- c(1, 4, 8, 9)
experiment_group <- c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Initialize a Data Frame to Store Statistical Summaries
stat_summary <- data.frame()

# ============================================================================================================
# Load Survival Data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF", ]

# Ensure 'patient_id' is Numeric
# Adjust column name as per your CSV structure
if ("patient_id" %in% colnames(survival_df)) {
  survival_df$patient_id <- as.numeric(as.character(survival_df$patient_id))
} else if ("Patient" %in% colnames(survival_df)) {
  survival_df$patient_id <- as.numeric(as.character(survival_df$Patient))
} else {
  stop("Cannot find 'patient_id' or 'Patient' column in survival data.")
}

# ============================================================================================================
# Define Columns and Method

patient_col <- "Patient"       # Updated back to original column name
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "ratio"  # Options: "difference" or "ratio"

# Define Timepoints and Pair Combinations
timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

# ============================================================================================================
# Load Seurat Object and Extract Metadata

seurat_object <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
seurat_metadata <- seurat_object@meta.data

# ============================================================================================================
# Define Cell Type to Cluster Mapping

mapping <- list(
  "Classical_Monocytes" = c(1, 3, 16, 12),
  "Non_Classical_Monocytes" = c(9, 33),
  "cDC" = c(31),
  "pDC" = c(36),
  "NK" = c(0)
)

# Convert the List to a Data Frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# ============================================================================================================
# Iterate Over Each Cell Type

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  
  # Check cluster_list
  print(paste("Processing Cell Type:", celltype))
  print(paste("Clusters of Interest:", paste(cluster_list, collapse = ", ")))
  
  # Initialize Lists to Store Plots
  plot_list_line <- list()
  plot_list_box <- list()
  
  # Iterate Over Each Timepoint Pair
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    print(paste("Comparing Timepoints:", timepoint_a, "vs", timepoint_b))
    
    # Calculate Cluster Proportions for Both Timepoints
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
    
    # Check if data frames are empty
    if(nrow(timepoint_a_cluster_proportion_df) == 0 || nrow(timepoint_b_cluster_proportion_df) == 0){
      warning(paste("Skipping comparison between", timepoint_a, "and", timepoint_b, "due to insufficient data."))
      next  # Skip to next iteration
    }
    
    # Merge Proportion Data Frames
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
    
    # Merge with Survival Data to Get 'Arm' Information
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
    
    # Ensure patient_id is numeric
    merged_proportion_df$patient_id <- as.numeric(as.character(merged_proportion_df$patient_id))
    
    # Calculate 'Proportion Change' Based on Selected Method
    if (method == "difference") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
      # Handle Inf or NaN Values
      merged_proportion_df$proportion_change[!is.finite(merged_proportion_df$proportion_change)] <- NA
    }
    
    # Remove Rows with NA in 'proportion_change'
    merged_proportion_df <- merged_proportion_df %>% filter(!is.na(proportion_change))
    
    # Assign Groups to Patients
    merged_proportion_df$Groups <- lapply(merged_proportion_df$patient_id, function(pid) {
      groups <- c()
      if(pid %in% control_group) {
        groups <- c(groups, 'Control')
      }
      if(pid %in% experiment_group) {
        groups <- c(groups, 'Experiment')
      }
      if(pid %in% short_term_survivor_group) {
        groups <- c(groups, 'Short-term')
      }
      if(pid %in% long_term_survivor_group) {
        groups <- c(groups, 'Long-term')
      }
      return(groups)
    })
    
    # Expand Data Frame to Have One Row per Group per Patient
    expanded_df <- merged_proportion_df %>% unnest(Groups)
    
    # Define Custom Colors
    custom_colors <- c("Control" = "blue", 
                       "Experiment" = "red",
                       "Short-term" = "orange", 
                       "Long-term" = "green")
    
    # Ensure 'Groups' is a Factor with Desired Order
    expanded_df$Groups <- factor(expanded_df$Groups, levels = c("Control", "Experiment", "Short-term", "Long-term"))
    
    # Verify Group Assignments
    print(table(expanded_df$Groups))
    
    # ============================================================================================================
    # Create Line Plot
    
    # For Line Plot, use Control, Short-term, Long-term groups
    line_plot_df <- merged_proportion_df
    line_plot_df$Group <- ifelse(line_plot_df$patient_id %in% control_group, 'Control',
                                 ifelse(line_plot_df$patient_id %in% short_term_survivor_group, 'Short-term',
                                        ifelse(line_plot_df$patient_id %in% long_term_survivor_group, 'Long-term', NA)))
    line_plot_df <- line_plot_df %>% filter(!is.na(Group))
    line_plot_df$Group <- factor(line_plot_df$Group, levels = c("Control", "Short-term", "Long-term"))
    
    # Convert Data to Long Format for ggplot
    data_long <- line_plot_df %>%
      pivot_longer(cols = c(paste("Cluster_Proportion", timepoint_a, sep = "_"),
                            paste("Cluster_Proportion", timepoint_b, sep = "_")),
                   names_to = "Time_Point", 
                   values_to = "Cluster_Proportion")
    
    # Rename Time Points
    data_long$Time_Point <- factor(gsub("Cluster_Proportion_", "", data_long$Time_Point), 
                                   levels = c(timepoint_a, timepoint_b))
    
    # Generate Line Plot
    p_line <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Group)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(title = "Cluster Proportion Change by Group",
           x = "Time Point",
           y = "Cluster Proportion") +
      theme_minimal() +
      theme(axis.line = element_line(color = "black"),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line())
    
    # Add Line Plot to Plot List
    plot_list_line[[paste0(timepoint_a, "_vs_", timepoint_b, "_lineplot")]] <- p_line
    
    # ============================================================================================================
    # Create Boxplot with Expanded Data Frame
    
    # Generate Boxplot
    p_box <- ggplot(expanded_df, aes(x = Groups, y = proportion_change, fill = Groups)) +
      geom_boxplot(width = 0.5) +
      geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
           x = NULL,  # Remove x-axis label
           y = "Proportion Change") +
      theme_minimal() +
      theme(axis.line = element_line(color = "black"),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_blank(),  # Remove x-axis text
            axis.ticks.x = element_blank(),  # Remove x-axis ticks
            legend.position = "none",        # Remove legend
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line()) +
      # Calculate Maximum y-value for Annotation Placement
      ylim(NA, max(expanded_df$proportion_change, na.rm = TRUE) * 1.3)  # Increased multiplier for more space
    
    # Define Comparisons: Control vs Experiment, Control vs Short-term, Control vs Long-term
    comparisons <- list(
      c("Control", "Experiment"),
      c("Control", "Short-term"),
      c("Control", "Long-term")
    )
    
    # Initialize a List to Store p-values
    p_values <- list()
    
    # Perform Wilcoxon Tests for Each Comparison
    for (comp in comparisons){
      group1 <- comp[1]
      group2 <- comp[2]
      
      data1 <- expanded_df$proportion_change[expanded_df$Groups == group1]
      data2 <- expanded_df$proportion_change[expanded_df$Groups == group2]
      
      # Remove duplicates if any
      data1 <- unique(data.frame(patient_id = expanded_df$patient_id[expanded_df$Groups == group1], proportion_change = data1))
      data2 <- unique(data.frame(patient_id = expanded_df$patient_id[expanded_df$Groups == group2], proportion_change = data2))
      
      # Ensure no overlapping patients between groups
      overlapping_patients <- intersect(data1$patient_id, data2$patient_id)
      if(length(overlapping_patients) > 0){
        data2 <- data2[!data2$patient_id %in% overlapping_patients, ]
      }
      
      # Perform Wilcoxon Test if Both Groups Have At Least 2 Observations
      if (nrow(data1) >= 2 & nrow(data2) >= 2){
        test_result <- wilcox.test(data1$proportion_change, data2$proportion_change)
        p_value <- test_result$p.value
      } else {
        p_value <- NA_real_
      }
      
      # Store p-value
      p_values[[paste(group1, "vs", group2, sep = "_")]] <- p_value
      
      # Append to Statistical Summary Table
      stat_summary <- rbind(stat_summary, data.frame(
        celltype = celltype,
        timepoint_a = timepoint_a,
        timepoint_b = timepoint_b,
        comparison = paste(group1, "vs", group2, sep = "_"),
        p_value = p_value
      ))
    }
    
    # Define y-position increments to avoid overlapping annotations
    stat_max <- max(expanded_df$proportion_change, na.rm = TRUE)
    n_comparisons <- length(comparisons)
    y_positions <- seq(stat_max * 1.05, stat_max * 1.3, length.out = n_comparisons)
    
    # Add p-value Annotations
    for (i in 1:length(comparisons)){
      comp <- comparisons[[i]]
      group1 <- comp[1]
      group2 <- comp[2]
      p_value <- p_values[[paste(group1, "vs", group2, sep = "_")]]
      
      # Determine Significance Level and Label
      if (!is.na(p_value)){
        if (p_value < 0.001){
          p_label <- "***"
          p_color <- "red"
        } else if (p_value < 0.01){
          p_label <- "**"
          p_color <- "red"
        } else if (p_value < 0.05){
          p_label <- "*"
          p_color <- "red"
        } else{
          p_label <- "ns"
          p_color <- "black"
        }
      } else {
        p_label <- ""
        p_color <- "black"
      }
      
      # Calculate Midpoint for Comparison Lines
      x_start <- which(levels(expanded_df$Groups) == group1)
      x_end <- which(levels(expanded_df$Groups) == group2)
      
      # Check if both groups exist in the current data
      if(length(x_start) == 0 | length(x_end) == 0){
        warning(paste("One of the groups for comparison", paste(comp, collapse = " vs "), "is missing in the data."))
        next
      }
      
      # Add Comparison Line
      p_box <- p_box +
        geom_segment(aes(x = x_start, xend = x_end, y = y_positions[i], yend = y_positions[i]), inherit.aes = FALSE) +
        geom_text(aes(x = (x_start + x_end)/2, y = y_positions[i], label = p_label), 
                  inherit.aes = FALSE, color = p_color, size = 5)
    }
    
    # Add Boxplot to Plot List
    plot_list_box[[paste0(timepoint_a, "_vs_", timepoint_b, "_boxplot")]] <- p_box
  }
  
  # ============================================================================================================
  # Combine All Plots into a Grid
  
  # Combine Line Plots and Boxplots
  final_plot <- plot_grid(plotlist = c(plot_list_line, plot_list_box), 
                          ncol = length(pairs_list), 
                          nrow = 2, 
                          labels = "AUTO", 
                          align = 'v')
  
  # Define Plot Directory Based on Method
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  
  # Create Plot Directory if It Doesn't Exist
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Define File Name
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_All_Cells_Comparison.pdf")
  
  # Save the Grid Plot as a PDF
  ggsave(filename = file_name, 
         plot = final_plot, 
         device = "pdf", 
         width = 6 * length(pairs_list), 
         height = 6 * 2, 
         limitsize = FALSE)
}

# ============================================================================================================
# Save Statistical Summary as CSV

# Define Output CSV Path
output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material/ITT/Cluster_Proportion_Comparison_Results/ratio/statistical_summary.csv"

# Format p-values with Significance Asterisks
stat_summary <- stat_summary %>%
  mutate(
    significance = case_when(
      !is.na(p_value) & p_value < 0.001 ~ "***",
      !is.na(p_value) & p_value < 0.01  ~ "**",
      !is.na(p_value) & p_value < 0.05  ~ "*",
      !is.na(p_value)                    ~ "ns",
      TRUE                                ~ ""
    )
  )

# Save the Summary Table as CSV
write.csv(stat_summary, file = output_csv_path, row.names = FALSE)

# ============================================================================================================
# End of Script













# 
# ##################################################################################################################################################################################################################
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(Seurat)
# library(cowplot)
# 
# # ==================================================================================================================================================================
# # cluster proportion change comparison
# 
# get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
#   # Filter metadata for the cluster of interest and the specified timepoint
#   cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest & seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
#   timepoint_data <- seurat_metadata[seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
#   # Group the filtered data by patient
#   cluster_timpoint_grouped_data <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
#   timpoint_grouped_data <- split(timepoint_data, timepoint_data[[patient_col]])
#   
#   # Initialize lists to store patient_id and Cluster_Proportion_X
#   subject_ids <- c()
#   cluster_proportions <- c()
#   
#   # Loop through each patient
#   for (key in names(cluster_timpoint_grouped_data)) {
#     patient_cluster_timepoint_data <- cluster_timpoint_grouped_data[[key]]
#     patient_timepoint_data <- timpoint_grouped_data[[key]]
#     
#     # Calculate the total number of cells for the patient
#     patient_cluster_timepoint_cells <- nrow(patient_cluster_timepoint_data)
#     patient_timepoint_cells <- nrow(patient_timepoint_data)
#     
#     # Add the total number of cells to the results if there are any cells
#     if (patient_cluster_timepoint_cells > 0) {
#       # Calculate the cluster proportion
#       cluster_proportion <- patient_cluster_timepoint_cells / patient_timepoint_cells
#       
#       # Append the patient_id and Cluster_Proportion_X to the lists
#       subject_ids <- c(subject_ids, key)
#       cluster_proportions <- c(cluster_proportions, cluster_proportion)
#     }
#   }
#   
#   custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
#   # Create a data frame with patient_id and custom Cluster_Proportion_X column
#   result_df <- data.frame(patient_id = subject_ids)
#   result_df[[custom_column_name]] <- cluster_proportions
#   
#   return(result_df)
# }
# 
# # Helper function to determine significance based on p-value
# get_significance <- function(p) {
#   if (is.na(p)) {
#     return("ns")
#   } else if (p < 0.001) {
#     return("***")
#   } else if (p < 0.01) {
#     return("**")
#   } else if (p < 0.05) {
#     return("*")
#   } else {
#     return("ns")
#   }
# }
# 
# # Initialize an empty data frame to store comparison results
# comparison_results <- data.frame(
#   celltype = character(),
#   timepoint_A = character(),
#   timepoint_B = character(),
#   comparison = character(),
#   p_value = numeric(),
#   significance = character(),
#   stringsAsFactors = FALSE
# )
# 
# survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
# survival_df <- read.csv(survival_data)
# survival_df <- survival_df[survival_df$site == "UF",]
# 
# patient_col <- "Patient"
# timepoint_col <- "TimePoint"
# cluster_col <- "seurat_clusters"
# 
# method = "ratio"
# 
# timepoints <- c("Pre", "C1", "C2")
# pair_combinations <- combn(timepoints, 2)
# pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])
# 
# seurat_object <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
# seurat_metadata <- seurat_object@meta.data
# 
# # Define the mapping
# mapping <- list(
#   "Classical_Monocytes" = c(1, 3, 16, 12),
#   "Non_Classical_Monocytes" = c(9,33),
#   "cDC" = c(31),
#   "pDC" = c(36),
#   "NK" = c(0)
# )
# 
# # Convert the list to a data frame
# celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
#   data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
# }))
# 
# for (celltype in unique(celltype_to_cluster$celltype)){
#   cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
#   seurat_metadata <- seurat_metadata
#   
#   # Initialize a list to store the plots
#   plot_list_1 <- list()
#   plot_list_2 <- list()
#   
#   for (timepoint_pair in pairs_list){
#     timepoint_a <- timepoint_pair[1]
#     timepoint_b <- timepoint_pair[2]
#     
#     timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
#     timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
#     merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
#     merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], by = "patient_id")
#     merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
#     # Define custom colors
#     custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
#     
#     if (method == "difference") {
#       merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] - merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
#     } else if (method == "ratio") {
#       merged_proportion_df$proportion_change = merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] / merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
#     }
#     
#     if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
#         sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
#         sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
#       # Perform Wilcoxon test
#       test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
#       p_value = test_result$p.value
#     } else {
#       p_value = NA_real_  # Use NA_real_ for numerical context missing values
#     }
#     
#     # Determine significance based on p-value
#     significance <- get_significance(p_value)
#     
#     # Append the results to the comparison_results data frame
#     comparison_results <- rbind(comparison_results, data.frame(
#       celltype = celltype,
#       timepoint_A = timepoint_a,
#       timepoint_B = timepoint_b,
#       comparison = "control vs experiment",
#       p_value = p_value,
#       significance = significance,
#       stringsAsFactors = FALSE
#     ))
#     
#     # Remove the 'pre_' prefix from column names
#     colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
#     
#     # Convert data to long format for ggplot
#     data_long <- merged_proportion_df %>%
#       pivot_longer(cols = c(timepoint_a, timepoint_b), names_to = "Time_Point", values_to = "Cluster_Proportion")
#     
#     data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
#     
#     # Create the plot
#     p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
#       geom_line(size = 1) +
#       geom_point(size = 3) +
#       scale_color_manual(values = custom_colors) +
#       labs(title = "Cluster Proportion Change by Arm",
#            x = "Time Point",
#            y = "Cluster Proportion") +
#       theme_minimal() +
#       theme(axis.line = element_line(color = "black"),
#             axis.title.x = element_text(size = 20),  # Adjust x axis label size
#             axis.title.y = element_text(size = 20),  # Adjust y axis label size
#             plot.title = element_text(size = 20),
#             axis.text.y = element_text(size = 15),
#             axis.text.x = element_text(size = 12),  # Add axis lines
#             panel.grid.major = element_blank(),  # Remove major grid lines
#             panel.grid.minor = element_blank(),  # Remove minor grid lines
#             axis.ticks.y = element_line())
#     
#     # Add the plot to the list
#     plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p
#     
#     # Create boxplot
#     p_box <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
#       geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
#       # geom_jitter(position = position_dodge(width = 1.2)) +
#       geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), size = 4, alpha = 0.8) +
#       scale_fill_manual(values = custom_colors) +
#       labs(title = paste("Change between", timepoint_a, "and", timepoint_b),
#            x = NULL,  # Remove x-axis label
#            y = "Proportion Change") +
#       theme_minimal() +
#       theme(axis.line = element_line(color = "black"),
#             axis.title.x = element_blank(),  # Remove x axis label
#             axis.title.y = element_text(size = 20),  # Adjust y axis label size
#             plot.title = element_text(size = 20),
#             axis.text.y = element_text(size = 15),
#             axis.text.x = element_blank(),  # Remove x-axis text (tick labels)
#             panel.grid.major = element_blank(),  # Remove major grid lines
#             panel.grid.minor = element_blank(),  # Remove minor grid lines
#             axis.ticks.y = element_line(),
#             legend.position = "none") +  # Remove legend
#       annotate("text", x = 1.5, y = max(merged_proportion_df$proportion_change, na.rm = TRUE), 
#                label = paste("p-value =", ifelse(is.na(p_value), "NA", round(p_value, 4))), size = 4, vjust = -0.5)
#     
#     # Add the boxplot to the list
#     plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
#   }
#   
#   # Combine all patient plots into a grid
#   plot_list <- c(plot_list_1, plot_list_2)
#   final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), nrow = 2, align = 'v')
#   
#   # Define plot directory based on method
#   if (method == "difference") {
#     plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
#   } else {
#     plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
#   }
#   
#   # Create the plot directory if it doesn't exist
#   if (!dir.exists(plot_dir)) {
#     dir.create(plot_dir, recursive = TRUE)
#   }
#   
#   # Define the file name and save the plot
#   file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_All_Cells_Comparison.pdf")
#   ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
# }
# 
# # After processing all celltypes and timepoint pairs, save the comparison results to CSV
# output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/comparison_results.csv"
# write.csv(comparison_results, file = output_csv_path, row.names = FALSE)
# 
# # Optional: Print a message upon successful completion
# cat("Comparison results have been saved to:", output_csv_path, "\n")









#########################################################################################################################
# adding logFC column

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# =============================================================================================================
# cluster proportion change comparison
get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for the cluster of interest and the specified timepoint
  cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest & 
                                              seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
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

# Helper function to determine significance based on p-value
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

# Initialize an empty data frame to store comparison results
comparison_results <- data.frame(
  celltype = character(),
  timepoint_A = character(),
  timepoint_B = character(),
  comparison = character(),
  p_value = numeric(),
  significance = character(),
  fold_change = numeric(),    # <-- NEW COLUMN FOR FOLD CHANGE
  stringsAsFactors = FALSE
)

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

method = "difference"

timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- lapply(1:ncol(pair_combinations), function(i) pair_combinations[, i])

seurat_object <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
seurat_metadata <- seurat_object@meta.data

# Define the mapping
mapping <- list(
  "Classical_Monocytes" = c(1, 3, 16, 12),
  "Non_Classical_Monocytes" = c(9,33),
  "cDC" = c(31),
  "pDC" = c(36),
  "NK" = c(0)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

for (celltype in unique(celltype_to_cluster$celltype)){
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  seurat_metadata <- seurat_metadata
  
  # Initialize lists to store the plots
  plot_list_1 <- list()
  plot_list_2 <- list()
  
  for (timepoint_pair in pairs_list){
    timepoint_a <- timepoint_pair[1]
    timepoint_b <- timepoint_pair[2]
    
    # Get proportions at timepoint A and B
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, 
                                                                patient_col, 
                                                                timepoint_col, 
                                                                cluster_col, 
                                                                cluster_list, 
                                                                timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, 
                                                                patient_col, 
                                                                timepoint_col, 
                                                                cluster_col, 
                                                                cluster_list, 
                                                                timepoint_b)
    
    merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, 
                                  timepoint_b_cluster_proportion_df, 
                                  by = "patient_id")
    merged_proportion_df <- merge(merged_proportion_df, survival_df[, c("patient_id", "Arm")], 
                                  by = "patient_id")
    merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
    
    # Define custom colors
    custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
    
    # ---------------------------------------------------------------------------------
    # KEEP the original proportion_change calculation for plotting:
    if (method == "difference") {
      merged_proportion_df$proportion_change <- 
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] -
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    } else if (method == "ratio") {
      merged_proportion_df$proportion_change <- 
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")] /
        merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    }
    
    # ---------------------------------------------------------------------------------
    # NEW STEP: Calculate a "fold_change" statistic across all patients for CSV output.
    # We'll define log2 fold change as an overall measure to represent
    # both direction and magnitude. We add a small epsilon to avoid NaN/Inf.
    proportion_a_vals <- merged_proportion_df[, paste("Cluster_Proportion", timepoint_a, sep = "_")]
    proportion_b_vals <- merged_proportion_df[, paste("Cluster_Proportion", timepoint_b, sep = "_")]
    
    eps <- 1e-9
    log2_ratios <- log2((proportion_b_vals + eps) / (proportion_a_vals + eps))
    median_fold_change <- median(log2_ratios, na.rm = TRUE)
    if (is.na(median_fold_change)) {
      median_fold_change <- 0  # fallback in case everything is NA
    }
    
    # ---------------------------------------------------------------------------------
    # Perform statistical test if both arms have >= 2 samples
    if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
        sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
        sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
      test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
      p_value <- test_result$p.value
    } else {
      p_value <- NA_real_
    }
    
    # Determine significance based on p-value
    significance <- get_significance(p_value)
    
    # Append the results to the comparison_results data frame, now with fold_change
    comparison_results <- rbind(
      comparison_results, 
      data.frame(
        celltype = celltype,
        timepoint_A = timepoint_a,
        timepoint_B = timepoint_b,
        comparison = "control vs experiment",
        p_value = p_value,
        significance = significance,
        fold_change = median_fold_change,  # log2 fold change
        stringsAsFactors = FALSE
      )
    )
    
    # ---------------------------------------------------------------------------------
    # Remove the 'Cluster_Proportion_' prefix from column names for plotting
    colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
    
    # Convert data to long format for ggplot
    data_long <- merged_proportion_df %>%
      pivot_longer(cols = c(timepoint_a, timepoint_b), 
                   names_to = "Time_Point", 
                   values_to = "Cluster_Proportion")
    
    data_long$Time_Point <- factor(data_long$Time_Point, 
                                   levels = c(timepoint_a, timepoint_b))
    
    # -----------------------
    # PLOT 1: Dot/line plot
    p <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, 
                               group = patient_id, color = Arm)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
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
    
    # -----------------------
    # PLOT 2: Boxplot of proportion_change
    p_box <- ggplot(merged_proportion_df, 
                    aes(x = Arm, y = proportion_change, fill = Arm)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +
      geom_jitter(position = position_jitterdodge(jitter.width = 1, 
                                                  dodge.width = 1.2), 
                  size = 4, alpha = 0.8) +
      scale_fill_manual(values = custom_colors) +
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
               size = 4, vjust = -0.5)
    
    plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
  }
  
  # Combine all patient plots into a grid
  plot_list <- c(plot_list_1, plot_list_2)
  final_plot <- plot_grid(plotlist = plot_list, ncol = length(pairs_list), 
                          nrow = 2, align = 'v')
  
  # Define plot directory based on method
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  
  # Create the plot directory if it doesn't exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Define the file name and save the plot
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_All_Cells_Comparison.pdf")
  ggsave(filename = file_name, plot = final_plot, device = "pdf", 
         width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}

# After processing all celltypes and timepoint pairs, save the comparison results to CSV.
# Note: The file path below is set to the "difference" folder, adjust as needed.
output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/comparison_results.csv"
write.csv(comparison_results, file = output_csv_path, row.names = FALSE)

cat("Comparison results have been saved to:", output_csv_path, "\n")











# ######################################################################################################################################################################################################################
# # new annotation
# # Load Required Libraries
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(Seurat)
# library(cowplot)
# 
# # ==================================================================================================================================================================
# # Cluster Proportion Change Comparison
# 
# # Define the helper function to calculate cluster proportions
# get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
#   # Filter metadata for the clusters of interest and the specified timepoint
#   cluster_timepoint_data <- seurat_metadata %>%
#     filter((!!sym(cluster_col)) %in% clusters_of_interest & (!!sym(timepoint_col)) == timepoint_of_interest)
#   
#   timepoint_data <- seurat_metadata %>%
#     filter((!!sym(timepoint_col)) == timepoint_of_interest)
#   
#   # Group the filtered data by patient
#   cluster_timepoint_grouped <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
#   timepoint_grouped <- split(timepoint_data, timepoint_data[[patient_col]])
#   
#   # Initialize vectors to store results
#   subject_ids <- c()
#   cluster_proportions <- c()
#   
#   # Loop through each patient to calculate proportions
#   for (patient_id in names(cluster_timepoint_grouped)) {
#     patient_cluster_data <- cluster_timepoint_grouped[[patient_id]]
#     patient_timepoint_data <- timepoint_grouped[[patient_id]]
#     
#     patient_cluster_cells <- nrow(patient_cluster_data)
#     patient_total_cells <- nrow(patient_timepoint_data)
#     
#     if (patient_cluster_cells > 0) {
#       cluster_proportion <- patient_cluster_cells / patient_total_cells
#       subject_ids <- c(subject_ids, patient_id)
#       cluster_proportions <- c(cluster_proportions, cluster_proportion)
#     }
#   }
#   
#   # Create a data frame with the results
#   custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
#   result_df <- data.frame(patient_id = subject_ids, stringsAsFactors = FALSE)
#   result_df[[custom_column_name]] <- cluster_proportions
#   
#   return(result_df)
# }
# 
# # Define the helper function to determine significance
# get_significance <- function(p) {
#   if (is.na(p)) {
#     return("ns")
#   } else if (p < 0.001) {
#     return("***")
#   } else if (p < 0.01) {
#     return("**")
#   } else if (p < 0.05) {
#     return("*")
#   } else {
#     return("ns")
#   }
# }
# 
# # Initialize an empty data frame to store comparison results
# comparison_results <- data.frame(
#   celltype = character(),
#   timepoint_A = character(),
#   timepoint_B = character(),
#   comparison = character(),
#   p_value = numeric(),
#   significance = character(),
#   stringsAsFactors = FALSE
# )
# 
# # ==================================================================================================================================================================
# # Load Survival Data
# survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
# survival_df <- read.csv(survival_data)
# survival_df <- survival_df %>% filter(site == "UF")
# 
# # Define Columns
# patient_col <- "Patient"
# timepoint_col <- "TimePoint"
# cluster_col <- "seurat_clusters"
# 
# # Define Method
# method <- "difference"  # Options: "difference" or "ratio"
# 
# # Define Timepoints and Pair Combinations
# timepoints <- c("Pre", "C1", "C2")
# pair_combinations <- combn(timepoints, 2)
# pairs_list <- split(pair_combinations, col(pair_combinations))
# 
# # ==================================================================================================================================================================
# # Load Seurat Object for T Cells
# seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
# seurat_metadata_t_cells <- seurat_object_t_cells@meta.data
# 
# # ==================================================================================================================================================================
# # Define the Mapping from Celltypes to Clusters
# mapping <- list(
#   "Activated_CD4" = c(0),
#   "Effector_CD8" = c(1),
#   "Effector_Memory_Precursor_CD8" = c(2),
#   "Exhausted_T" = c(3),
#   "Gamma_Delta_T" = c(4),
#   "Active_CD4" = c(5),
#   "Naive_CD4" = c(6, 9, 18),
#   "Memory_CD4" = c(7),
#   "Stem_Like_CD8" = c(8),
#   "Effector_Memory_CD8" = c(10),
#   "Central_Memory_CD8" = c(12),
#   "GZMK_Effector_Memory_CD8" = c(13),
#   "Proliferating_Effector" = c(14, 16, 17)
# )
# 
# # Convert the Mapping to a Data Frame
# celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
#   data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]], stringsAsFactors = FALSE)
# }))
# 
# # ==================================================================================================================================================================
# # Iterate Over Each Celltype
# for (celltype in unique(celltype_to_cluster$celltype)) {
#   # Get the list of clusters for the current celltype
#   cluster_list <- celltype_to_cluster %>%
#     filter(celltype == !!celltype) %>%
#     pull(cluster)
#   
#   # Assign the T cell metadata
#   seurat_metadata <- seurat_metadata_t_cells
#   
#   # Initialize lists to store plots
#   plot_list_1 <- list()  # Line Plots (Dotplots)
#   plot_list_2 <- list()  # Boxplots
#   
#   # Iterate Over Each Pair of Timepoints
#   for (i in 1:ncol(pair_combinations)) {
#     timepoint_pair <- pair_combinations[, i]
#     timepoint_a <- timepoint_pair[1]
#     timepoint_b <- timepoint_pair[2]
#     
#     # Calculate Cluster Proportions for Both Timepoints
#     timepoint_a_cluster_proportion_df <- get_cluster_proportion(
#       seurat_metadata, patient_col, timepoint_col, cluster_col,
#       clusters_of_interest = cluster_list, timepoint_of_interest = timepoint_a
#     )
#     
#     timepoint_b_cluster_proportion_df <- get_cluster_proportion(
#       seurat_metadata, patient_col, timepoint_col, cluster_col,
#       clusters_of_interest = cluster_list, timepoint_of_interest = timepoint_b
#     )
#     
#     # Merge Proportion Data Frames by Patient ID
#     merged_proportion_df <- merge(timepoint_a_cluster_proportion_df, timepoint_b_cluster_proportion_df, by = "patient_id")
#     
#     # Merge with Survival Data to Get Arm Information
#     merged_proportion_df <- merge(merged_proportion_df, survival_df %>% select(patient_id, Arm), by = "patient_id")
#     merged_proportion_df$Arm <- factor(merged_proportion_df$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
#     
#     # Define Custom Colors
#     custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
#     
#     # Calculate Proportion Change Based on Method
#     if (method == "difference") {
#       merged_proportion_df$proportion_change <- merged_proportion_df[[paste("Cluster_Proportion", timepoint_b, sep = "_")]] -
#         merged_proportion_df[[paste("Cluster_Proportion", timepoint_a, sep = "_")]]
#     } else if (method == "ratio") {
#       merged_proportion_df$proportion_change <- merged_proportion_df[[paste("Cluster_Proportion", timepoint_b, sep = "_")]] /
#         merged_proportion_df[[paste("Cluster_Proportion", timepoint_a, sep = "_")]]
#     }
#     
#     # Perform Wilcoxon Test if Conditions are Met
#     if (all(c("MK-3475 Alone", "MK-3475 + MLA") %in% unique(merged_proportion_df$Arm)) &&
#         sum(merged_proportion_df$Arm == "MK-3475 Alone") >= 2 &&
#         sum(merged_proportion_df$Arm == "MK-3475 + MLA") >= 2) {
#       test_result <- wilcox.test(proportion_change ~ Arm, data = merged_proportion_df)
#       p_value <- test_result$p.value
#     } else {
#       p_value <- NA_real_
#     }
#     
#     # Determine Significance
#     significance <- get_significance(p_value)
#     
#     # Append Results to comparison_results Data Frame
#     comparison_results <- rbind(comparison_results, data.frame(
#       celltype = celltype,
#       timepoint_A = timepoint_a,
#       timepoint_B = timepoint_b,
#       comparison = "control vs experiment",
#       p_value = p_value,
#       significance = significance,
#       stringsAsFactors = FALSE
#     ))
#     
#     # Remove 'Cluster_Proportion_' Prefix from Column Names for Plotting
#     colnames(merged_proportion_df) <- sub("^Cluster_Proportion_", "", colnames(merged_proportion_df))
#     
#     # Convert Data to Long Format for ggplot (Line Plots)
#     data_long <- merged_proportion_df %>%
#       pivot_longer(
#         cols = c(timepoint_a, timepoint_b),
#         names_to = "Time_Point",
#         values_to = "Cluster_Proportion"
#       )
#     
#     data_long$Time_Point <- factor(data_long$Time_Point, levels = c(timepoint_a, timepoint_b))
#     
#     # Create the Line Plot (Dotplot) without x-axis label
#     p_line <- ggplot(data_long, aes(x = Time_Point, y = Cluster_Proportion, group = patient_id, color = Arm)) +
#       geom_line(size = 1) +
#       geom_point(size = 3) +
#       scale_color_manual(values = custom_colors) +
#       labs(
#         title = "Cluster Proportion Change by Arm",
#         x = NULL,  # Remove x-axis label
#         y = "Cluster Proportion"
#       ) +
#       theme_minimal() +
#       theme(
#         axis.line = element_line(color = "black"),
#         axis.title.x = element_blank(),  # Remove x-axis title
#         axis.title.y = element_text(size = 20),
#         plot.title = element_text(size = 20),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.y = element_line()
#       )
#     
#     # Add the Line Plot to the List
#     plot_list_1[[paste0(timepoint_a, "_vs_", timepoint_b, "_dotplot")]] <- p_line
#     
#     # Create the Boxplot without x-axis label, x-axis text, and legend
#     p_box <- ggplot(merged_proportion_df, aes(x = Arm, y = proportion_change, fill = Arm)) +
#       geom_boxplot(width = 0.5, position = position_dodge(width = 1.2), outlier.shape = NA) +
#       # geom_jitter(position = position_dodge(width = 1.2)) +
#       geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), size = 4, alpha = 0.8) +
#       scale_fill_manual(values = custom_colors) +
#       labs(
#         title = paste("Change between", timepoint_a, "and", timepoint_b),
#         x = NULL,  # Remove x-axis label
#         y = "Proportion Change"
#       ) +
#       theme_minimal() +
#       theme(
#         axis.line = element_line(color = "black"),
#         axis.title.x = element_blank(),  # Remove x-axis title
#         axis.title.y = element_text(size = 20),
#         plot.title = element_text(size = 20),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_blank(),  # Remove x-axis text (tick labels)
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.y = element_line(),
#         legend.position = "none"  # Remove legend
#       ) +
#       annotate(
#         "text",
#         x = 1.5,
#         y = max(merged_proportion_df$proportion_change, na.rm = TRUE),
#         label = ifelse(is.na(p_value), "p-value = NA", paste("p-value =", round(p_value, 4))),
#         size = 4,
#         vjust = -0.5
#       )
#     
#     # Add the Boxplot to the List
#     plot_list_2[[paste(timepoint_a, timepoint_b, sep = "_vs_")]] <- p_box
#   }
#   
#   # Combine All Plots into a Grid
#   plot_list <- c(plot_list_1, plot_list_2)
#   final_plot <- plot_grid(
#     plotlist = plot_list,
#     ncol = length(pairs_list),
#     nrow = 2,
#     align = 'v'
#   )
#   
#   # Define Plot Directory Based on Method
#   if (method == "difference") {
#     plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
#   } else {
#     plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
#   }
#   
#   # Create the Plot Directory if it Doesn't Exist
#   if (!dir.exists(plot_dir)) {
#     dir.create(plot_dir, recursive = TRUE)
#   }
#   
#   # Define the File Name and Save the Plot
#   file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
#   ggsave(
#     filename = file_name,
#     plot = final_plot,
#     device = "pdf",
#     width = 6 * length(pairs_list),
#     height = 6 * 2,
#     limitsize = FALSE
#   )
# }
# 
# # ==================================================================================================================================================================
# # Save the Comparison Results to CSV
# output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/comparison_results_T_Cells.csv"
# write.csv(comparison_results, file = output_csv_path, row.names = FALSE)
# 
# # Optional: Print a Completion Message
# cat("Comparison results have been saved to:", output_csv_path, "\n")









######################################################################################################################################################################
# adding logFC column in the comaprison results

# Load Required Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

# =============================================================================================================
# Helper Function: Calculate Cluster Proportion
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
# Initialize a Data Frame to Store Comparison Results
# Added a "fold_change" column
comparison_results <- data.frame(
  celltype     = character(),
  timepoint_A  = character(),
  timepoint_B  = character(),
  comparison   = character(),
  p_value      = numeric(),
  significance = character(),
  fold_change  = numeric(),  # <-- NEW COLUMN
  stringsAsFactors = FALSE
)

# =============================================================================================================
# Load Survival Data
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df %>% filter(site == "UF")

# Define Columns
patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

# Define Method ("difference" or "ratio")
method <- "ratio"

# Define Timepoints and Pair Combinations
timepoints <- c("Pre", "C1", "C2")
pair_combinations <- combn(timepoints, 2)
pairs_list <- split(pair_combinations, col(pair_combinations))

# =============================================================================================================
# Load Seurat Object for T Cells
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

# =============================================================================================================
# Define the Mapping from Celltypes to Clusters
mapping <- list(
  "Activated_CD4"                = c(0),
  "Effector_CD8"                 = c(1),
  "Effector_Memory_Precursor_CD8"= c(2),
  "Exhausted_T"                  = c(3),
  "Gamma_Delta_T"                = c(4),
  "Active_CD4"                   = c(5),
  "Naive_CD4"                    = c(6, 9, 18),
  "Memory_CD4"                   = c(7),
  "Stem_Like_CD8"                = c(8),
  "Effector_Memory_CD8"          = c(10),
  "Central_Memory_CD8"           = c(12),
  "GZMK_Effector_Memory_CD8"     = c(13),
  "Proliferating_Effector"       = c(14, 16, 17)
)

# Convert the Mapping to a Data Frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), 
             cluster  = mapping[[celltype]], 
             stringsAsFactors = FALSE)
}))

# =============================================================================================================
# Iterate Over Each Celltype
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
  
  # Combine All Plots into a Grid
  plot_list <- c(plot_list_1, plot_list_2)
  final_plot <- plot_grid(
    plotlist = plot_list,
    ncol = length(pairs_list),
    nrow = 2,
    align = 'v'
  )
  
  # Define Plot Directory Based on Method
  if (method == "difference") {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/"
  }
  
  # Create the Plot Directory if it Doesn't Exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Define the File Name and Save the Plot
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

# =============================================================================================================
# Save the Comparison Results to CSV (Note: Currently points to the 'ratio' folder)
output_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/ratio/comparison_results_T_Cells.csv"
write.csv(comparison_results, file = output_csv_path, row.names = FALSE)

# Optional: Print a Completion Message
cat("Comparison results have been saved to:", output_csv_path, "\n")




#############################################################################################
# exact cluster proportion and not change
#############################################################################################
# =====================================================================
# 0.  Load libraries ---------------------------------------------------
# =====================================================================
library(dplyr)      # data wrangling
library(ggplot2)    # plotting
library(tidyr)      # (not strictly needed but handy) :contentReference[oaicite:5]{index=5}
library(Seurat)     # single‑cell object handling :contentReference[oaicite:6]{index=6}
library(cowplot)    # assembling multi‑panel figures
# =====================================================================
# 1.  Define survivor cohorts -----------------------------------------
# =====================================================================
control_group            <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

survival_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_path) |>
  mutate(patient_id = as.integer(patient_id)) |>
  filter(site == "UF") |>
  mutate(
    SurvivalGroup = case_when(                     # vectorised if‑else
      patient_id %in% short_term_survivor_group ~ "ShortTerm",
      patient_id %in% long_term_survivor_group  ~ "LongTerm",
      patient_id %in% control_group             ~ "Control",
      TRUE                                      ~ NA_character_
    )
  ) |>
  filter(!is.na(SurvivalGroup)) |>
  arrange(OS.months.)

# =====================================================================
# 2.  Load Seurat object ----------------------------------------------
# =====================================================================
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")  # :contentReference[oaicite:7]{index=7}
meta <- seurat_obj@meta.data

# =====================================================================
# 3.  Mapping: cell‑type → cluster IDs --------------------------------
# =====================================================================
mapping <- list(
  "Activated_CD4"                 = c(0),
  "Effector_CD8"                  = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Exhausted_T"                   = c(3),
  "Gamma_Delta_T"                 = c(4),
  "Active_CD4"                    = c(5),
  "Naive_CD4"                     = c(6, 9, 18),
  "Memory_CD4"                    = c(7),
  "Stem_Like_CD8"                 = c(8),
  "Effector_Memory_CD8"           = c(10),
  "Central_Memory_CD8"            = c(12),
  "GZMK_Effector_Memory_CD8"      = c(13),
  "Proliferating_Effector"        = c(14, 16, 17)
)

# =====================================================================
# 4.  Function to build one 3‑panel plot ------------------------------
# =====================================================================
plot_cluster_proportions <- function(seurat_meta,
                                     patient_col   = "Patient",
                                     timepoint_col = "TimePoint",
                                     cluster_col   = "seurat_clusters",
                                     clusters_of_interest,
                                     surv_df,
                                     cohort_col    = "SurvivalGroup",
                                     timepoints    = c("Pre", "C1", "C2"),
                                     cohort_levels = c("Control", "ShortTerm", "LongTerm")) {
  # -------------------------------------------------------------------
  # 1.  Calculate per‑patient proportions for the chosen cluster set
  # -------------------------------------------------------------------
  meta_sub <- seurat_meta |>
    filter((!!sym(timepoint_col)) %in% timepoints) |>
    mutate(InCluster = (!!sym(cluster_col)) %in% clusters_of_interest)
  
  prop_df <- meta_sub |>
    group_by(!!sym(patient_col), !!sym(timepoint_col)) |>
    summarise(Cluster_Proportion = mean(InCluster), .groups = "drop")
  
  # -------------------------------------------------------------------
  # 2.  Attach cohort label and set factor levels for ordering
  # -------------------------------------------------------------------
  prop_df <- prop_df |>
    left_join(surv_df |> select(patient_id, !!sym(cohort_col)),
              by = setNames("patient_id", patient_col)) |>
    filter(!is.na(!!sym(cohort_col))) |>
    mutate(
      !!timepoint_col := factor(.data[[timepoint_col]], levels = timepoints),
      !!cohort_col    := factor(.data[[cohort_col]],    levels = cohort_levels)
    )
  
  # -------------------------------------------------------------------
  # 3.  Make one plot per cohort
  # -------------------------------------------------------------------
  plotlist <- lapply(cohort_levels, function(cohort_name) {
    
    ggplot(prop_df |> filter(.data[[cohort_col]] == cohort_name),
           aes(x   = .data[[timepoint_col]],
               y   = Cluster_Proportion,
               fill = .data[[timepoint_col]])) +          # NO group aesthetic here
      geom_boxplot(outlier.shape = NA, width = 0.55) +    # 1 box per time‑point
      geom_line(aes(group = .data[[patient_col]]),        # connect each patient
                linewidth = 0.7, colour = "black", alpha = 0.6) +
      geom_point(aes(group = .data[[patient_col]]),       # patient dots
                 size = 2.2, alpha = 0.8) +
      labs(title = paste0(cohort_name, " Cohort"),
           x = "Time‑point",
           y = "Cluster proportion") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })
  
  # -------------------------------------------------------------------
  # 4.  Combine the three cohort panels side‑by‑side
  # -------------------------------------------------------------------
  cowplot::plot_grid(plotlist = plotlist, nrow = 1)
}


# =====================================================================
# 5.  Create output directory -----------------------------------------
# =====================================================================
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Cluster_Proportion_Comparison_Results/not_change/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =====================================================================
# 6.  Loop over ALL cell‑types and save PDFs ---------------------------
# =====================================================================
timepoints_to_plot <- c("Pre", "C1", "C2", "C4", "C6", "C9")  # change if needed

for (celltype in names(mapping)) {
  
  message("Processing: ", celltype)
  
  clusters_of_interest <- mapping[[celltype]]
  
  plot_obj <- plot_cluster_proportions(
    seurat_meta          = meta,
    clusters_of_interest = clusters_of_interest,
    surv_df              = survival_df,
    timepoints           = timepoints_to_plot
  )
  
  ggsave(filename = file.path(output_dir,
                              paste0(celltype, "_cluster_proportion_boxplots.pdf")),
         plot     = plot_obj,
         device   = "pdf",
         width    = 12,
         height   = 4,      # one row of three boxes
         limitsize = FALSE)
}

message("All plots saved to: ", output_dir)

