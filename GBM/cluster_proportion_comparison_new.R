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
survival_df <- survival_df[survival_df$site == "UF" & survival_df$IDH == "NEG",]

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
  "pDC" = c(36)
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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/ratio/"
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
survival_df <- survival_df[survival_df$site == "UF" & survival_df$IDH == "NEG",]

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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_Cluster_Proportion_in_Monocyte_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}





######################################################################################################################################################################################################################
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF" & survival_df$IDH == "NEG",]

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
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/difference/"
  } else {
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Cluster_Proportion_Comparison_Results/ratio/"
  }
  file_name <- paste0(plot_dir, celltype, "_T_Cells_Cluster_Proportion_in_T_Cells_Comparison.pdf")
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = final_plot, device = "pdf", width = 6 * length(pairs_list), height = 6 * 2, limitsize = FALSE)
}


