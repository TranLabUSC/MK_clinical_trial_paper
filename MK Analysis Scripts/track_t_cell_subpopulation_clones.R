# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)

# Define the mapping
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
# Convert the list to a data frame if needed
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# Define all timepoints
all_timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF" & Patient != 1)

# Function to create clone tracking UMAPs
createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, subpop_seurat_obj, clone_cell_barcode_file, sample_clone_file, all_timepoints, output_folder) {

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
  umap_coords <- Embeddings(full_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = full_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(full_seurat_obj@meta.data)
  )
  print(dim(umap_df))

  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1

  for (row_tp in all_timepoints) {
    message("Processing clones from timepoint: ", row_tp)
    # Step 2: Identify the cell barcodes in the clusters at timepoint row_tp.
    cells_in_subpop_timepoint <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters & TimePoint == row_tp)

    if (length(cells_in_subpop_timepoint) == 0) {
      # No cells in subpopulation at this timepoint
      message("no cells at timepoint ", row_tp)
      next
    }

    # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
    clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop_timepoint]
    clones_from_cells <- unique(clones_from_cells)
    clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]

    # Step 4: Identify clones with proportion > 0 in samples corresponding to timepoint under consideration in sample_clone_df
    # sample_ids corresponding to timepoint row_tp
    sample_ids <- unique(subpop_seurat_obj@meta.data$origin[subpop_seurat_obj@meta.data$TimePoint == row_tp])

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
      # No intersected clones at this timepoint
      message("no intersected clones at timepoint ", row_tp)
      next
    }

    # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
    all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
    # Remove NAs and duplicates
    all_cells_for_clones <- unique(all_cells_for_clones)
    all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]

    # These are the cells_to_plot
    cells_to_plot <- all_cells_for_clones

    # Now, for each col_tp in all_timepoints:
    for (col_tp in all_timepoints) {
      # Get cells_in_col_tp
      cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
      cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)

      # Create a copy of umap_df
      plot_df <- umap_df
      plot_df$highlight <- "gray"
      plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "red"
      print(dim(plot_df))

      # Create the plot
      p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = highlight)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = c("gray" = "gray", "red" = "red")) +
        theme_void() +
        theme(legend.position = "none",
              plot.title = element_blank(),
              axis.title = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

      p <- p + coord_flip()
      # Optionally, add titles
      # p <- p + ggtitle(paste0("Clones from ", row_tp, " in ", col_tp))

      # Store the plot
      plot_list[[plot_counter]] <- p
      plot_counter <- plot_counter + 1
    }
  }
  
  message("total number of UMAPs: ", plot_counter)

  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }

  # Arrange plots into grid
  n_rows <- length(all_timepoints)
  n_cols <- length(all_timepoints)

  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)

  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap.pdf"))

  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)

  message("Saved UMAP grid for subpopulation: ", subpop_name)
}


# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Tracking_Results/new"
# Sample clone file path
sample_clone_file <- paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv")

# Loop over each subpopulation and perform the analysis
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Subset the seurat_obj to include only the clusters of the subpopulation
  subpop_seurat_obj <- subset(full_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found for subpopulation: ", subpop_name)
    next
  }
  
  # Call the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, subpop_seurat_obj, clone_cell_barcode_file, sample_clone_file, all_timepoints, output_folder)
  }, error = function(e) {
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}











##############################################################################################################################################################################################################################################################
# collapsing all rows into a single row
# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)

# Define the mapping with unique subpopulation names
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

# Convert the list to a data frame if needed (optional)
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# Define all timepoints
all_timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF" & Patient != 1)

# Function to create clone tracking UMAPs with single row per subpopulation
createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, subpop_seurat_obj, clone_cell_barcode_file, sample_clone_file, all_timepoints, output_folder) {
  
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
  umap_coords <- Embeddings(full_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = full_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(full_seurat_obj@meta.data)
  )
  print(dim(umap_df))
  
  
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
  # sample_ids corresponding to subpopulation
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
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in all_timepoints:
  for (col_tp in all_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create a copy of umap_df
    plot_df <- umap_df
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "red"
    print(dim(plot_df))
    
    # Create the plot with cluster labels
    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = highlight), size = 0.5) +
      scale_color_manual(values = c("gray" = "gray", "red" = "red")) +
      # geom_text_repel(data = cluster_centroids, aes(x = UMAP_1, y = UMAP_2, label = cluster),
      #                 color = "black", size = 3, segment.size = 0.2, max.overlaps = Inf) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))
    
    p <- p + coord_flip()
    # Optionally, add titles
    # p <- p + ggtitle(paste0(subpop_name, " clones in ", col_tp))
    
    # Store the plot
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over row_tp
  n_cols <- length(all_timepoints)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name)
}

# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Tracking_Results/new"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"

# Loop over each subpopulation and perform the analysis
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Subset the seurat_obj to include only the clusters of the subpopulation
  subpop_seurat_obj <- subset(full_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found for subpopulation: ", subpop_name)
    next
  }
  
  # Call the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, subpop_seurat_obj, clone_cell_barcode_file, sample_clone_file, all_timepoints, output_folder)
  }, error = function(e) {
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}









##############################################################################################################################################################################################################################################################
# segragating on the basis of cohort
# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the cohorts
# experiment_group <- c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21)
short_term_survivor_group <- c(5, 7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 13, 19, 20, 21)
control_group <- c(1, 4, 8, 9)

# Convert patient IDs to character if necessary
# experiment_group <- as.character(experiment_group)
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# cohorts <- list(experiment = experiment_group, control = control_group)
cohorts <- list(control = control_group, short_term_survivor = short_term_survivor_group, long_term_survivor = long_term_survivor_group)

# Function to create clone tracking UMAPs without cluster labels
createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, cohort_name, patient_ids) {
  
  # Subset full_seurat_obj to include only patients in the cohort
  cohort_seurat_obj <- subset(full_seurat_obj, subset = Patient %in% patient_ids)
  
  # Check if there are any cells left
  if (ncol(cohort_seurat_obj) == 0) {
    message("No cells found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Define all timepoints present in the cohort data
  cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))
  
  # Proceed only if there are timepoints
  if (length(cohort_timepoints) == 0) {
    message("No timepoints found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Create subpop_seurat_obj
  subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
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
    cell_id = rownames(cohort_seurat_obj@meta.data)
  )
  
  # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    # No cells in subpopulation
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
  # sample_ids corresponding to subpopulation
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
    message("No intersected clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in cohort_timepoints:
  for (col_tp in cohort_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create data frames for gray and red cells
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "red"
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    red_cells <- plot_df[plot_df$highlight == "red", ]
    
    # Create the plot without cluster labels
    p <- ggplot() +
      geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
      geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
      coord_flip()
    
    # Optionally, add titles
    # p <- p + ggtitle(paste0(subpop_name, " clones in ", col_tp))
    
    # Store the plot
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
  n_cols <- length(cohort_timepoints)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_", cohort_name, ".pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
}

# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/new_grouping"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"

# Loop over each subpopulation and perform the analysis for each cohort
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  for (cohort_name in names(cohorts)) {
    patient_ids <- cohorts[[cohort_name]]
    
    output_folder <- file.path(base_output_dir, subpop_name)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Call the function
    tryCatch({
      createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, cohort_name, patient_ids)
    }, error = function(e) {
      message("Error in processing subpopulation: ", subpop_name, " for cohort: ", cohort_name)
      message("Error message: ", e$message)
    })
  }
}









#################################################################################################################################################
# plotting for each patient

# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the patient groups
experiment_group <- c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21)
control_group <- c(1, 4, 8, 9, 11)  # Excluding Patient 1 as per previous subset

# Convert patient IDs to character if necessary
experiment_group <- as.character(experiment_group)
control_group <- as.character(control_group)

patient_groups <- list(experiment = experiment_group, control = control_group)

createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, patient_id) {
  
  # Subset full_seurat_obj to include only cells from the patient
  patient_seurat_obj <- subset(full_seurat_obj, subset = Patient == patient_id)
  
  # Check if there are any cells left
  if (ncol(patient_seurat_obj) == 0) {
    message("No cells found for patient: ", patient_id, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Define all timepoints present for the patient
  patient_timepoints <- intersect(all_timepoints, unique(patient_seurat_obj$TimePoint))
  
  # Proceed only if there are timepoints
  if (length(patient_timepoints) == 0) {
    message("No timepoints found for patient: ", patient_id, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Create subpop_seurat_obj
  subpop_seurat_obj <- subset(patient_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for patient: ", patient_id)
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
  umap_coords <- Embeddings(patient_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = patient_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(patient_seurat_obj@meta.data)
  )
  
  # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    # No cells in subpopulation
    message("No cells found in subpopulation: ", subpop_name, " for patient: ", patient_id)
    return()
  }
  
  # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
  # sample_ids corresponding to subpopulation
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
    message("No intersected clones found for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in patient_timepoints:
  for (col_tp in patient_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create data frames for gray and red cells
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "red"
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    red_cells <- plot_df[plot_df$highlight == "red", ]
    
    # Create the plot with red cells on top
    p <- ggplot() +
      geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
      geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
      coord_flip()
    
    # Optionally, add titles
    # p <- p + ggtitle(paste0("Patient ", patient_id, " - ", subpop_name, " in ", col_tp))
    
    # Store the plot
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
  n_cols <- length(patient_timepoints)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0("Patient_", patient_id, "_", subpop_name, "_clone_tracking_umap.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name, " for patient: ", patient_id)
}


# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"

# Create the patient_wise directory
patient_wise_dir <- file.path(base_output_dir, "patient_wise")
if (!dir.exists(patient_wise_dir)) {
  dir.create(patient_wise_dir)
}

# Loop over each patient group and patient, perform the analysis
for (group_name in names(patient_groups)) {
  patient_ids <- patient_groups[[group_name]]
  
  # Create group directory (experiment or control)
  group_dir <- file.path(patient_wise_dir, group_name)
  if (!dir.exists(group_dir)) {
    dir.create(group_dir)
  }
  
  for (patient_id in patient_ids) {
    message("Processing patient: ", patient_id, " in group: ", group_name)
    
    # Loop over each subpopulation
    for (subpop_name in names(mapping)) {
      clusters <- mapping[[subpop_name]]
      
      # Create subpopulation directory within group directory
      output_folder <- file.path(group_dir, subpop_name)
      if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      # Call the function
      tryCatch({
        createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, patient_id)
      }, error = function(e) {
        message("Error in processing subpopulation: ", subpop_name, " for patient: ", patient_id)
        message("Error message: ", e$message)
      })
    }
  }
}









#################################################################################################################################################
# still plotting for each patient, but only the clones which have more than 2 cells

# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
mapping <- list(
  "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the patient groups
experiment_group <- c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21)
control_group <- c(1, 4, 8, 9, 11)  # Excluding Patient 1 as per previous subset

# Convert patient IDs to character if necessary
experiment_group <- as.character(experiment_group)
control_group <- as.character(control_group)

patient_groups <- list(experiment = experiment_group, control = control_group)

createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, patient_id) {
  
  # Subset full_seurat_obj to include only cells from the patient
  patient_seurat_obj <- subset(full_seurat_obj, subset = Patient == patient_id)
  
  # Check if there are any cells left
  if (ncol(patient_seurat_obj) == 0) {
    message("No cells found for patient: ", patient_id, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Define all timepoints present for the patient
  patient_timepoints <- intersect(all_timepoints, unique(patient_seurat_obj$TimePoint))
  
  # Proceed only if there are timepoints
  if (length(patient_timepoints) == 0) {
    message("No timepoints found for patient: ", patient_id, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Create subpop_seurat_obj
  subpop_seurat_obj <- subset(patient_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for patient: ", patient_id)
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
  umap_coords <- Embeddings(patient_seurat_obj, "umap")
  umap_df <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    TimePoint = patient_seurat_obj@meta.data$TimePoint,
    cell_id = rownames(patient_seurat_obj@meta.data)
  )
  
  # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
  cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
  if (length(cells_in_subpop) == 0) {
    # No cells in subpopulation
    message("No cells found in subpopulation: ", subpop_name, " for patient: ", patient_id)
    return()
  }
  
  # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
  # sample_ids corresponding to subpopulation
  sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
  
  # Get clones with proportion > 0 in these samples
  clones_from_samples <- c()
  for (sample_id in sample_ids) {
    if (sample_id %in% colnames(sample_clone_df)) {
      clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 1]
      clones_from_samples <- c(clones_from_samples, clones_in_sample)
    }
  }
  clones_from_samples <- unique(clones_from_samples)
  
  # Step 5: Take intersection of the two
  intersected_clones <- intersect(clones_from_cells, clones_from_samples)
  
  if (length(intersected_clones) == 0) {
    # No intersected clones
    message("No intersected clones found for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in patient_timepoints:
  for (col_tp in patient_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create data frames for gray and red cells
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "red"
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    red_cells <- plot_df[plot_df$highlight == "red", ]
    
    # Create the plot with red cells on top
    p <- ggplot() +
      geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
      geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
      coord_flip()
    
    # Optionally, add titles
    # p <- p + ggtitle(paste0("Patient ", patient_id, " - ", subpop_name, " in ", col_tp))
    
    # Store the plot
    plot_list[[plot_counter]] <- p
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name, " in patient: ", patient_id)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
  n_cols <- length(patient_timepoints)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0("Patient_", patient_id, "_", subpop_name, "_clone_tracking_umap.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name, " for patient: ", patient_id)
}


# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/multiple_clones"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_absolute_all_T_cells.csv"

# Create the patient_wise directory
patient_wise_dir <- file.path(base_output_dir, "patient_wise")
if (!dir.exists(patient_wise_dir)) {
  dir.create(patient_wise_dir)
}

# Loop over each patient group and patient, perform the analysis
for (group_name in names(patient_groups)) {
  patient_ids <- patient_groups[[group_name]]
  
  # Create group directory (experiment or control)
  group_dir <- file.path(patient_wise_dir, group_name)
  if (!dir.exists(group_dir)) {
    dir.create(group_dir)
  }
  
  for (patient_id in patient_ids) {
    message("Processing patient: ", patient_id, " in group: ", group_name)
    
    # Loop over each subpopulation
    for (subpop_name in names(mapping)) {
      clusters <- mapping[[subpop_name]]
      
      # Create subpopulation directory within group directory
      output_folder <- file.path(group_dir, subpop_name)
      if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      # Call the function
      tryCatch({
        createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, patient_id)
      }, error = function(e) {
        message("Error in processing subpopulation: ", subpop_name, " for patient: ", patient_id)
        message("Error message: ", e$message)
      })
    }
  }
}



##############################################################################################################################################################################################################################################################
# segragating on the basis of cohort
# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the cohorts
control_group <- c(1, 4, 8, 9)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character if necessary
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Combine experiment group
experiment_group <- c(short_term_survivor_group, long_term_survivor_group)

# Create a combined cohorts list
cohorts <- list(
  control = control_group,
  experiment = experiment_group
)

# Add SurvivorGroup metadata to the full_seurat_obj
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term", NA)
)

# # Function to create clone tracking UMAPs with different colors for survivor groups
# createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, cohort_name, patient_ids) {
#   
#   # Subset full_seurat_obj to include only patients in the cohort
#   cohort_seurat_obj <- subset(full_seurat_obj, subset = Patient %in% patient_ids)
#   
#   # Check if there are any cells left
#   if (ncol(cohort_seurat_obj) == 0) {
#     message("No cells found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Define all timepoints present in the cohort data
#   cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))
#   
#   # Proceed only if there are timepoints
#   if (length(cohort_timepoints) == 0) {
#     message("No timepoints found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Create subpop_seurat_obj
#   subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
#   
#   # Check if there are any cells in subpopulation
#   if (ncol(subpop_seurat_obj) == 0) {
#     message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Read in clone_cell_barcode_file
#   clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
#   clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
#   
#   # Ensure 'barcode' and 'cdr3' columns exist
#   if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
#     stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
#   }
#   
#   # Read in sample_clone_file
#   sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
#   colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))
#   
#   # Ensure sample IDs are in columns
#   if (ncol(sample_clone_df) == 0) {
#     stop("sample_clone_file must contain sample IDs as columns.")
#   }
#   
#   # Prepare umap_df
#   umap_coords <- Embeddings(cohort_seurat_obj, "umap")
#   umap_df <- data.frame(
#     UMAP_1 = umap_coords[,1],
#     UMAP_2 = umap_coords[,2],
#     TimePoint = cohort_seurat_obj@meta.data$TimePoint,
#     cell_id = rownames(cohort_seurat_obj@meta.data),
#     SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup
#   )
#   
#   # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
#   cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
#   if (length(cells_in_subpop) == 0) {
#     # No cells in subpopulation
#     message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
#     return()
#   }
#   
#   # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
#   clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
#   clones_from_cells <- unique(clones_from_cells)
#   clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
#   
#   if (length(clones_from_cells) == 0) {
#     # No clones found for the subpopulation
#     message("No clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
#     return()
#   }
#   
#   # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
#   # sample_ids corresponding to subpopulation
#   sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
#   
#   # Get clones with proportion > 0 in these samples
#   clones_from_samples <- c()
#   for (sample_id in sample_ids) {
#     if (sample_id %in% colnames(sample_clone_df)) {
#       clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
#       clones_from_samples <- c(clones_from_samples, clones_in_sample)
#     }
#   }
#   clones_from_samples <- unique(clones_from_samples)
#   
#   # Step 5: Take intersection of the two
#   intersected_clones <- intersect(clones_from_cells, clones_from_samples)
#   
#   if (length(intersected_clones) == 0) {
#     # No intersected clones
#     message("No intersected clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
#     return()
#   }
#   
#   # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
#   all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
#   # Remove NAs and duplicates
#   all_cells_for_clones <- unique(all_cells_for_clones)
#   all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
#   
#   # These are the cells_to_plot
#   cells_to_plot <- all_cells_for_clones
#   
#   # Initialize plot list
#   plot_list <- list()
#   plot_counter <- 1
#   
#   # Now, for each col_tp in cohort_timepoints:
#   for (col_tp in cohort_timepoints) {
#     # Get cells_in_col_tp
#     cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
#     cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
#     
#     # Create data frames for gray and highlighted cells
#     plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
#     plot_df$highlight <- "gray"
#     plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"
#     
#     gray_cells <- plot_df[plot_df$highlight == "gray", ]
#     
#     if (cohort_name == "experiment") {
#       # For highlighted cells, assign color based on SurvivorGroup
#       highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
#       highlighted_cells$color <- ifelse(highlighted_cells$SurvivorGroup == "short_term", "blue", "red")
#       
#       blue_cells <- highlighted_cells[highlighted_cells$color == "blue", ]
#       red_cells <- highlighted_cells[highlighted_cells$color == "red", ]
#       
#       # Create the plot
#       p <- ggplot() +
#         geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
#         geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
#         geom_point(data = blue_cells, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5, alpha = 1) +
#         theme_void() +
#         theme(legend.position = "none",
#               plot.title = element_blank(),
#               axis.title = element_blank(),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               panel.background = element_blank(),
#               axis.line = element_line(colour = "black"),
#               plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
#         coord_flip()
#     } else {
#       # For control cohort, proceed as before
#       highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
#       
#       p <- ggplot() +
#         geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
#         geom_point(data = highlighted_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
#         theme_void() +
#         theme(legend.position = "none",
#               plot.title = element_blank(),
#               axis.title = element_blank(),
#               panel.border = element_rect(colour = "black", fill=NA, size=1),
#               panel.background = element_blank(),
#               axis.line = element_line(colour = "black"),
#               plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
#         coord_flip()
#     }
#     
#     # Optionally, add titles
#     # p <- p + ggtitle(paste0(subpop_name, " clones in ", col_tp))
#     
#     # Store the plot
#     plot_list[[plot_counter]] <- p
#     plot_counter <- plot_counter + 1
#   }
#   
#   message("Total number of UMAPs: ", plot_counter - 1)
#   
#   # Check if we have any plots
#   if (length(plot_list) == 0) {
#     message("No plots to save for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
#     return()
#   }
#   
#   # Arrange plots into grid
#   n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
#   n_cols <- length(cohort_timepoints)
#   
#   plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
#   
#   # Save the grid to a file
#   output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_", cohort_name, ".pdf"))
#   
#   ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
#   
#   message("Saved UMAP grid for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
# }

createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, cohort_name, patient_ids) {
  
  # Subset full_seurat_obj to include only patients in the cohort
  cohort_seurat_obj <- subset(full_seurat_obj, subset = Patient %in% patient_ids)
  
  # Check if there are any cells left
  if (ncol(cohort_seurat_obj) == 0) {
    message("No cells found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Define all timepoints present in the cohort data
  cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))
  
  # Proceed only if there are timepoints
  if (length(cohort_timepoints) == 0) {
    message("No timepoints found for cohort: ", cohort_name, " in subpopulation: ", subpop_name)
    return()
  }
  
  # Create subpop_seurat_obj
  subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
  
  # Check if there are any cells in subpopulation
  if (ncol(subpop_seurat_obj) == 0) {
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
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
    message("No cells found in subpopulation: ", subpop_name, " for cohort: ", cohort_name)
    return()
  }
  
  # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
  clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
  clones_from_cells <- unique(clones_from_cells)
  clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
  
  if (length(clones_from_cells) == 0) {
    # No clones found for the subpopulation
    message("No clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
  # sample_ids corresponding to subpopulation
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
    message("No intersected clones found for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in cohort_timepoints:
  for (col_tp in cohort_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create data frames for gray and highlighted cells
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    
    if (cohort_name == "experiment") {
      # For highlighted cells, assign color based on SurvivorGroup
      highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
      highlighted_cells$color <- ifelse(highlighted_cells$SurvivorGroup == "short_term", "blue", "red")
      
      # Determine the number of cells in each group
      num_blue_cells <- sum(highlighted_cells$color == "blue")
      num_red_cells <- sum(highlighted_cells$color == "red")
      
      # Find the minimum number of cells between the two groups
      min_cells <- min(num_blue_cells, num_red_cells)
      
      # Downsample the larger group
      if (min_cells == 0) {
        # If one group has zero cells, skip downsampling
        blue_cells <- highlighted_cells[highlighted_cells$color == "blue", ]
        red_cells <- highlighted_cells[highlighted_cells$color == "red", ]
      } else if (num_blue_cells > num_red_cells) {
        blue_cells <- highlighted_cells[highlighted_cells$color == "blue", ]
        red_cells <- highlighted_cells[highlighted_cells$color == "red", ]
        
        # Randomly sample blue cells to match red cells
        set.seed(123)  # For reproducibility
        blue_cells <- blue_cells[sample(nrow(blue_cells), min_cells), ]
      } else if (num_red_cells > num_blue_cells) {
        blue_cells <- highlighted_cells[highlighted_cells$color == "blue", ]
        red_cells <- highlighted_cells[highlighted_cells$color == "red", ]
        
        # Randomly sample red cells to match blue cells
        set.seed(123)  # For reproducibility
        red_cells <- red_cells[sample(nrow(red_cells), min_cells), ]
      } else {
        # If equal, no downsampling needed
        blue_cells <- highlighted_cells[highlighted_cells$color == "blue", ]
        red_cells <- highlighted_cells[highlighted_cells$color == "red", ]
      }
      
      # Create the plot
      p <- ggplot() +
        geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
        geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
        geom_point(data = blue_cells, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5, alpha = 1) +
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
      
    } else {
      # For control cohort, proceed as before
      highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
      
      p <- ggplot() +
        geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5) +
        geom_point(data = highlighted_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1) +
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
    }
    
    # Store the plot
    plot_list[[plot_counter]] <- p + coord_flip()
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
  n_cols <- length(cohort_timepoints)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_", cohort_name, "_equal_numbers_fixed_axes.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name, " in cohort: ", cohort_name)
}

# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/temp_equal"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"

# Loop over each subpopulation and perform the analysis for each cohort
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  for (cohort_name in names(cohorts)) {
    patient_ids <- cohorts[[cohort_name]]
    
    output_folder <- file.path(base_output_dir, subpop_name)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Call the function
    tryCatch({
      createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder, cohort_name, patient_ids)
    }, error = function(e) {
      message("Error in processing subpopulation: ", subpop_name, " for cohort: ", cohort_name)
      message("Error message: ", e$message)
    })
  }
}




##############################################################################################################################################################################################################################################################
# all cohorts in a single UMAPs denoted by different colors (equal numbers)


# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the cohorts
control_group <- c(1, 4, 8, 9)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character if necessary
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Add SurvivorGroup metadata to the full_seurat_obj
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term",
         ifelse(full_seurat_obj$Patient %in% control_group, "control", NA))
)

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
  # sample_ids corresponding to subpopulation
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
  
  # Step 6: Identify all cell barcodes corresponding to the intersected set of clones from clone_cell_barcode_df
  all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
  # Remove NAs and duplicates
  all_cells_for_clones <- unique(all_cells_for_clones)
  all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
  
  # These are the cells_to_plot
  cells_to_plot <- all_cells_for_clones
  
  # Initialize plot list
  plot_list <- list()
  plot_counter <- 1
  
  # Now, for each col_tp in cohort_timepoints:
  for (col_tp in cohort_timepoints) {
    # Get cells_in_col_tp
    cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
    cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
    
    # Create data frames for gray and highlighted cells
    plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
    plot_df$highlight <- "gray"
    plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    
    # For highlighted cells, assign color based on SurvivorGroup
    highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
    highlighted_cells$color <- highlighted_cells$SurvivorGroup
    
    # Handle NAs (if any)
    highlighted_cells$color[is.na(highlighted_cells$color)] <- "unknown"
    
    # Separate cells by group
    red_cells <- highlighted_cells[highlighted_cells$color == "long_term", ]
    blue_cells <- highlighted_cells[highlighted_cells$color == "short_term", ]
    yellow_cells <- highlighted_cells[highlighted_cells$color == "control", ]
    
    # Determine the number of cells in each group
    num_red_cells <- nrow(red_cells)
    num_blue_cells <- nrow(blue_cells)
    num_yellow_cells <- nrow(yellow_cells)
    
    # Get counts of groups with cells
    cell_counts <- c(num_red_cells, num_blue_cells, num_yellow_cells)
    cell_counts_nonzero <- cell_counts[cell_counts > 0]
    
    # Proceed if at least one group has cells
    if (length(cell_counts_nonzero) > 0) {
      # Find the minimum number of cells among the groups with cells
      min_cells <- min(cell_counts_nonzero)
      
      set.seed(123)  # For reproducibility
      if (num_red_cells > min_cells) {
        red_cells <- red_cells[sample(nrow(red_cells), min_cells), ]
      }
      if (num_blue_cells > min_cells) {
        blue_cells <- blue_cells[sample(nrow(blue_cells), min_cells), ]
      }
      if (num_yellow_cells > min_cells) {
        yellow_cells <- yellow_cells[sample(nrow(yellow_cells), min_cells), ]
      }
    } else {
      # If no group has cells, set data frames to have zero rows but same columns
      red_cells <- red_cells[FALSE, ]
      blue_cells <- blue_cells[FALSE, ]
      yellow_cells <- yellow_cells[FALSE, ]
    }
    
    # Create the plot
    p <- ggplot()
    
    if (nrow(gray_cells) > 0) {
      p <- p + geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5)
    }
    if (nrow(red_cells) > 0) {
      p <- p + geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1)
    }
    if (nrow(blue_cells) > 0) {
      p <- p + geom_point(data = blue_cells, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5, alpha = 1)
    }
    if (nrow(yellow_cells) > 0) {
      p <- p + geom_point(data = yellow_cells, aes(x = UMAP_1, y = UMAP_2), color = "yellow", size = 0.5, alpha = 1)
    }
    
    # If no layers were added, skip this plot
    if (length(p$layers) == 0) {
      message("No data to plot for subpopulation: ", subpop_name, " at timepoint: ", col_tp)
      next
    }
    
    p <- p + coord_fixed() +
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
    
    # Store the plot
    plot_list[[plot_counter]] <- p + coord_flip()
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1  # Only one row since we are not looping over timepoints for clone identification
  n_cols <- length(plot_list)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_equal_numbers_fixed_axes.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name)
}


# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/temp_equal"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"


# Loop over each subpopulation and perform the analysis
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Call the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }, error = function(e) {
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}





##############################################################################################################################################################################################################################################################
# all cohorts in a single UMAPs denoted by different colors (not equal numbers)


# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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

# Read the Seurat object (replace with your actual file path)
full_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF")

# Define the cohorts
control_group <- c(1, 4, 8, 9)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character if necessary
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Add SurvivorGroup metadata to the full_seurat_obj
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term",
         ifelse(full_seurat_obj$Patient %in% control_group, "control", NA))
)

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
  
  # Initialize plot list
  plot_list <- list()
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
    
    gray_cells <- plot_df[plot_df$highlight == "gray", ]
    
    # Assign colors based on SurvivorGroup
    highlighted_cells <- plot_df[plot_df$highlight == "highlighted", ]
    highlighted_cells$color <- highlighted_cells$SurvivorGroup
    highlighted_cells$color[is.na(highlighted_cells$color)] <- "unknown"
    
    # Separate cells by group
    red_cells <- highlighted_cells[highlighted_cells$color == "long_term", ]
    blue_cells <- highlighted_cells[highlighted_cells$color == "short_term", ]
    yellow_cells <- highlighted_cells[highlighted_cells$color == "control", ]
    
    # Create the plot
    p <- ggplot()
    
    if (nrow(gray_cells) > 0) {
      p <- p + geom_point(data = gray_cells, aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5)
    }
    if (nrow(red_cells) > 0) {
      p <- p + geom_point(data = red_cells, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.5, alpha = 1)
    }
    if (nrow(blue_cells) > 0) {
      p <- p + geom_point(data = blue_cells, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5, alpha = 1)
    }
    if (nrow(yellow_cells) > 0) {
      p <- p + geom_point(data = yellow_cells, aes(x = UMAP_1, y = UMAP_2), color = "yellow", size = 0.5, alpha = 1)
    }
    
    # If no layers were added, skip this plot
    if (length(p$layers) == 0) {
      message("No data to plot for subpopulation: ", subpop_name, " at timepoint: ", col_tp)
      next
    }
    
    p <- p + coord_fixed() +
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
    
    # Store the plot
    plot_list[[plot_counter]] <- p + coord_flip()
    plot_counter <- plot_counter + 1
  }
  
  message("Total number of UMAPs: ", plot_counter - 1)
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }
  
  # Arrange plots into grid
  n_rows <- 1
  n_cols <- length(plot_list)
  
  plot_grid_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
  
  # Save the grid to a file
  output_file <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf"))
  
  ggsave(filename = output_file, plot = plot_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
  
  message("Saved UMAP grid for subpopulation: ", subpop_name)
}



# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/temp"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"


# Loop over each subpopulation and perform the analysis
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Call the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }, error = function(e) {
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}




##############################################################################################################################################################################################################################################################
# density plots
# all cohorts in a single UMAPs denoted by different colors (not equal numbers)


# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# Define the mapping with unique subpopulation names
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
full_seurat_obj <- subset(full_seurat_obj, subset = Site == "UF" & Patient %in% all_patients)
full_seurat_obj <- subset(full_seurat_obj, subset = seurat_clusters != 13)

# Define the cohorts
control_group <- c(1, 4, 8, 9)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group <- c(2, 3, 5, 13, 19, 20, 21)

# Convert patient IDs to character if necessary
control_group <- as.character(control_group)
short_term_survivor_group <- as.character(short_term_survivor_group)
long_term_survivor_group <- as.character(long_term_survivor_group)

# Add SurvivorGroup metadata to the full_seurat_obj
full_seurat_obj$SurvivorGroup <- ifelse(
  full_seurat_obj$Patient %in% short_term_survivor_group, "short_term",
  ifelse(full_seurat_obj$Patient %in% long_term_survivor_group, "long_term",
         ifelse(full_seurat_obj$Patient %in% control_group, "control", NA))
)

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

    # Create UMAP plot
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

    # Create Density Plot with Gray Background
    # Remove gray cells for density plot
    density_df <- plot_df[plot_df$color != "gray", ]

    # Start with gray background
    p_density <- ggplot() +
      geom_point(data = plot_df[plot_df$color == "gray", ], aes(x = UMAP_1, y = UMAP_2), color = "gray", size = 0.5, alpha = 0.5)

    # Add density contours if there are highlighted cells
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

    # Store the density plot
    density_plot_list[[plot_counter]] <- p_density + coord_flip()
    plot_counter <- plot_counter + 1
  }

  message("Total number of UMAPs and density plots: ", plot_counter - 1)

  # Check if we have any plots
  if (length(umap_plot_list) == 0 || length(density_plot_list) == 0) {
    message("No plots to save for subpopulation: ", subpop_name)
    return()
  }

  # Arrange UMAP plots into grid
  n_rows <- 1
  n_cols <- length(umap_plot_list)

  umap_grid_combined <- plot_grid(plotlist = umap_plot_list, ncol = n_cols)

  # Save the UMAP grid to a file
  output_file_umap <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf"))

  ggsave(filename = output_file_umap, plot = umap_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)

  message("Saved UMAP grid for subpopulation: ", subpop_name)

  # Arrange Density plots into grid
  density_grid_combined <- plot_grid(plotlist = density_plot_list, ncol = n_cols)

  # Save the Density grid to a file
  output_file_density <- file.path(output_folder, paste0(subpop_name, "_density_plots_with_background_all_patients_all_cells_fixed_axes.pdf"))

  ggsave(filename = output_file_density, plot = density_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)

  message("Saved density plot grid for subpopulation: ", subpop_name)
}

# createCloneTrackingUMAPs <- function(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder) {
#   
#   # Use full_seurat_obj to include all patients
#   cohort_seurat_obj <- full_seurat_obj
#   
#   # Define all timepoints present in the data
#   cohort_timepoints <- intersect(all_timepoints, unique(cohort_seurat_obj$TimePoint))
#   
#   # Create subpop_seurat_obj
#   subpop_seurat_obj <- subset(cohort_seurat_obj, seurat_clusters %in% clusters)
#   
#   # Check if there are any cells in subpopulation
#   if (ncol(subpop_seurat_obj) == 0) {
#     message("No cells found in subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Read in clone_cell_barcode_file
#   clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
#   clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
#   
#   # Ensure 'barcode' and 'cdr3' columns exist
#   if (!all(c('barcode', 'cdr3') %in% colnames(clone_cell_barcode_df))) {
#     stop("clone_cell_barcode_file must contain 'barcode' and 'cdr3' columns.")
#   }
#   
#   # Read in sample_clone_file
#   sample_clone_df <- read.csv(sample_clone_file, row.names = 1, check.names = FALSE)
#   colnames(sample_clone_df) <- sub("_sc$", "", colnames(sample_clone_df))
#   
#   # Ensure sample IDs are in columns
#   if (ncol(sample_clone_df) == 0) {
#     stop("sample_clone_file must contain sample IDs as columns.")
#   }
#   
#   # Prepare umap_df
#   umap_coords <- Embeddings(cohort_seurat_obj, "umap")
#   umap_df <- data.frame(
#     UMAP_1 = umap_coords[,1],
#     UMAP_2 = umap_coords[,2],
#     TimePoint = cohort_seurat_obj@meta.data$TimePoint,
#     cell_id = rownames(cohort_seurat_obj@meta.data),
#     SurvivorGroup = cohort_seurat_obj@meta.data$SurvivorGroup
#   )
#   
#   # Get overall x and y limits
#   x_limits <- range(umap_df$UMAP_1)
#   y_limits <- range(umap_df$UMAP_2)
#   
#   # Step 2: Identify the cell barcodes in the clusters (subpopulation) across all timepoints
#   cells_in_subpop <- WhichCells(subpop_seurat_obj, expression = seurat_clusters %in% clusters)
#   if (length(cells_in_subpop) == 0) {
#     # No cells in subpopulation
#     message("No cells found in subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Step 3: Get clones corresponding to these cell barcodes from clone_cell_barcode_df
#   clones_from_cells <- clone_cell_barcode_df$'cdr3'[clone_cell_barcode_df$'barcode' %in% cells_in_subpop]
#   clones_from_cells <- unique(clones_from_cells)
#   clones_from_cells <- clones_from_cells[!is.na(clones_from_cells)]
#   
#   if (length(clones_from_cells) == 0) {
#     # No clones found for the subpopulation
#     message("No clones found for subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Step 4: Identify clones with proportion > 0 in samples corresponding to the subpopulation
#   sample_ids <- unique(subpop_seurat_obj@meta.data$origin)
#   
#   # Get clones with proportion > 0 in these samples
#   clones_from_samples <- c()
#   for (sample_id in sample_ids) {
#     if (sample_id %in% colnames(sample_clone_df)) {
#       clones_in_sample <- rownames(sample_clone_df)[sample_clone_df[, sample_id] > 0]
#       clones_from_samples <- c(clones_from_samples, clones_in_sample)
#     }
#   }
#   clones_from_samples <- unique(clones_from_samples)
#   
#   # Step 5: Take intersection of the two
#   intersected_clones <- intersect(clones_from_cells, clones_from_samples)
#   
#   if (length(intersected_clones) == 0) {
#     # No intersected clones
#     message("No intersected clones found for subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Step 6: Identify all cell barcodes corresponding to the intersected clones
#   all_cells_for_clones <- clone_cell_barcode_df$'barcode'[clone_cell_barcode_df$'cdr3' %in% intersected_clones]
#   all_cells_for_clones <- unique(all_cells_for_clones)
#   all_cells_for_clones <- all_cells_for_clones[!is.na(all_cells_for_clones)]
#   
#   # Cells to plot
#   cells_to_plot <- all_cells_for_clones
#   
#   # Initialize plot lists
#   umap_plot_list <- list()
#   density_plot_list <- list()
#   plot_counter <- 1
#   
#   # Loop over timepoints
#   for (col_tp in cohort_timepoints) {
#     # Cells in current timepoint
#     cells_in_col_tp <- umap_df$cell_id[umap_df$TimePoint == col_tp]
#     cells_to_highlight <- intersect(cells_to_plot, cells_in_col_tp)
#     
#     # Prepare data frames
#     plot_df <- umap_df[umap_df$cell_id %in% cells_in_col_tp, ]
#     plot_df$highlight <- "gray"
#     plot_df$highlight[plot_df$cell_id %in% cells_to_highlight] <- "highlighted"
#     
#     # Assign colors based on SurvivorGroup
#     plot_df$color <- "gray"
#     plot_df$color[plot_df$highlight == "highlighted"] <- plot_df$SurvivorGroup[plot_df$highlight == "highlighted"]
#     plot_df$color[is.na(plot_df$color)] <- "unknown"
#     
#     # Map colors to actual colors
#     color_mapping <- c("gray" = "gray", "long_term" = "red", "short_term" = "blue", "control" = "yellow", "unknown" = "black")
#     plot_df$color <- color_mapping[plot_df$color]
#     
#     # Create UMAP plot
#     p_umap <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
#       geom_point(aes(color = color), size = 0.5, alpha = ifelse(plot_df$color == "gray", 0.5, 1)) +
#       scale_color_identity() +
#       coord_fixed() +
#       scale_x_continuous(limits = x_limits) +
#       scale_y_continuous(limits = y_limits) +
#       theme_void() +
#       theme(legend.position = "none",
#             plot.title = element_blank(),
#             axis.title = element_blank(),
#             panel.border = element_rect(colour = "black", fill=NA, size=1),
#             panel.background = element_blank(),
#             axis.line = element_line(colour = "black"),
#             plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))
#     
#     # Store the UMAP plot
#     umap_plot_list[[plot_counter]] <- p_umap + coord_flip()
#     
#     # Create Density Plot with Gray Background
#     # Remove gray cells for density plot
#     density_df <- plot_df[plot_df$color != "gray", ]
#     
#     # Start with gray background
#     p_density <- ggplot() +
#       geom_point(data = plot_df[plot_df$color == "gray", ], 
#                  aes(x = UMAP_1, y = UMAP_2), 
#                  color = "gray", size = 0.5, alpha = 0.5)
#     
#     # Add filled density plots for each group
#     groups_present <- unique(density_df$SurvivorGroup)
#     groups_present <- groups_present[!is.na(groups_present)]
#     
#     for (grp in groups_present) {
#       grp_data <- density_df[density_df$SurvivorGroup == grp, ]
#       grp_color <- color_mapping[[grp]]
# 
#       if (nrow(grp_data) > 2) {  # Need at least 3 points for density estimation
#         p_density <- p_density +
#           stat_density_2d(data = grp_data,
#                           aes(x = UMAP_1, y = UMAP_2, fill = ..level..),
#                           geom = "polygon",
#                           bins = 10,
#                           alpha = 0.5,
#                           color = NA,
#                           fill = grp_color,
#                           contour = TRUE,
#                           inherit.aes = FALSE)
#       }
#     }
#     
#     p_density <- p_density +
#       coord_fixed() +
#       scale_x_continuous(limits = x_limits) +
#       scale_y_continuous(limits = y_limits) +
#       theme_void() +
#       theme(
#         legend.position = "none",
#         plot.title = element_blank(),
#         axis.title = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")
#       )
#     
#     # Store the density plot
#     density_plot_list[[plot_counter]] <- p_density + coord_flip()
#     plot_counter <- plot_counter + 1
#   }
#   
#   message("Total number of UMAPs and density plots: ", plot_counter - 1)
#   
#   # Check if we have any plots
#   if (length(umap_plot_list) == 0 || length(density_plot_list) == 0) {
#     message("No plots to save for subpopulation: ", subpop_name)
#     return()
#   }
#   
#   # Arrange UMAP plots into grid
#   n_rows <- 1
#   n_cols <- length(umap_plot_list)
#   
#   umap_grid_combined <- plot_grid(plotlist = umap_plot_list, ncol = n_cols)
#   
#   # Save the UMAP grid to a file
#   output_file_umap <- file.path(output_folder, paste0(subpop_name, "_clone_tracking_umap_all_patients_all_cells_fixed_axes.pdf"))
#   
#   ggsave(filename = output_file_umap, plot = umap_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
#   
#   message("Saved UMAP grid for subpopulation: ", subpop_name)
#   
#   # Arrange Density plots into grid
#   density_grid_combined <- plot_grid(plotlist = density_plot_list, ncol = n_cols)
#   
#   # Save the Density grid to a file
#   output_file_density <- file.path(output_folder, paste0(subpop_name, "_filled_density_plots_with_background_all_patients_all_cells_fixed_axes.pdf"))
#   
#   ggsave(filename = output_file_density, plot = density_grid_combined, device = "pdf", width = 4 * n_cols, height = 4 * n_rows, limitsize = FALSE)
#   
#   message("Saved density plot grid for subpopulation: ", subpop_name)
# }



# Paths to files (replace with your actual file paths)
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
base_output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/new/temp"
# Sample clone file path
sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/clonotype_df_proportion_all_T_cells.csv"


# Loop over each subpopulation and perform the analysis
for (subpop_name in names(mapping)) {
  clusters <- mapping[[subpop_name]]
  
  output_folder <- file.path(base_output_dir, subpop_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Call the function
  tryCatch({
    createCloneTrackingUMAPs(subpop_name, clusters, full_seurat_obj, clone_cell_barcode_file, sample_clone_file, output_folder)
  }, error = function(e) {
    message("Error in processing subpopulation: ", subpop_name)
    message("Error message: ", e$message)
  })
}