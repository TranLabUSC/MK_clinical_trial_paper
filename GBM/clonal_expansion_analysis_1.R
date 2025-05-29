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


# read seurat object
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
# only considering UF data as WUSTL TCR data is not good quality
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF" & Patient != 1)
seurat_metadata <- seurat_object_t_cells@meta.data
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]

# celltype_to_cluster <- data.frame(
#   celltype = c("Activated_CD4", "Activated_CD4", "Activated_CD4", "Activated_CD4", "Activated_CD4", "Effector_CD8", "Memory_CD4", "Memory_CD4", "Memory_CD4", "Naive_CD4", "Transitional_CD8", "Transitional_CD8", "Memory_CD8", "Naive_CD8", "Anergic_CD8", "Transitional_CD4", "Exhausted_T", "Proliferating_Effector_CD8", "Proliferating_Effector_CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8"),
#   cluster = c(0, 1, 11, 14, 22, 2, 3, 4, 7, 5, 6, 16, 8, 9, 10, 12, 15, 17, 18, 2, 6, 16, 8, 9, 10, 17, 18, 21)
# )

# Define the mapping
mapping <- list(
  "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
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
  print(celltype)
  cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
  seurat_object_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
  barcodes <- rownames(seurat_object_subset@meta.data)
  clone_cell_barcode_df_subset <- clone_cell_barcode_df[clone_cell_barcode_df$barcode %in% barcodes,]

  # generate a list of list of lists.
  final_list <- list()
  for (i in 1:nrow(clone_cell_barcode_df_subset)) {
    current_row <- clone_cell_barcode_df_subset[i, ]
    barcode <- current_row$barcode
    aa_sequence <- current_row$cdr3
    origin <- current_row$origin
    if (origin %in% names(final_list)){
      print("sample exists")
    } else {
      final_list[[origin]] <- list()
      final_list[[origin]][["barcodes"]] <- c()
    }
    if (aa_sequence %in% names(final_list[[origin]])){
      print("AA sequence exists")
    } else {
      final_list[[origin]][[aa_sequence]] <- 0
    }

    # update values in list
    final_list[[origin]][["barcodes"]] <- c(final_list[[origin]][["barcodes"]], barcode)
    final_list[[origin]][[aa_sequence]] <- final_list[[origin]][[aa_sequence]] + 1
  }

  # Function to convert nested list to a dataframe
  convert_to_df <- function(final_list) {
    clone_df <- data.frame()
    for (sample_name in names(final_list)) {
      print(sample_name)
      sample_list <- final_list[[sample_name]]
      sample_list$barcodes <- NULL
      temp_df <- as.data.frame(as.list(sample_list), stringsAsFactors = FALSE)
      temp_df <- as.data.frame(t(temp_df))
      colnames(temp_df) <- sample_name
      clone_df <- merge(clone_df, temp_df, by="row.names", all = TRUE)
      rownames(clone_df) <- clone_df$Row.names
      clone_df$Row.names <- NULL
    }
    return(clone_df)
  }

  clonotype_df_absolute <- convert_to_df(final_list)
  clonotype_df_absolute[is.na(clonotype_df_absolute)] <- 0
  clonotype_df_absolute$CDR3.aa <- rownames(clonotype_df_absolute)
  write.csv(clonotype_df_absolute, paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Expansion_Results/clonotype_df_absolute_", celltype, "_T_cells.csv"))

  numeric_cols <- sapply(clonotype_df_absolute, is.numeric)
  print(numeric_cols)
  print(colnames(clonotype_df_absolute)[numeric_cols])
  # Calculate sums of numeric columns and check if they are less than 50
  cols_to_remove <- sapply(clonotype_df_absolute[, numeric_cols, drop = FALSE], sum) < 1
  # Subset the dataframe to remove those columns
  clonotype_df_absolute <- clonotype_df_absolute[, !cols_to_remove]

  # Function to calculate proportions for all numeric columns
  calculate_proportions_all_numeric <- function(df) {
    numeric_cols <- sapply(df, is.numeric)

    # Copy the original dataframe to keep non-numeric data intact
    new_df <- df

    # Divide each numeric column by its sum
    new_df[, numeric_cols] <- sapply(df[, numeric_cols, drop = FALSE], function(x) x / sum(x))

    # Check if the sum of all numeric columns is now 1
    print(sapply(new_df[, numeric_cols, drop = FALSE], sum))
    return(new_df)
  }

  # Calculate proportions for all numeric columns
  proportion_df <- calculate_proportions_all_numeric(clonotype_df_absolute)
  write.csv(proportion_df, paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Expansion_Results/clonotype_df_proportion_", celltype, "_T_cells.csv"))
}


# ###########################################################################################################################################
# # clonal expansion for different clones (instead of cell types)
# 
# # subsetting for CD8 T cells
# seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
# seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters %in% c(2, 6, 16, 8, 9, 10, 17, 18, 21))
# # only considering UF data as WUSTL TCR data is not good quality
# seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
# seurat_metadata <- seurat_object_t_cells@meta.data
# clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate/outs/vdj_t/filtered_contig_annotations.csv"
# clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
# clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
# clone_cell_barcode_df <- clone_cell_barcode_df[, c("barcode", "cdr3")]
# 
# 
# # Run the table command on the specified column
# value_freq_table <- table(clone_cell_barcode_df$cdr3)
# # Convert the table output to a dataframe
# value_freq_df <- as.data.frame(value_freq_table)
# # Rename the columns for clarity
# colnames(value_freq_df) <- c("cdr3", "Frequency")
# 
# clone_cell_barcode_df <- merge(clone_cell_barcode_df, value_freq_df, by = "cdr3")
# # Filter the dataframe to keep the row with higher frequency for duplicated barcodes
# clone_cell_barcode_df <- clone_cell_barcode_df %>%
#   group_by(barcode) %>%
#   arrange(desc(Frequency)) %>%
#   slice(1) %>%
#   ungroup()
# 
# seurat_object_tcr_cells <- subset(seurat_object_t_cells, cells = clone_cell_barcode_df$barcode)
# 
# seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, clone_cell_barcode_df, by.x = "row.names", by.y = "barcode")
# rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells@meta.data
# 
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
#   group_by(origin, cdr3) %>%
#   mutate(clone_frequency = n()) %>%
#   ungroup()
# 
# # Calculate the total number of cells in each origin
# total_cells_in_origin <- seurat_object_tcr_cells_metadata %>%
#   group_by(origin) %>%
#   summarise(total_cells = n())
# 
# # Join the total cell count back to the metadata
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
#   left_join(total_cells_in_origin, by = "origin")
# 
# # Calculate clone proportion
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
#   mutate(clone_proportion = clone_frequency / total_cells) %>%
#   select(-total_cells)  # Optionally remove the total_cells column
# 
# seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, seurat_object_tcr_cells_metadata[, c("Row.names", "clone_frequency", "clone_proportion")], by.x = "row.names", by.y = "Row.names")
# rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
# seurat_object_tcr_cells_metadata <- as.data.frame(seurat_object_tcr_cells_metadata)
# rownames(seurat_object_tcr_cells_metadata) <- seurat_object_tcr_cells_metadata$Row.names
# 
# # # Create a histogram to visualize the distribution of clone proportions
# # ggplot(seurat_object_tcr_cells_metadata, aes(x = clone_proportion)) +
# #   geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
# #   labs(title = "Distribution of Clone Proportions",
# #        x = "Clone Proportion",
# #        y = "Count") +
# #   theme_minimal() +
# #   scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
# 
# 
# # # 0.01 is the threshold for our frequent clones
# # high_proportion_clone_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$clone_proportion > 0.01,])
# # 
# # # Remove graphs from the Seurat object
# # seurat_object_tcr_cells@graphs <- list()
# # 
# # # Remove neighbor and other related slots if present
# # seurat_object_tcr_cells@neighbors <- list()
# # seurat_object_tcr_cells@misc <- list()
# # # Subset the Seurat object with the ordered cell names
# # seurat_object_tcr_cells_subset <- subset(seurat_object_tcr_cells, cells = high_proportion_clone_cells)
# 
# # subset seurat object for cell barcodes where barcodes match the TCR barcodes and tag cell barcodes as "shared" or "not shared"
# # map TCR barcodes
# blood_clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate/outs/vdj_t/filtered_contig_annotations.csv"
# blood_clone_cell_barcode_df <- read.csv(blood_clone_cell_barcode_file)
# cell_barcodes <- unique(blood_clone_cell_barcode_df$barcode)
# 
# seurat_object_tcr_barcodes <- subset(seurat_object_t_cells, cells = cell_barcodes)
# 
# # add an attribute to the seurat object metadata, "clone_status"
# # Initialize the clone_status column in metadata
# seurat_object_tcr_barcodes@meta.data$clone_proportion_status <- "non_dominant"
# 
# # Get the list of donors
# donors <- unique(seurat_object_tcr_barcodes@meta.data$donor)
# 
# # 0.01 is the threshold for our frequent clones
# high_proportion_clone_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$clone_proportion > 0.01,])
# 
# # Iterate over each donor
# for (donor in donors) {
#   # Get the cell barcodes for the current donor
#   donor_cells <- rownames(seurat_object_tcr_barcodes@meta.data[seurat_object_tcr_barcodes@meta.data$donor == donor, ])
#   patient_id <- donor
#   # Find the shared cell barcodes for the current donor
#   hp_cells <- intersect(donor_cells, high_proportion_clone_cells)
#   
#   # Update the clone_status for shared cells
#   seurat_object_tcr_barcodes@meta.data[hp_cells, "clone_proportion_status"] <- "dominant"
# }
# 
# 
# get_tumor_tcr_clone_cells <- function(tumor_sample_clone_file, blood_sample_clone_file, tumor_tcr_metadata_df, blood_tcr_metadata_file, blood_clone_cell_barcode_file, seurat_obj) {
#   
#   tumor_tcr_clone_cells <- list()
#   # Read files
#   blood_tcr_metadata_df <- read.csv(blood_tcr_metadata_file)
#   tumor_sample_clone_df <- read.csv(tumor_sample_clone_file, sep = "\t", check.names = FALSE)
#   blood_sample_clone_df <- read.csv(blood_sample_clone_file, sep = "\t", check.names = FALSE)
#   blood_clone_cell_barcode_df <- read.csv(blood_clone_cell_barcode_file)
#   
#   # Create the 'timepoint' column by removing the donor prefix from the origin
#   blood_tcr_metadata_df$Timepoint <- mapply(function(donor, origin) {
#     sub(paste0("^", donor), "", origin)
#   }, blood_tcr_metadata_df$donor, blood_tcr_metadata_df$origin)
#   # Replace "S" with "Pre" in the timepoint column
#   blood_tcr_metadata_df$Timepoint <- gsub("S", "Pre", blood_tcr_metadata_df$Timepoint)
#   
#   ordered_timepoints <- c("R1", "R2", "Pre", "C1", "C2", "C4",  "C6", "C9", "C18", "C36")
#   
#   # Iterate over each patient
#   for (patient in unique(tumor_tcr_metadata_df$Patient)) {
#     print(patient)
#     # for (patient in c("p9", "p19")) {
#     patient_metadata_tumor <- tumor_tcr_metadata_df[tumor_tcr_metadata_df$Patient == patient, ]
#     patient_metadata_blood <- blood_tcr_metadata_df[blood_tcr_metadata_df$Patient == patient, ]
#     
#     if (nrow(patient_metadata_blood) == 0){
#       next
#     }
#     
#     if (nrow(patient_metadata_tumor) == 0){
#       next
#     }
#     
#     patient_specific_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$donor == as.integer(patient),])
#     if (length(patient_specific_cells) == 0) {
#       next
#     }
#     patient_seurat_obj <- subset(seurat_obj, cells = patient_specific_cells)
#     
#     patient_timepoint_tumor <- intersect(ordered_timepoints, patient_metadata_tumor$Timepoint)
#     patient_timepoint_blood <- intersect(ordered_timepoints, patient_metadata_blood$Timepoint)
#     
#     tumor_tcr_clone_cells[[patient]] <- list()
#     
#     # Iterate over each timepoint
#     for (timepoint in patient_timepoint_tumor) {
#       print(timepoint)
#       # print(patient_timepoint_tumor[patient_timepoint_tumor$integrated_sample_type == timepoint])
#       sample_ids <- patient_metadata_tumor$Sample[patient_metadata_tumor$Timepoint == timepoint]
#       
#       if (length(sample_ids) > 0) {
#         # Identify all clones for the sample_id
#         all_clones <- which(tumor_sample_clone_df[, sample_ids] > 0)
#         clone_names <- tumor_sample_clone_df$CDR3.aa[all_clones]
#         
#         # Identify cell barcodes for these clones
#         cell_barcodes <- blood_clone_cell_barcode_df$barcode[blood_clone_cell_barcode_df$cdr3 %in% clone_names]
#         
#         tumor_tcr_clone_cells[[patient]][[timepoint]] <- list()
#         # Iterate over each timepoint in timepoint_mapping
#         for (tp in patient_timepoint_blood) {
#           print(tp)
#           tp_cells <- WhichCells(patient_seurat_obj, expression = TimePoint == tp)
#           # Plot UMAP
#           if (length(tp_cells) > 0) {
#             
#             patient_tp_seurat_obj <- subset(patient_seurat_obj, cells = tp_cells)
#             tcr_clone_cells <- intersect(cell_barcodes, rownames(patient_tp_seurat_obj@meta.data))
#             tumor_tcr_clone_cells[[patient]][[timepoint]][[tp]] <- tcr_clone_cells
#           }
#         }
#       }
#     }
#   }
#   return(tumor_tcr_clone_cells)
# }
# 
# tumor_sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt"
# blood_sample_clone_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt"
# tumor_tcr_metadata_df <- data.frame(
#   Patient = c("1", "1", "3", "4", "5", 
#               "9", "9", "10", "11", "13", 
#               "14", "18", "19", "20"),
#   Sample = c("MK1_TCRB", "MK2_TCRB", "MK3_TCRB", "MK4_TCRB", "M5F_TCRB", "M5F_TCRB", "MK7_TCRB", "MK8_TCRB", "MK9_TCRB", "M10F_TCRB", 
#              "MK11_TCRB", "MK12_TCRB", "M13_TCRB", "M14_TCRB"),
#   Timepoint = c("R1", "R2", "R1", "R1", "R1", "R1", "R2", "R1", "R1", "R1", 
#                 "R1", "R1", "R1", "R1")
# )
# blood_tcr_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate.csv"
# blood_clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate/outs/vdj_t/filtered_contig_annotations.csv"
# tumor_tcr_clone_cells <- get_tumor_tcr_clone_cells(tumor_sample_clone_file, blood_sample_clone_file, tumor_tcr_metadata_df, blood_tcr_metadata_file, blood_clone_cell_barcode_file, seurat_object_t_cells)
# 
# # add an attribute to the seurat object metadata, "clone_status"
# # Initialize the clone_status column in metadata
# seurat_object_tcr_barcodes@meta.data$clone_share_status <- "non_shared"
# 
# # Get the list of donors
# donors <- unique(seurat_object_tcr_barcodes@meta.data$donor)
# 
# get_shared_clone_cell_barcodes_for_patient <- function(patient_id, tumor_tcr_clone_cells){
#   cells <- c()
#   patient_tumor_tcr_clone_cells <- tumor_tcr_clone_cells[[patient_id]]
#   for (tumor_timepoint in names(patient_tumor_tcr_clone_cells)){
#     tumor_timepoint_tcr_clone_cells <- patient_tumor_tcr_clone_cells[[tumor_timepoint]]
#     for (blood_timepoint in names(tumor_timepoint_tcr_clone_cells)){
#       cells <- c(cells, tumor_timepoint_tcr_clone_cells[[blood_timepoint]])
#     }
#   }
#   return(unique(cells))
# }
# 
# # Iterate over each donor
# for (donor in donors) {
#   # Get the cell barcodes for the current donor
#   donor_cells <- rownames(seurat_object_tcr_barcodes@meta.data[seurat_object_tcr_barcodes@meta.data$donor == donor, ])
#   patient_id <- donor
#   shared_clone_cell_barcodes <- get_shared_clone_cell_barcodes_for_patient(patient_id, tumor_tcr_clone_cells)
#   # Find the shared cell barcodes for the current donor
#   shared_cells <- intersect(donor_cells, shared_clone_cell_barcodes)
#   
#   # Update the clone_status for shared cells
#   seurat_object_tcr_barcodes@meta.data[shared_cells, "clone_share_status"] <- "shared"
# }
# 
# # Create the new column 'clone_status' separately
# clone_status <- paste(seurat_object_tcr_barcodes@meta.data$clone_proportion_status, seurat_object_tcr_barcodes@meta.data$clone_share_status, sep = "_")
# 
# # Assign the new column back to the dataframe
# seurat_object_tcr_barcodes@meta.data$clone_status <- clone_status
# 
# 
# clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate/outs/vdj_t/filtered_contig_annotations.csv"
# clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
# clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
# 
# 
# for (clone_type in unique(seurat_object_tcr_barcodes@meta.data$clone_status)){
#   print(clone_type)
#   seurat_object_subset <- subset(seurat_object_tcr_barcodes, subset = clone_status == clone_type)
#   barcodes <- rownames(seurat_object_subset@meta.data)
#   clone_cell_barcode_df_subset <- clone_cell_barcode_df[clone_cell_barcode_df$barcode %in% barcodes,]
#   
#   # generate a list of list of lists.
#   final_list <- list()
#   for (i in 1:nrow(clone_cell_barcode_df_subset)) {
#     current_row <- clone_cell_barcode_df_subset[i, ]
#     barcode <- current_row$barcode
#     aa_sequence <- current_row$cdr3
#     origin <- current_row$origin
#     if (origin %in% names(final_list)){
#       print("sample exists")
#     } else {
#       final_list[[origin]] <- list()
#       final_list[[origin]][["barcodes"]] <- c()
#     }
#     if (aa_sequence %in% names(final_list[[origin]])){
#       print("AA sequence exists")
#     } else {
#       final_list[[origin]][[aa_sequence]] <- 0
#     }
#     
#     # update values in list
#     final_list[[origin]][["barcodes"]] <- c(final_list[[origin]][["barcodes"]], barcode)
#     final_list[[origin]][[aa_sequence]] <- final_list[[origin]][[aa_sequence]] + 1
#   }
#   
#   # Function to convert nested list to a dataframe
#   convert_to_df <- function(final_list) {
#     clone_df <- data.frame()
#     for (sample_name in names(final_list)) {
#       print(sample_name)
#       sample_list <- final_list[[sample_name]]
#       sample_list$barcodes <- NULL
#       temp_df <- as.data.frame(as.list(sample_list), stringsAsFactors = FALSE)
#       temp_df <- as.data.frame(t(temp_df))
#       colnames(temp_df) <- sample_name
#       clone_df <- merge(clone_df, temp_df, by="row.names", all = TRUE)
#       rownames(clone_df) <- clone_df$Row.names
#       clone_df$Row.names <- NULL
#     }
#     return(clone_df)
#   }
#   
#   clonotype_df_absolute <- convert_to_df(final_list)
#   clonotype_df_absolute[is.na(clonotype_df_absolute)] <- 0
#   clonotype_df_absolute$CDR3.aa <- rownames(clonotype_df_absolute)
#   write.csv(clonotype_df_absolute, paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/latest_nomenclature/clonotype_df_absolute_", clone_type, "_T_cells.csv"))
#   
#   numeric_cols <- sapply(clonotype_df_absolute, is.numeric)
#   print(numeric_cols)
#   print(colnames(clonotype_df_absolute)[numeric_cols])
#   # Calculate sums of numeric columns and check if they are less than 50
#   cols_to_remove <- sapply(clonotype_df_absolute[, numeric_cols, drop = FALSE], sum) < 50
#   # Subset the dataframe to remove those columns
#   clonotype_df_absolute <- clonotype_df_absolute[, !cols_to_remove]
#   
#   # Function to calculate proportions for all numeric columns
#   calculate_proportions_all_numeric <- function(df) {
#     numeric_cols <- sapply(df, is.numeric)
#     
#     # Copy the original dataframe to keep non-numeric data intact
#     new_df <- df
#     
#     # Divide each numeric column by its sum
#     new_df[, numeric_cols] <- sapply(df[, numeric_cols, drop = FALSE], function(x) x / sum(x))
#     
#     # Check if the sum of all numeric columns is now 1
#     print(sapply(new_df[, numeric_cols, drop = FALSE], sum))
#     return(new_df)
#   }
#   
#   # Calculate proportions for all numeric columns
#   proportion_df <- calculate_proportions_all_numeric(clonotype_df_absolute)
#   write.csv(proportion_df, paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/latest_nomenclature/clonotype_df_proportion_", clone_type, "_T_cells.csv"))
# }
