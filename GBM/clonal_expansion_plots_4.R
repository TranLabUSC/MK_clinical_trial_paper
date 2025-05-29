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

  
make_comparison_value_plot = function(df,
                                      variable_name,
                                      diversity_index = "simpson",
                                      alternative = "l",
                                      plot_dir,
                                      arm_type)
{
  #alternative	: used in t test function. a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
  print(df)
  # Filter out infinite values
  finite_values <- df[, variable_name][is.finite(df[, variable_name])]
  p_val = t.test(x = finite_values,
                 mu = 1,
                 alternative = alternative)[["p.value"]]
  p_label = paste("p =", round(p_val, 3))
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  plot_label = paste("Relative",
                     variable_name)
  x_label = paste("relative", variable_name)
  df_plot <- df[is.finite(df[, variable_name]), ]
  if (diversity_index == "shannon"){
    p = ggplot(df_plot, aes(x = x_label , y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.2) +
      geom_boxplot(width = 0.1) +
      geom_point() +
      annotate(
        "text",
        label = p_label,
        x = 1,
        y = max(df_plot[, variable_name], na.rm = TRUE) * 1.5,
        size = 10  ,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label)  +
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
        size = 10  ,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label)  +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 2
      )
  }
  
  p = p + theme_bw() +
    theme(
      legend.position = "none" ,
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
  file_name <- paste0(diversity_index, "_Clonal_Diversity_", arm_type, "_", variable_name, ".pdf")
  ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 4, height = 6, dpi = 300, limitsize = FALSE)
}

patient_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
patient_metadata_df <- read.csv(patient_metadata_file)
# celltypes = celltype = c("Activated_CD4", "Activated_CD4", "Activated_CD4", "Activated_CD4", "Activated_CD4", "Effector_CD8", "Memory_CD4", "Memory_CD4", "Memory_CD4", "Naive_CD4", "Transitional_CD8", "Transitional_CD8", "Memory_CD8", "Naive_CD8", "Anergic_CD8", "Transitional_CD4", "Exhausted_T", "Proliferating_Effector_CD8", "Proliferating_Effector_CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "CD8", "all")
# celltypes <- unique(celltypes)

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

celltypes <- unique(celltype_to_cluster$celltype)

diversity_indices = c("shannon", "simpson")
dir = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Clonal_Expansion_Results/C1_vs_C2"
for (celltype in celltypes) {
  print(celltype)
  for (diversity_index in diversity_indices){
    print(diversity_index)
    plot_dir = paste0(dir, "/clonal_expansion_", diversity_index, "_", celltype, "_T_cells_division/")
    df <- read.delim(paste0(plot_dir, "C1_vs_C2.txt"))
    df_merged <- merge(df, patient_metadata_df[, c("patient_id", "Arm")], by.x = "row.names", by.y = "patient_id", all.x = TRUE, all.y = FALSE)
    for (arm_type in unique(df_merged$Arm)){
      df <- df_merged[df_merged$Arm == arm_type, ]
      df <- df[!is.na(df$clonal_expansion), ]
      if (nrow(df) > 1){
        make_comparison_value_plot(df, "clonal_expansion", diversity_index = diversity_index, alternative = "l", plot_dir, arm_type)
      } else {
        print("=======================================================================================================================")
        print(celltype)
        print(df)
      }
    }
  }
}

######################################################################################################################################################################################################################
library(immunarch)  # Load the package into R
data(immdata)
# generate alluvial plot for Pre vs C1 for effector CD8 cells (experiment group)
clonotype_df_proportion_Effector_CD8_T_cells <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/latest_nomenclature/clonotype_df_proportion_Effector_CD8_T_cells.csv", row.names = 1, check.names = FALSE)

clonotype_df_proportion_Effector_CD8_no_threshold_T_cells <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/latest_nomenclature/clonotype_df_proportion_Effector_CD8_no_threshold_T_cells.csv", row.names = 1, check.names = FALSE)

clonotype_df_absolute_Effector_CD8_no_threshold_T_cells <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/latest_nomenclature/clonotype_df_absolute_Effector_CD8_no_threshold_T_cells.csv", row.names = 1, check.names = FALSE)

data <- clonotype_df_proportion_Effector_CD8_no_threshold_T_cells
data$Clone <- rownames(data)

abs_data <- clonotype_df_absolute_Effector_CD8_no_threshold_T_cells
abs_data$Clone <- rownames(abs_data)

# Select the subset of patients (e.g., 7, 10, 3)
patients_of_interest <- c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7")

# Step 1: Remove patients missing either C1 or S timepoints
valid_patients <- c()  # Initialize list for valid patients

for (patient in patients_of_interest) {
  c1_col <- paste0(patient, "C1")
  s_col <- paste0(patient, "S")

  # Check if both C1 and S timepoints exist for this patient
  if (c1_col %in% colnames(data) && s_col %in% colnames(data)) {
    valid_patients <- c(valid_patients, patient)
  }
}

# Update patients_of_interest to only include valid patients
patients_of_interest <- valid_patients

# patients_of_interest <- c("9")

# Filter only columns that correspond to C1 timepoints for the valid patients
c1_columns <- paste0(patients_of_interest, "C1")
c1_data <- data[, c("Clone", c1_columns)]
c1_data_abs <- abs_data[, c("Clone", c1_columns)]

# Step 2: Identify the top 1% clones by proportion for each valid patient
top_clones <- data.frame()  # Initialize an empty dataframe to store results

for (patient in patients_of_interest) {
  print(patient)
  # Get the column name for this patient's C1 timepoint
  patient_c1_col <- paste0(patient, "C1")

  # Get the 1% threshold
  # threshold <- quantile(c1_data[[patient_c1_col]], 0.99)
  proportion_threshold <- 0.01
  absolute_threshold <- 1

  # Identify the clones above this threshold
  top_proportion_patient_clones <- c1_data[c1_data[[patient_c1_col]] > proportion_threshold, c("Clone", patient_c1_col)]
  colnames(top_proportion_patient_clones) <- c("Clone", "Proportion_C1")
  top_absolute_patient_clones <- c1_data_abs[c1_data_abs[[patient_c1_col]] > absolute_threshold, c("Clone", patient_c1_col)]
  colnames(top_absolute_patient_clones) <- c("Clone", "Absolute_C1")
  top_clone_names <- intersect(top_proportion_patient_clones$Clone, top_absolute_patient_clones$Clone)
  top_patient_clones <- top_proportion_patient_clones[top_proportion_patient_clones$Clone %in% top_clone_names, ]
  print(top_patient_clones)

  if(nrow(top_patient_clones) == 0){
    next
  }
  # Add a column indicating the patient ID
  top_patient_clones$Patient <- patient

  # Append to the final result
  top_clones <- rbind(top_clones, top_patient_clones)
}

# Step 3: Get the proportion of these clones at timepoint S by looking at the "PatientS" column
top_clones$Proportion_S <- NA  # Create a new column to store the S proportions

for (i in 1:nrow(top_clones)) {
  # Get the patient ID and clone name for the current row
  patient <- top_clones$Patient[i]
  clone_name <- top_clones$Clone[i]

  # Get the column name for this patient's S timepoint
  patient_s_col <- paste0(patient, "S")

  # Find the proportion of the clone at timepoint S
  proportion_s <- data[data$Clone == clone_name, patient_s_col]

  # Store the proportion at S
  top_clones$Proportion_S[i] <- proportion_s
}

num_unique_patients <- length(unique(top_clones$Patient))

# Data for C1
c1_data <- data.frame(
  CDR3.aa = top_clones$Clone,
  Clones = top_clones$Proportion_C1*100,
  Proportion = top_clones$Proportion_C1,
  Sample = "C1"
  # Sample = paste0(top_clones$Patient, "_C1")  # Labeling each sample with patient and timepoint
)

# New row to add
new_row <- data.frame(
  CDR3.aa = 'Others',
  Clones = (num_unique_patients - sum(c1_data$Proportion))*100,
  Proportion = num_unique_patients - sum(c1_data$Proportion),
  Sample = "C1"
)

# Add the new row to the existing data frame
c1_data <- rbind(c1_data, new_row)

# Data for S
s_data <- data.frame(
  CDR3.aa = top_clones$Clone,
  Clones = top_clones$Proportion_S*100,
  Proportion = top_clones$Proportion_S,
  Sample = "S"
  # Sample = paste0(top_clones$Patient, "_S")  # Labeling each sample with patient and timepoint
)

# New row to add
new_row <- data.frame(
  CDR3.aa = 'Others',
  Clones = (num_unique_patients - sum(s_data$Proportion))*100,
  Proportion = num_unique_patients - sum(s_data$Proportion),
  Sample = "S"
)

# Add the new row to the existing data frame
s_data <- rbind(s_data, new_row)

# Combine both into a single data frame
immunarch_data <- rbind(s_data, c1_data)

# Step 2: Convert to the immunarch data format
# Create a named list with sample names, where each entry is a sample-specific repertoire
immdata_new <- split(immunarch_data[, c("Clones", "CDR3.aa", "Proportion")], immunarch_data$Sample)

# Step 3: Prepare the metadata as required by immunarch (you can omit patient info)
meta_data <- data.frame(
  Sample = unique(immunarch_data$Sample),
  Timepoint = unique(immunarch_data$Sample)  # Only "C1" and "S" timepoints
)

# Construct the final `immdata` object required by immunarch
immdata_new <- list(
  data = immdata_new,  # List of clone proportions by sample
  meta = meta_data  # Metadata including timepoint info
)

target <- s_data$CDR3.aa[s_data$CDR3.aa != "Others"]
tc <- trackClonotypes(immdata_new$data, target, .col = "aa")

vis(tc, .order = c(2,1), .plot = "smooth")
#######################################################################################################################################################################################################################
# read seurat object
  seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
  # only considering UF data as WUSTL TCR data is not good quality
  seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
  seurat_metadata <- seurat_object_t_cells@meta.data
  clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/Aggregate/MK_aggregate/outs/vdj_t/filtered_contig_annotations.csv"
  clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
  clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
  clone_cell_barcode_df <- clone_cell_barcode_df[, c("barcode", "cdr3")]
  
  
  # Run the table command on the specified column
  value_freq_table <- table(clone_cell_barcode_df$cdr3)
  # Convert the table output to a dataframe
  value_freq_df <- as.data.frame(value_freq_table)
  # Rename the columns for clarity
  colnames(value_freq_df) <- c("cdr3", "Frequency")
  
  clone_cell_barcode_df <- merge(clone_cell_barcode_df, value_freq_df, by = "cdr3")
  # Filter the dataframe to keep the row with higher frequency for duplicated barcodes
  clone_cell_barcode_df <- clone_cell_barcode_df %>%
    group_by(barcode) %>%
    arrange(desc(Frequency)) %>%
    slice(1) %>%
    ungroup()
  
  seurat_object_tcr_cells <- subset(seurat_object_t_cells, cells = clone_cell_barcode_df$barcode)
  
  seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, clone_cell_barcode_df, by.x = "row.names", by.y = "barcode")
  rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
  seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells@meta.data
  
  seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
    group_by(origin, cdr3) %>%
    mutate(clone_frequency = n()) %>%
    ungroup()
  
  # Calculate the total number of cells in each origin
  total_cells_in_origin <- seurat_object_tcr_cells_metadata %>%
    group_by(origin) %>%
    summarise(total_cells = n())
  
  # Join the total cell count back to the metadata
  seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
    left_join(total_cells_in_origin, by = "origin")
  
  # Calculate clone proportion
  seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
    mutate(clone_proportion = clone_frequency / total_cells) %>%
    select(-total_cells)  # Optionally remove the total_cells column
  
  seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, seurat_object_tcr_cells_metadata[, c("Row.names", "clone_frequency", "clone_proportion")], by.x = "row.names", by.y = "Row.names")
  rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
  seurat_object_tcr_cells_metadata <- as.data.frame(seurat_object_tcr_cells_metadata)
  rownames(seurat_object_tcr_cells_metadata) <- seurat_object_tcr_cells_metadata$Row.names
  
  
  control_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4"), ])
  exp_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7"), ])
  sampled_exp_cells <- sample(exp_cells, length(control_cells))
  
  control_cells_C1 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4") & seurat_object_tcr_cells_metadata$TimePoint == "C1", ])
  exp_cells_C1 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7") & seurat_object_tcr_cells_metadata$TimePoint == "C1", ])
  sampled_control_cells_C1 <- sample(control_cells_C1, 2183)
  sampled_exp_cells_C1 <- sample(exp_cells_C1, 2183)
  
  control_cells_C2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4") & seurat_object_tcr_cells_metadata$TimePoint == "C2", ])
  exp_cells_C2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7") & seurat_object_tcr_cells_metadata$TimePoint == "C2", ])
  sampled_exp_cells_C2 <- sample(exp_cells_C2, 2183)

# # Generate UMAP plot colored by the 'frequency' column
# umap_plot <- FeaturePlot(seurat_object_tcr_cells, features = "Frequency", reduction = "umap", cells = sampled_exp_cells_C2) +
#   scale_color_gradientn(colors = c("blue", "yellow", "red"), 
#                         values = scales::rescale(c(1, 5, 10, 20, 50)),
#                         breaks = c(1, 5, 10, 20, 50))

seurat_object_tcr_cells@meta.data$clone_proportion_capped <- pmin(seurat_object_tcr_cells@meta.data$clone_proportion, 0.2)
# Generate UMAP plot colored by the 'frequency' column
umap_plot <- FeaturePlot(seurat_object_tcr_cells, features = "clone_proportion_capped", reduction = "umap", cells = sampled_exp_cells) +
  scale_color_gradientn(colors = c("blue", "yellow", "red"), 
                        values = scales::rescale(c(0, 0.05, 0.10, 0.15, 0.20)),
                        breaks = c(0, 0.05, 0.10, 0.15, 0.20))

umap_plot + coord_flip()

umap_plot <- FeaturePlot(seurat_object_tcr_cells, features = "clone_proportion_capped", reduction = "umap", cells = sampled_exp_cells_C2, pt.size =2) +
  scale_color_gradientn(colors = c("blue", "yellow", "red"), 
                        values = scales::rescale(c(0, 0.05, 0.10, 0.15, 0.20)),
                        breaks = c(0, 0.05, 0.10, 0.15, 0.20))

umap_plot + coord_flip()


FeaturePlot(seurat_object_tcr_cells, features = "TRAV12-2", reduction = "umap", cells = exp_cells_C2, pt.size =2, cols = c("white", "red")) + coord_flip()

control_cells_C2_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4") & seurat_object_tcr_cells_metadata$TimePoint == "C2" & seurat_object_tcr_cells_metadata$seurat_clusters == "2", ])
exp_cells_C2_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7") & seurat_object_tcr_cells_metadata$TimePoint == "C2" & seurat_object_tcr_cells_metadata$seurat_clusters == "2", ])
# sampled_exp_cells_C2 <- sample(exp_cells_C2, 2260)
seurat_object_tcr_cells@meta.data$group <- NA
# Annotate group1 cells
seurat_object_tcr_cells@meta.data$group[rownames(seurat_object_tcr_cells@meta.data) %in% control_cells_C2_cluster_2] <- "group1"

# Annotate group2 cells
seurat_object_tcr_cells@meta.data$group[rownames(seurat_object_tcr_cells@meta.data) %in% exp_cells_C2_cluster_2] <- "group2"

table(seurat_object_tcr_cells@meta.data$group)

seurat_object <- seurat_object_tcr_cells
seurat_object <- DietSeurat(seurat_object, graphs = NULL)
Idents(seurat_object_tcr_cells) <- seurat_object_tcr_cells@meta.data$group

# Get the cell names for group1 and group2
group1_cells <- WhichCells(seurat_object, ident = "group1")
group2_cells <- WhichCells(seurat_object, ident = "group2")

# Subset the Seurat object using the cell names
seurat_object_subset <- subset(seurat_object, cells = c(group1_cells, group2_cells))

Idents(seurat_object_subset) <- seurat_object_subset@meta.data$group

# Perform differential expression analysis
markers <- FindMarkers(seurat_object_tcr_cells, ident.1 = "group1", ident.2 = "group2")

# View the top significant genes
head(markers)



library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# # Function to create violin plot, calculate p-value, and annotate the plot
# create_violin_plot <- function(seurat_object, cell_list_1, cell_list_2, metadata_column, output_folder, cluster) {
#   # Check for duplicate column names and make them unique
#   colnames(seurat_object@meta.data) <- make.names(colnames(seurat_object@meta.data), unique = TRUE)
#   # Initialize the group column with NA
#   seurat_object@meta.data$group <- NA
#   
#   # Annotate group1 cells
#   seurat_object@meta.data$group[rownames(seurat_object@meta.data) %in% cell_list_1] <- "group1"
#   
#   # Annotate group2 cells
#   seurat_object@meta.data$group[rownames(seurat_object@meta.data) %in% cell_list_2] <- "group2"
#   
#   # Verify the annotations
#   print(table(seurat_object@meta.data$group))
#   
#   # Subset the Seurat object to include only the annotated cells
#   # seurat_object <- subset(seurat_object, subset = !is.na(group))
#   
#   # Extract the metadata values for group1 and group2
#   group1_values <- seurat_object@meta.data %>% filter(group == "group1") %>% pull(!!sym(metadata_column))
#   group2_values <- seurat_object@meta.data %>% filter(group == "group2") %>% pull(!!sym(metadata_column))
#   
#   # Perform Wilcoxon test
#   p_value <- wilcox.test(group1_values, group2_values)$p.value
#   
#   # Create the violin plot, filtering out NA values
#   violin_plot <- ggplot(seurat_object@meta.data %>% filter(!is.na(group)), aes(x = group, y = !!sym(metadata_column))) +
#     geom_violin(aes(fill = group), trim = FALSE) +
#     geom_jitter(shape = 16, position = position_jitter(0.2), size = 1) +
#     theme_minimal() +
#     labs(y = metadata_column, x = "Group") +
#     ggtitle(paste("Violin Plot of", metadata_column, "in Group 1 and Group 2")) +
#     annotate("text", x = 1.5, y = max(c(group1_values, group2_values), na.rm = TRUE), 
#              label = paste("P-value:", format(p_value, digits = 3)), size = 4, vjust = -1)
#   
#   # Save the plot as a PDF
#   output_file <- file.path(output_folder, paste0("violin_plot_", metadata_column, ".pdf"))
#   ggsave(output_file, plot = violin_plot, device = "pdf")
#   
#   # Print the p-value
#   print(paste("P-value:", p_value))
# }

create_violin_plot <- function(seurat_object, cell_list_1, cell_list_2, metadata_column, output_folder, cluster, patient_metadata_df) {
  # Check for duplicate column names and make them unique
  colnames(seurat_object@meta.data) <- make.names(colnames(seurat_object@meta.data), unique = TRUE)
  
  # Initialize the group column with NA
  seurat_object@meta.data$group <- NA
  
  # Annotate group1 cells
  seurat_object@meta.data$group[rownames(seurat_object@meta.data) %in% cell_list_1] <- "group1"
  
  # Annotate group2 cells
  seurat_object@meta.data$group[rownames(seurat_object@meta.data) %in% cell_list_2] <- "group2"
  
  # Verify the annotations
  print(table(seurat_object@meta.data$group))
  
  # Merge with patient_metadata_df to get OS.months
  seurat_object@meta.data <- seurat_object@meta.data %>%
    left_join(patient_metadata_df, by = c("Patient" = "patient_id"))
  
  # Order patients by decreasing OS.months
  patients_ordered <- seurat_object@meta.data %>%
    distinct(Patient, .keep_all = TRUE) %>%
    arrange(desc(OS.months.)) %>%
    pull(Patient)
  
  # Create a list to store individual plots
  plot_list <- list()
  
  seurat_metadata_df <- as.data.frame(seurat_object@meta.data)
  print(patients)
  # Iterate over each patient and create a violin plot
  for (patient in patients_ordered) {
    print(patient)
    # Subset the data for the current patient
    patient_data <- seurat_metadata_df %>% filter(Patient == patient, !is.na(group))
    if (nrow(patient_data) == 0) next
    print(head(patient_data))
    
    # Identify common cdr3 values in both groups
    common_cdr3 <- intersect(
      patient_data %>% filter(group == "group1") %>% pull(cdr3),
      patient_data %>% filter(group == "group2") %>% pull(cdr3)
    )
    
    # Add a new column for coloring based on cdr3 presence in both groups
    patient_data <- patient_data %>%
      mutate(cdr3_color = ifelse(cdr3 %in% common_cdr3, cdr3, "black"))
    
    # Extract the metadata values for group1 and group2
    group1_values <- patient_data %>% filter(group == "group1") %>% pull(!!sym(metadata_column))
    group2_values <- patient_data %>% filter(group == "group2") %>% pull(!!sym(metadata_column))
    
    # Perform Wilcoxon test
    p_value <- wilcox.test(group1_values, group2_values)$p.value
    
    # Create the violin plot
    violin_plot <- ggplot(patient_data, aes(x = group, y = !!sym(metadata_column))) +
      geom_violin(aes(fill = group), trim = FALSE) +
      geom_jitter(aes(color = cdr3_color), shape = 16, position = position_jitter(0.2), size = 1) +
      theme_minimal() +
      labs(y = metadata_column, x = "Group") +
      ggtitle(paste("Patient:", patient)) +
      annotate("text", x = 1.5, y = max(c(group1_values, group2_values), na.rm = TRUE), 
               label = paste("P-value:", format(p_value, digits = 3)), size = 4, vjust = -1) +
      scale_color_manual(values = setNames(c("black", rainbow(length(common_cdr3))), c("black", common_cdr3)))
    
    # Add the plot to the list
    plot_list[[patient]] <- violin_plot
  }
  
  # Combine all plots into a single PDF
  output_file <- file.path(output_folder, paste0("violin_plots_C1_vs_C2_Cluster_2_control_group_", metadata_column, ".pdf"))
  pdf(output_file, width = 12*length(plot_list), height = 7)
  print(plot_grid(plotlist = plot_list, ncol = length(plot_list), nrow = 1))
  dev.off()
  
  # Print the p-values for each patient
  print(paste("P-values for each patient:"))
  for (patient in names(plot_list)) {
    print(paste(patient, ":", wilcox.test(
      seurat_object@meta.data %>% filter(Patient == patient, group == "group1") %>% pull(!!sym(metadata_column)),
      seurat_object@meta.data %>% filter(Patient == patient, group == "group2") %>% pull(!!sym(metadata_column))
    )$p.value))
  }
}

# Example usage
# Replace seurat_object, cell_list_1, cell_list_2, "feature_name", and "output_folder" with your specific values
cluster = "2"
seurat_object <- seurat_object_tcr_cells
control_cells_C1_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4") & seurat_object_tcr_cells_metadata$TimePoint == "C1" & seurat_object_tcr_cells_metadata$seurat_clusters == cluster, ])
control_cells_C2_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4") & seurat_object_tcr_cells_metadata$TimePoint == "C2" & seurat_object_tcr_cells_metadata$seurat_clusters == cluster, ])
exp_cells_C1_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7") & seurat_object_tcr_cells_metadata$TimePoint == "C1" & seurat_object_tcr_cells_metadata$seurat_clusters == cluster, ])
exp_cells_C2_cluster_2 <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7") & seurat_object_tcr_cells_metadata$TimePoint == "C2" & seurat_object_tcr_cells_metadata$seurat_clusters == cluster, ])
sampled_exp_cells_C2_cluster_2 <- sample(exp_cells_C2_cluster_2, 436)
metadata_column = "clone_proportion"
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/new_nomenclature/violin_plots/" 
patient_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
patient_metadata_df <- read.csv(patient_metadata_file)
create_violin_plot(seurat_object, control_cells_C1_cluster_2, control_cells_C2_cluster_2, metadata_column, output_folder, cluster, patient_metadata_df)






plot_top_n_clones <- function(seurat_object, top_n, output_folder) {
  # Check for duplicate column names and make them unique
  colnames(seurat_object@meta.data) <- make.names(colnames(seurat_object@meta.data), unique = TRUE)
  # Extract unique patient IDs
  patients <- unique(seurat_object@meta.data$Patient)
  
  # Create a list to store individual plots
  plot_list <- list()
  
  # Iterate over each patient
  for (patient in patients) {
    # Define origins for C1 and C2
    origin_C2 <- paste0(patient, "C2")
    origin_C1 <- paste0(patient, "C1")
    print(origin_C1)
    print(origin_C2)
    # Subset data for the current patient and origin C2
    patient_data_C2 <- seurat_object@meta.data %>% filter(Patient == patient & origin == origin_C2)
    print("here")
    if(nrow(patient_data_C2) == 0) {
      next
    }
    
    # Identify top_n clones with the highest clone_proportion in origin C2
    top_clones <- patient_data_C2 %>%
      group_by(cdr3) %>%
      summarize(clone_proportion = max(clone_proportion)) %>%
      top_n(n = top_n, wt = clone_proportion) %>%
      pull(cdr3)
    
    # Subset data for the current patient and origin C1
    patient_data_C1 <- seurat_object@meta.data %>% filter(Patient == patient & origin == origin_C1)
    if(nrow(patient_data_C1) == 0) {
      next
    }
    
    # Create a data frame to store clone proportions for C1 and C2
    clone_proportions <- data.frame(
      cdr3 = top_clones,
      clone_proportion_C1 = sapply(top_clones, function(clone) {
        prop <- patient_data_C1 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
        if (length(prop) == 0) 0 else prop[1]
      }),
      clone_proportion_C2 = sapply(top_clones, function(clone) {
        prop <- patient_data_C2 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
        if (length(prop) == 0) 0 else prop[1]
      })
    )
    
    # Reshape data for plotting
    clone_proportions_long <- clone_proportions %>%
      pivot_longer(cols = starts_with("clone_proportion"), names_to = "timepoint", values_to = "clone_proportion") %>%
      mutate(timepoint = recode(timepoint, "clone_proportion_C1" = "C1", "clone_proportion_C2" = "C2"))
    
    # Create the line plot
    line_plot <- ggplot(clone_proportions_long, aes(x = timepoint, y = clone_proportion, group = cdr3, color = cdr3)) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(y = "Clone Proportion", x = "Timepoint", title = paste("Patient:", patient)) +
      theme(legend.position = "none")
    
    # Add the plot to the list
    plot_list[[patient]] <- line_plot
  }
  
  # Combine all plots into a single PDF
  output_file <- file.path(output_folder, paste0("top_", top_n, "_clones_plots.pdf"))
  pdf(output_file, width = 7*length(plot_list), height = 7)
  print(plot_grid(plotlist = plot_list, ncol = length(plot_list), nrow = 1))
  dev.off()
}

# Example usage
plot_top_n_clones(seurat_object = seurat_object_tcr_cells, top_n = 5, output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/clonal_replacement_plots")





library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

plot_top_n_clones <- function(seurat_object_tcr_cells_metadata, top_n, output_folder, patient_metadata_df) {
  # Check for duplicate column names and make them unique
  colnames(seurat_object_tcr_cells_metadata) <- make.names(colnames(seurat_object_tcr_cells_metadata), unique = TRUE)
  # Subset patient metadata
  patient_metadata_df <- patient_metadata_df %>%
    filter(IDH != "POS" & site == "UF")
  
  # Get unique arms
  arms <- unique(patient_metadata_df$Arm)
  
  # Iterate over each arm
  for (arm in arms) {
    # Subset patient metadata for the current arm
    arm_metadata <- patient_metadata_df %>% filter(Arm == arm)
    
    # Order patients by decreasing OS.months
    patients_ordered <- arm_metadata %>% arrange(desc(OS.months.)) %>% pull(patient_id)
    
    # Create a list to store individual plots
    plot_list <- list()
    
    # Iterate over each patient in the ordered list and create plots
    for (patient in patients_ordered) {
      # Define origins for C1 and C2
      origin_C2 <- paste0(patient, "C2")
      origin_C1 <- paste0(patient, "C1")
      
      # Subset data for the current patient and origin C2
      patient_data_C2 <- seurat_object_tcr_cells_metadata %>% filter(Patient == patient & origin == origin_C2)
      
      # Identify top_n clones with the highest clone_proportion in origin C2
      top_clones <- patient_data_C2 %>%
        group_by(cdr3) %>%
        summarize(clone_proportion = max(clone_proportion)) %>%
        top_n(n = top_n, wt = clone_proportion) %>%
        pull(cdr3)
      
      # Subset data for the current patient and origin C1
      patient_data_C1 <- seurat_object_tcr_cells_metadata %>% filter(Patient == patient & origin == origin_C1)
      
      # Create a data frame to store clone proportions for C1 and C2
      clone_proportions <- data.frame(
        cdr3 = top_clones,
        clone_proportion_C1 = sapply(top_clones, function(clone) {
          prop <- patient_data_C1 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
          if (length(prop) == 0) 0 else prop[1]
        }),
        clone_proportion_C2 = sapply(top_clones, function(clone) {
          prop <- patient_data_C2 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
          if (length(prop) == 0) 0 else prop[1]
        })
      )
      
      # Reshape data for plotting
      clone_proportions_long <- clone_proportions %>%
        pivot_longer(cols = starts_with("clone_proportion"), names_to = "timepoint", values_to = "clone_proportion") %>%
        mutate(timepoint = recode(timepoint, "clone_proportion_C1" = "C1", "clone_proportion_C2" = "C2"))
      
      # Create the line plot
      line_plot <- ggplot(clone_proportions_long, aes(x = timepoint, y = clone_proportion, group = cdr3, color = cdr3)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(y = "Clone Proportion", x = "Timepoint", title = paste("Patient:", patient)) +
        theme(legend.position = "none") +
        scale_y_continuous(limits = c(0, 0.4))
      
      # Add the plot to the list
      plot_list[[patient]] <- line_plot
    }
    
    # Combine all plots into a single PDF for the current arm
    output_file <- file.path(output_folder, paste0("top_", top_n, "_clones_plots_", arm, ".pdf"))
    pdf(output_file, width = 7*length(plot_list), height = 7)
    print(plot_grid(plotlist = plot_list, ncol = length(plot_list), nrow = 1))
    dev.off()
  }
}


patient_metadata_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
patient_metadata_df <- read.csv(patient_metadata_file)
plot_top_n_clones(seurat_object = seurat_object_tcr_cells, top_n = 5, output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/clonal_replacement_plots", patient_metadata_df)




plot_top_n_clones_aggregated <- function(seurat_object_tcr_cells_metadata, top_n, output_folder, patient_metadata_df) {
  # Check for duplicate column names and make them unique
  colnames(seurat_object_tcr_cells_metadata) <- make.names(colnames(seurat_object_tcr_cells_metadata), unique = TRUE)
  # Subset patient metadata
  patient_metadata_df <- patient_metadata_df %>%
    filter(IDH != "POS" & site == "UF")
  
  # Get unique arms
  arms <- unique(patient_metadata_df$Arm)
  
  # Create a data frame to store aggregated data
  aggregated_data <- data.frame(Patient = character(), Arm = character(), Timepoint = character(), Mean_Clone_Proportion = numeric(), stringsAsFactors = FALSE)
  
  # Iterate over each arm
  for (arm in arms) {
    # Subset patient metadata for the current arm
    arm_metadata <- patient_metadata_df %>% filter(Arm == arm)
    
    # Order patients by decreasing OS.months
    patients_ordered <- arm_metadata %>% arrange(desc(OS.months.)) %>% pull(patient_id)
    
    # Iterate over each patient in the ordered list
    for (patient in patients_ordered) {
      # Define origins for C1 and C2
      origin_C2 <- paste0(patient, "C2")
      origin_C1 <- paste0(patient, "C1")
      
      # Subset data for the current patient and origin C2
      patient_data_C2 <- seurat_object_tcr_cells_metadata %>% filter(Patient == patient & origin == origin_C2)
      
      # Identify top_n clones with the highest clone_proportion in origin C2
      top_clones <- patient_data_C2 %>%
        group_by(cdr3) %>%
        summarize(clone_proportion = max(clone_proportion)) %>%
        top_n(n = top_n, wt = clone_proportion) %>%
        pull(cdr3)
      
      # Subset data for the current patient and origin C1
      patient_data_C1 <- seurat_object_tcr_cells_metadata %>% filter(Patient == patient & origin == origin_C1)
      
      # Create a data frame to store clone proportions for C1 and C2
      clone_proportions <- data.frame(
        cdr3 = top_clones,
        clone_proportion_C1 = sapply(top_clones, function(clone) {
          prop <- patient_data_C1 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
          if (length(prop) == 0) 0 else prop[1]
        }),
        clone_proportion_C2 = sapply(top_clones, function(clone) {
          prop <- patient_data_C2 %>% filter(cdr3 == clone) %>% pull(clone_proportion)
          if (length(prop) == 0) 0 else prop[1]
        })
      )
      
      # Calculate mean clone proportions for C1 and C2
      mean_proportion_C1 <- mean(clone_proportions$clone_proportion_C1, na.rm = TRUE)
      mean_proportion_C2 <- mean(clone_proportions$clone_proportion_C2, na.rm = TRUE)
      
      # Add data to the aggregated data frame
      aggregated_data <- rbind(aggregated_data, data.frame(Patient = patient, Arm = arm, Timepoint = "C1", Mean_Clone_Proportion = mean_proportion_C1))
      aggregated_data <- rbind(aggregated_data, data.frame(Patient = patient, Arm = arm, Timepoint = "C2", Mean_Clone_Proportion = mean_proportion_C2))
    }
  }
  
  print(aggregated_data)
  aggregated_data$Arm <- factor(aggregated_data$Arm, levels = c("MK-3475 Alone", "MK-3475 + MLA"))
  # Define custom colors
  custom_colors <- c("MK-3475 + MLA" = "red", "MK-3475 Alone" = "blue")
  
  # Create the aggregated line plot
  aggregated_plot <- ggplot(aggregated_data, aes(x = Timepoint, y = Mean_Clone_Proportion, group = Patient, color = Arm)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    scale_color_manual(values = custom_colors) +
    labs(y = "Mean Clone Proportion", x = "Timepoint", title = "Aggregated Clone Proportion by Timepoint and Arm")
    # scale_y_continuous(limits = c(0, 0.4))
  
  # Save the plot to a PDF
  output_file <- file.path(output_folder, paste0("aggregated_top_", top_n, "_clones_plots_17_18_shared.pdf"))
  ggsave(output_file, plot = aggregated_plot, width = 14, height = 7)
  
  # Print the plot
  print(aggregated_plot)
  return(aggregated_data)
}

# # Example usage
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells@meta.data
# 
# #########################################
# # only shared clones
# seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/MK_T_Cells_reclustered_res_1_seurat_object_GBM.rds")
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
# 
# ##########################################################################################################################################################
# # subset seurat object for cell barcodes where barcodes match the TCR barcodes and tag cell barcodes as "shared" or "not shared"
# # map TCR barcodes
# blood_clone_cell_barcode_df <- read.csv(blood_clone_cell_barcode_file)
# cell_barcodes <- unique(blood_clone_cell_barcode_df$barcode)
# 
# seurat_object_tcr_barcodes <- subset(seurat_object_t_cells, cells = cell_barcodes)
# 
# # add an attribute to the seurat object metadata, "clone_status"
# # Initialize the clone_status column in metadata
# seurat_object_tcr_barcodes@meta.data$clone_status <- "not shared"
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
#   seurat_object_tcr_barcodes@meta.data[shared_cells, "clone_status"] <- "shared"
# }

#########################################
seurat_object_shared_non_shared_metadata <- seurat_object_tcr_barcodes@meta.data
seurat_object_shared_non_shared_metadata <- seurat_object_shared_non_shared_metadata[seurat_object_shared_non_shared_metadata$clone_status == "shared",]
seurat_object_shared_non_shared_metadata$barcode <- rownames(seurat_object_shared_non_shared_metadata)
seurat_object_shared_non_shared_metadata <- merge(seurat_object_shared_non_shared_metadata[, c("barcode", "clone_status")], seurat_object_tcr_cells_metadata, by.x = "barcode", by.y="row.names")
seurat_object_shared_non_shared_metadata <- seurat_object_shared_non_shared_metadata[seurat_object_shared_non_shared_metadata$seurat_clusters %in% c(17, 18), ]
# seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$seurat_clusters %in% c(18), ]
temp <- plot_top_n_clones_aggregated(seurat_object_shared_non_shared_metadata, top_n = 10, output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/Finale/GBM/clonal_expansion/clonal_replacement_plots", patient_metadata_df)
