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



#######################################################################################################################################################################################################################
# read seurat object
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
# only considering UF data as WUSTL TCR data is not good quality
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF" & Patient != 1)
seurat_metadata <- seurat_object_t_cells@meta.data
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
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

seurat_object_tcr_cells@meta.data$clone_proportion_capped <- pmin(seurat_object_tcr_cells@meta.data$clone_proportion, 0.2)
# Generate UMAP plot colored by the 'frequency' column
umap_plot <- FeaturePlot(seurat_object_tcr_cells, features = "clone_proportion_capped", reduction = "umap") +
  scale_color_gradientn(colors = c("blue", "yellow", "red"), 
                        values = scales::rescale(c(0, 0.05, 0.10, 0.15, 0.20)),
                        breaks = c(0, 0.05, 0.10, 0.15, 0.20))
umap_plot + coord_flip()


control_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("8", "11", "9", "4"), ])
exp_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$Patient %in% c("21", "13", "2", "3", "18", "19", "20", "5", "14", "12", "10", "7"), ])
sampled_exp_cells <- sample(exp_cells, length(control_cells))

umap_plot <- FeaturePlot(seurat_object_tcr_cells, features = "clone_proportion_capped", reduction = "umap", cells = sampled_exp_cells) +
  scale_color_gradientn(colors = c("blue", "yellow", "red"), 
                        values = scales::rescale(c(0, 0.05, 0.10, 0.15, 0.20)),
                        breaks = c(0, 0.05, 0.10, 0.15, 0.20))
umap_plot + coord_flip()

