library(Seurat)
library(ggplot2)
library(dplyr)

plotMarkersSeurat <- function(markers_csv, output_folder, seurat_obj) {
  print(output_folder)
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  p <- DimPlot(seurat_obj, label = TRUE, raster = FALSE)
  ggsave(file.path(output_folder, "raw.png"), plot = p, device = "png", width = 10, height = 8)
  # Read the markers.csv file
  markers <- read.csv(markers_csv, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  markers$'Cell type' <- trimws(markers$'Cell type')
  
  # Loop over unique cell types in the markers file
  for (cell_type in unique(markers$'Cell type')) {
    
    # Enclose the plotting and saving in a tryCatch block
    tryCatch({
      print(paste("Plotting for cell type:", cell_type))
      
      # Get all marker genes for the cell type
      marker_genes <- markers$Name[markers$'Cell type' == cell_type]
      
      # Plot FeaturePlot for each cell type
      p <- FeaturePlot(seurat_obj, features = marker_genes, label = TRUE, cols = c("lightblue","red"), raster = FALSE)
      
      # Save the plot as a PNG file with the cell type name in the output directory
      png_filename <- file.path(output_folder, paste0(gsub(" ", "_", cell_type), ".png"))
      ggsave(png_filename, plot = p, device = "png", width = 20, height = 16)
      
    }, error = function(e) {
      # Print an error message and skip to the next cell type
      print(paste("Error plotting for cell type:", cell_type, "Error:", e$message))
    })
  }
}

markers_csv = "/project/dtran642_927/SonLe/USC_Source/source/Single_Cell/markers.csv"

# choosing the resolution

resolution_0.1 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.1/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.1/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_0.1$seurat_obj)

resolution_0.3 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.3/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.3/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_0.3$seurat_obj)

resolution_1 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_1$seurat_obj)

resolution_3 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_3/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_3/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_3$seurat_obj)



# chosen resolution 1, now taking out T cell clusters to recluster
# T Cell clusters - (5,8,16,23,24,25,27,28,32,46)

MK_T_Cells_seurat_obj <- subset(resolution_1$seurat_obj, subset = seurat_clusters %in% c(5,8,16,23,24,25,27,28,32,46))

# reclustering
# 2. Identify highly variable features (genes)
MK_T_Cells_seurat_obj <- FindVariableFeatures(MK_T_Cells_seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_T_Cells_seurat_obj)
MK_T_Cells_seurat_obj <- ScaleData(MK_T_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_T_Cells_seurat_obj <- RunPCA(MK_T_Cells_seurat_obj, features = VariableFeatures(object = MK_T_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_T_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_T_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_T_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_T_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:10  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_T_Cells_seurat_obj <- FindNeighbors(MK_T_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_T_Cells_seurat_obj <- FindClusters(MK_T_Cells_seurat_obj, resolution = 1)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_T_Cells_seurat_obj <- RunUMAP(MK_T_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

saveRDS(MK_T_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# creating plots for t cell subpopulation
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_sub_population_annotation_helper_plots/res_1"
plotMarkersSeurat(markers_csv, output_folder, MK_T_Cells_seurat_obj)

# Extract UMAP coordinates
umap_coords <- Embeddings(MK_T_Cells_seurat_obj, "umap")
# Extract cluster annotations
cluster_annotations <- Idents(MK_T_Cells_seurat_obj)
# Combine data into a data frame
export_data <- data.frame(Barcode = rownames(umap_coords), umap_coords, Cluster = as.integer(cluster_annotations))
# Write to CSV
write.csv(export_data, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_sub_population_annotation_helper_plots/res_1/loupe.csv", row.names = FALSE, quote = FALSE)


# preparing T cell seurat object for downstream analysis
MK_T_Cells_seurat_obj = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# removing irrelevant clusters from T cells
irrelevant_clusters <- c(11,15,19,20)
relevant_clusters <- unique(MK_T_Cells_seurat_obj@meta.data$seurat_clusters)[!unique(MK_T_Cells_seurat_obj@meta.data$seurat_clusters) %in% irrelevant_clusters]
MK_T_Cells_seurat_obj <- subset(MK_T_Cells_seurat_obj, subset = seurat_clusters %in% relevant_clusters)
saveRDS(MK_T_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# removing irrelevant clusters from All cells
irrelevant_clusters <- c(19,20,27, 15, 17, 18, 21, 41, 42)
relevant_clusters <- unique(resolution_1$seurat_obj@meta.data$seurat_clusters)[!unique(resolution_1$seurat_obj@meta.data$seurat_clusters) %in% irrelevant_clusters]
MK_Cells_seurat_obj <- subset(resolution_1$seurat_obj, subset = seurat_clusters %in% relevant_clusters)
saveRDS(MK_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")

#########################################################################################################################################################################################
# preparing Monocyte seurat object
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
MK_Monocyte_Cells_seurat_obj <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(1, 3, 6, 12, 9, 33))

# reclustering
# 2. Identify highly variable features (genes)
MK_Monocyte_Cells_seurat_obj <- FindVariableFeatures(seurat_object_monocyte_cells, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_Monocyte_Cells_seurat_obj)
MK_Monocyte_Cells_seurat_obj <- ScaleData(MK_Monocyte_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_Monocyte_Cells_seurat_obj <- RunPCA(MK_Monocyte_Cells_seurat_obj, features = VariableFeatures(object = MK_Monocyte_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_Monocyte_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_Monocyte_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_Monocyte_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_Monocyte_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:15  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_Monocyte_Cells_seurat_obj <- FindNeighbors(MK_Monocyte_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_Monocyte_Cells_seurat_obj <- FindClusters(MK_Monocyte_Cells_seurat_obj, resolution = 0.3)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_Monocyte_Cells_seurat_obj <- RunUMAP(MK_Monocyte_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:https://ondemand.carc.usc.edu/rnode/a01-14.hpc.usc.edu/26879/graphics/2af3ac3b-a1c6-4c4a-9f0a-de7d0f268b04.png
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

saveRDS(MK_Monocyte_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Monocyte_Cells_seurat_obj.RDS")




#########################################################################################################################################################################################
# preparing Monocyte seurat object
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
MK_NC_Monocyte_Cells_seurat_obj <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(9, 33))

# reclustering
# 2. Identify highly variable features (genes)
MK_NC_Monocyte_Cells_seurat_obj <- FindVariableFeatures(MK_NC_Monocyte_Cells_seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_NC_Monocyte_Cells_seurat_obj)
MK_NC_Monocyte_Cells_seurat_obj <- ScaleData(MK_NC_Monocyte_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_NC_Monocyte_Cells_seurat_obj <- RunPCA(MK_NC_Monocyte_Cells_seurat_obj, features = VariableFeatures(object = MK_NC_Monocyte_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_NC_Monocyte_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_NC_Monocyte_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_NC_Monocyte_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_NC_Monocyte_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:20  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_NC_Monocyte_Cells_seurat_obj <- FindNeighbors(MK_NC_Monocyte_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_NC_Monocyte_Cells_seurat_obj <- FindClusters(MK_NC_Monocyte_Cells_seurat_obj, resolution = 1)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_NC_Monocyte_Cells_seurat_obj <- RunUMAP(MK_NC_Monocyte_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:https://ondemand.carc.usc.edu/rnode/a01-14.hpc.usc.edu/26879/graphics/2af3ac3b-a1c6-4c4a-9f0a-de7d0f268b04.png
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

saveRDS(MK_NC_Monocyte_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")






#############################################################################################
FeaturePlot(resolution_1$seurat_obj, features = c("SV40-large-T-Antigene"), raster = FALSE, split.by = "Patient", label = TRUE, cols = c("lightblue", "red"))


######################################################################################################################################################
# creating h5ad object of the Non-Classical Monocytes for experiment patients
library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")

# Extract UMAP embeddings
umap_embeddings <- Embeddings(seurat_obj, reduction = "umap")

# Identify cells within the specified UMAP axis limits
cells_to_keep <- rownames(umap_embeddings)[
  umap_embeddings[, "UMAP_1"] >= -7 & umap_embeddings[, "UMAP_1"] <= 7 &
    umap_embeddings[, "UMAP_2"] >= -6 & umap_embeddings[, "UMAP_2"] <= 8
]

# Subset the Seurat object to keep only the cells within the limits
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

# remove unwanted clusters
seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16) & Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21) & TimePoint %in% c("Pre", "C1"))
# Set the default assay to 'RNA' to include all genes
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
# Save as h5Seurat
SaveH5Seurat(seurat_obj, filename = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.h5Seurat")
# Convert to h5ad
Convert("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.h5Seurat", dest = "h5ad")

# Read the mapping CSV
mapping_df <- read.csv('/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/predecessor_mapping.csv', stringsAsFactors = FALSE)

# Inspect the first few rows
head(mapping_df)

# Extract all cell barcodes from the Seurat object
all_barcodes <- colnames(seurat_obj)

# Initialize the new metadata column with NA
seurat_obj$predecessor_barcode <- NA

# Ensure that 'target_barcode' in mapping_df corresponds to cells in the Seurat object
# You may need to adjust based on your specific barcode naming conventions

# Create a named vector for mapping
predecessor_vector <- setNames(mapping_df$predecessor_barcode, mapping_df$target_barcode)

# Identify cells in Timepoint 1
# Replace 'timepoint' with the actual metadata column name for timepoints in your Seurat object
# and ensure that Timepoint 1 is labeled as '1' or as per your data
time1_cells <- WhichCells(seurat_obj, expression = TimePoint == "C1")

# Assign predecessors to Timepoint 1 cells
# Only update cells that are present in the mapping
common_cells <- intersect(time1_cells, mapping_df$target_barcode)
seurat_obj$predecessor_barcode[common_cells] <- predecessor_vector[common_cells]

# Save updated Seurat object
saveRDS(seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")

######################################################################################################################################################
# creating h5ad object of the Non-Classical Monocytes for control patients
library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")

# Extract UMAP embeddings
umap_embeddings <- Embeddings(seurat_obj, reduction = "umap")

# Identify cells within the specified UMAP axis limits
cells_to_keep <- rownames(umap_embeddings)[
  umap_embeddings[, "UMAP_1"] >= -7 & umap_embeddings[, "UMAP_1"] <= 7 &
    umap_embeddings[, "UMAP_2"] >= -6 & umap_embeddings[, "UMAP_2"] <= 8
]

# Subset the Seurat object to keep only the cells within the limits
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

# remove unwanted clusters
seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16) & Patient %in% c(1, 4, 8, 9, 11) & TimePoint %in% c("Pre", "C1"))
# Set the default assay to 'RNA' to include all genes
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
# Save as h5Seurat
SaveH5Seurat(seurat_obj, filename = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_Pre_C1.h5Seurat")
# Convert to h5ad
Convert("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_Pre_C1.h5Seurat", dest = "h5ad")

# Read the mapping CSV
mapping_df <- read.csv('/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/predecessor_mapping.csv', stringsAsFactors = FALSE)

# Inspect the first few rows
head(mapping_df)

# Extract all cell barcodes from the Seurat object
all_barcodes <- colnames(seurat_obj)

# Initialize the new metadata column with NA
seurat_obj$predecessor_barcode <- NA

# Ensure that 'target_barcode' in mapping_df corresponds to cells in the Seurat object
# You may need to adjust based on your specific barcode naming conventions

# Create a named vector for mapping
predecessor_vector <- setNames(mapping_df$predecessor_barcode, mapping_df$target_barcode)

# Identify cells in Timepoint 1
# Replace 'timepoint' with the actual metadata column name for timepoints in your Seurat object
# and ensure that Timepoint 1 is labeled as '1' or as per your data
time1_cells <- WhichCells(seurat_obj, expression = TimePoint == "C1")

# Assign predecessors to Timepoint 1 cells
# Only update cells that are present in the mapping
common_cells <- intersect(time1_cells, mapping_df$target_barcode)
seurat_obj$predecessor_barcode[common_cells] <- predecessor_vector[common_cells]

# Save updated Seurat object
saveRDS(seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_Pre_C1.RDS")



###########################################################################################################
# DimPlot All Cells for publication
MK_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
DimPlot(MK_Cells_seurat_obj, raster = FALSE, label = FALSE)
# Define the new cluster mappings
merged_clusters <- MK_Cells_seurat_obj@meta.data$seurat_clusters
# Convert to character for easier manipulation
merged_clusters <- as.character(merged_clusters)
# Define the mappings
merged_clusters[merged_clusters %in% c("1", "3", "6", "12")] <- "1"
merged_clusters[merged_clusters %in% c("9", "33")] <- "9"
MK_Cells_seurat_obj@meta.data$seurat_clusters <- merged_clusters
DimPlot(MK_Cells_seurat_obj, raster = FALSE, label = FALSE, group.by = "seurat_clusters")

dim_plot <- DimPlot(MK_Cells_seurat_obj, raster = FALSE, label = FALSE, group.by = "seurat_clusters") + 
  ggtitle("All Cells UMAP")

# Save the plot to a PDF file
ggsave(
  filename = "UMAP_all_cells.pdf",       # Name of the output file
  plot = dim_plot,                             # Plot to save
  device = "pdf",                              # File format
  path = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material",             # Optional: Specify path
  width = 16, height = 12,                        # Dimensions in inches
  units = "in",                                # Units for width and height
  dpi = 300                                    # Resolution (optional for PDF)
)


###########################################################################################################
# DimPlot T Cells for publication
MK_T_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
MK_T_Cells_seurat_obj <- subset(MK_T_Cells_seurat_obj, subset = seurat_clusters != 13)
# Define the new cluster mappings
merged_clusters <- MK_T_Cells_seurat_obj@meta.data$seurat_clusters
# Convert to character for easier manipulation
merged_clusters <- as.character(merged_clusters)
# Define the mappings
merged_clusters[merged_clusters %in% c("0", "5")] <- "0"
merged_clusters[merged_clusters %in% c("6", "9", "18")] <- "6"
merged_clusters[merged_clusters %in% c("14", "16", "17")] <- "14"
MK_T_Cells_seurat_obj@meta.data$seurat_clusters <- merged_clusters
DimPlot(MK_T_Cells_seurat_obj, raster = FALSE, label = FALSE, group.by = "seurat_clusters") + coord_flip()

# Ensure 'seurat_clusters' is a factor with levels ordered as desired
MK_T_Cells_seurat_obj@meta.data$seurat_clusters <- factor(MK_T_Cells_seurat_obj@meta.data$seurat_clusters, 
                                                          levels = sort(unique(merged_clusters)))

dim_plot <- DimPlot(MK_T_Cells_seurat_obj, raster = FALSE, label = FALSE, group.by = "seurat_clusters") + coord_flip() +
  ggtitle("T Cells UMAP")

# Build the ggplot object to access its data
plot_build <- ggplot_build(dim_plot)

# Inspect the scales to identify the color scale
print(dim_plot$scales$scales)  # Optional: to understand the scales present

# Correctly identify the color scale by checking if "colour" or "color" is in aesthetics
color_scale_index <- which(sapply(dim_plot$scales$scales, function(x) {
  "colour" %in% x$aesthetics || "color" %in% x$aesthetics
}))

# Check if a color scale was found
if(length(color_scale_index) == 0){
  stop("No color scale found in the plot.")
}

# Extract the color scale
color_scale <- dim_plot$scales$scales[[color_scale_index[1]]]

# Depending on the type of scale, extract the colors differently
# For discrete scales like scale_color_manual or scale_color_hue
if("scale" %in% class(color_scale)){
  # Extract the mapping from the scale
  cluster_labels <- color_scale$palette(length(color_scale$range$range))
  cluster_colors <- color_scale$range$values
  
  # If 'palette' is a function, generate colors based on the number of clusters
  if(is.function(color_scale$palette)){
    num_clusters <- length(levels(MK_T_Cells_seurat_obj@meta.data$seurat_clusters))
    cluster_colors <- color_scale$palette(num_clusters)
    cluster_labels <- levels(MK_T_Cells_seurat_obj@meta.data$seurat_clusters)
  } else {
    cluster_labels <- color_scale$levels
  }
  
} else {
  # For other scale types, you might need a different extraction method
  stop("Unsupported scale type for color extraction.")
}

# Create a named vector mapping clusters to colors
cluster_color_mapping <- setNames(cluster_colors, cluster_labels)

# Display the cluster-color mappings
print(cluster_color_mapping)

# (Optional) Save the color mappings to a CSV file
cluster_colors_df <- data.frame(
  Cluster = cluster_labels, 
  Color = cluster_colors, 
  stringsAsFactors = FALSE
)

# (Optional) Save the color mappings to a CSV file
write.csv(cluster_colors_df, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material/T_Cell_cluster_colors.csv", row.names = FALSE)

# Save the plot to a PDF file
ggsave(
  filename = "UMAP_T_cells.pdf",       # Name of the output file
  plot = dim_plot,                             # Plot to save
  device = "pdf",                              # File format
  path = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material",             # Optional: Specify path
  width = 8, height = 6,                        # Dimensions in inches
  units = "in",                                # Units for width and height
  dpi = 300                                    # Resolution (optional for PDF)
)




###########################################################################################################
# DimPlot NC Mono Cells for publication
MK_NC_Mono_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")
patient_specific_clusters <- c(7, 8, 2)
relevant_clusters <- setdiff(MK_NC_Mono_Cells_seurat_obj@meta.data$seurat_clusters, patient_specific_clusters)
MK_NC_Mono_Cells_seurat_obj <- subset(MK_NC_Mono_Cells_seurat_obj, subset = seurat_clusters %in% relevant_clusters)

dim_plot <- DimPlot(
  MK_NC_Mono_Cells_seurat_obj, 
  raster = FALSE, 
  label = FALSE, 
  group.by = "seurat_clusters"
) +
  scale_x_continuous(limits = c(-7, 7)) +  # Set x-axis limits for UMAP_1
  scale_y_continuous(limits = c(-6, 8)) +  # Set y-axis limits for UMAP_2
  ggtitle("NC Monocyte Cells UMAP")

# Save the plot to a PDF file
ggsave(
  filename = "UMAP_NC_Mono_cells.pdf",       
  plot = dim_plot,                             
  device = "pdf",                              
  path = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material",             
  width = 8, height = 6,                        
  units = "in",                                
  dpi = 300                                    
)




######################################################################################################################################################
# creating h5ad object of the Non-Classical Monocytes for control patients (1, 4)
library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")

# Extract UMAP embeddings
umap_embeddings <- Embeddings(seurat_obj, reduction = "umap")

# Identify cells within the specified UMAP axis limits
cells_to_keep <- rownames(umap_embeddings)[
  umap_embeddings[, "UMAP_1"] >= -7 & umap_embeddings[, "UMAP_1"] <= 7 &
    umap_embeddings[, "UMAP_2"] >= -6 & umap_embeddings[, "UMAP_2"] <= 8
]

# Subset the Seurat object to keep only the cells within the limits
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

# remove unwanted clusters
seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16) & Patient %in% c(1, 4) & TimePoint %in% c("Pre", "C1"))
# Set the default assay to 'RNA' to include all genes
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
# Save as h5Seurat
SaveH5Seurat(seurat_obj, filename = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_patient_1_4_Pre_C1.h5Seurat")
# Convert to h5ad
Convert("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_patient_1_4_Pre_C1.h5Seurat", dest = "h5ad")

# Read the mapping CSV
mapping_df <- read.csv('/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/predecessor_mapping.csv', stringsAsFactors = FALSE)

# Inspect the first few rows
head(mapping_df)

# Extract all cell barcodes from the Seurat object
all_barcodes <- colnames(seurat_obj)

# Initialize the new metadata column with NA
seurat_obj$predecessor_barcode <- NA

# Ensure that 'target_barcode' in mapping_df corresponds to cells in the Seurat object
# You may need to adjust based on your specific barcode naming conventions

# Create a named vector for mapping
predecessor_vector <- setNames(mapping_df$predecessor_barcode, mapping_df$target_barcode)

# Identify cells in Timepoint 1
# Replace 'timepoint' with the actual metadata column name for timepoints in your Seurat object
# and ensure that Timepoint 1 is labeled as '1' or as per your data
time1_cells <- WhichCells(seurat_obj, expression = TimePoint == "C1")

# Assign predecessors to Timepoint 1 cells
# Only update cells that are present in the mapping
common_cells <- intersect(time1_cells, mapping_df$target_barcode)
seurat_obj$predecessor_barcode[common_cells] <- predecessor_vector[common_cells]

# Save updated Seurat object
saveRDS(seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_patient_1_4_Pre_C1.RDS")






#########################################################################################################################################################################################
# preparing DC seurat object
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
MK_DC_Cells_seurat_obj <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(31, 36))

# reclustering
# 2. Identify highly variable features (genes)
MK_DC_Cells_seurat_obj <- FindVariableFeatures(MK_DC_Cells_seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_DC_Cells_seurat_obj)
MK_DC_Cells_seurat_obj <- ScaleData(MK_DC_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_DC_Cells_seurat_obj <- RunPCA(MK_DC_Cells_seurat_obj, features = VariableFeatures(object = MK_DC_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_DC_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_DC_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_DC_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_DC_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:20  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_DC_Cells_seurat_obj <- FindNeighbors(MK_DC_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_DC_Cells_seurat_obj <- FindClusters(MK_DC_Cells_seurat_obj, resolution = 1)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_DC_Cells_seurat_obj <- RunUMAP(MK_DC_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:https://ondemand.carc.usc.edu/rnode/a01-14.hpc.usc.edu/26879/graphics/2af3ac3b-a1c6-4c4a-9f0a-de7d0f268b04.png
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

# Subset for experiment patients and timepoints Pre an C1
MK_DC_Cells_seurat_obj <- subset(MK_DC_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21) & TimePoint %in% c("Pre", "C1"))


saveRDS(MK_DC_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_DC_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")



#########################################################################################################################################################################################
# prepare DC Seurat object and keep only CCR2-positive cells
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")

DefaultAssay(seurat_object_all_cells) <- "RNA"

MK_DC_Cells_seurat_obj <- subset(seurat_object_all_cells,
                                 subset = seurat_clusters %in% c(31, 36) & CCR2 > 0)
# reclustering
# 2. Identify highly variable features (genes)
MK_DC_Cells_seurat_obj <- FindVariableFeatures(MK_DC_Cells_seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_DC_Cells_seurat_obj)
MK_DC_Cells_seurat_obj <- ScaleData(MK_DC_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_DC_Cells_seurat_obj <- RunPCA(MK_DC_Cells_seurat_obj, features = VariableFeatures(object = MK_DC_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_DC_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_DC_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_DC_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_DC_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:20  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_DC_Cells_seurat_obj <- FindNeighbors(MK_DC_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_DC_Cells_seurat_obj <- FindClusters(MK_DC_Cells_seurat_obj, resolution = 1)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_DC_Cells_seurat_obj <- RunUMAP(MK_DC_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:https://ondemand.carc.usc.edu/rnode/a01-14.hpc.usc.edu/26879/graphics/2af3ac3b-a1c6-4c4a-9f0a-de7d0f268b04.png
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

# Subset for experiment patients and timepoints Pre an C1
MK_DC_Cells_seurat_obj <- subset(MK_DC_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21) & TimePoint %in% c("Pre", "C1"))


saveRDS(MK_DC_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_DC_CCR2_Pos_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")



#########################################################################################################################################################################################
# preparing Classical Monocyte seurat object
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
MK_C_Mono_Cells_seurat_obj <- subset(seurat_object_all_cells, subset = seurat_clusters %in% c(1,3,6,12))

# reclustering
# 2. Identify highly variable features (genes)
MK_C_Mono_Cells_seurat_obj <- FindVariableFeatures(MK_C_Mono_Cells_seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. Scale the data
all.genes <- rownames(MK_C_Mono_Cells_seurat_obj)
MK_C_Mono_Cells_seurat_obj <- ScaleData(MK_C_Mono_Cells_seurat_obj, features = all.genes)

# 4. Perform linear dimensional reduction (PCA)
MK_C_Mono_Cells_seurat_obj <- RunPCA(MK_C_Mono_Cells_seurat_obj, features = VariableFeatures(object = MK_C_Mono_Cells_seurat_obj))

# Optional: Visualize PCA results
print(MK_C_Mono_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_C_Mono_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_C_Mono_Cells_seurat_obj, reduction = "pca")
ElbowPlot(MK_C_Mono_Cells_seurat_obj)

# 5. Determine the dimensionality to use
# Based on ElbowPlot or other criteria, decide on the number of PCs
pcs_to_use <- 1:20  # Adjust based on your data

# 6. Construct a nearest neighbor graph
MK_C_Mono_Cells_seurat_obj <- FindNeighbors(MK_C_Mono_Cells_seurat_obj, dims = pcs_to_use)

# 7. Cluster the cells
MK_C_Mono_Cells_seurat_obj <- FindClusters(MK_C_Mono_Cells_seurat_obj, resolution = 1)  # Adjust resolution as needed

# 8. Run non-linear dimensional reduction (UMAP/t-SNE)
MK_C_Mono_Cells_seurat_obj <- RunUMAP(MK_C_Mono_Cells_seurat_obj, dims = pcs_to_use)
# Or, for t-SNE:https://ondemand.carc.usc.edu/rnode/a01-14.hpc.usc.edu/26879/graphics/2af3ac3b-a1c6-4c4a-9f0a-de7d0f268b04.png
# t_cells <- RunTSNE(t_cells, dims = pcs_to_use)

# Subset for experiment patients and timepoints Pre an C1
MK_C_Mono_Cells_seurat_obj <- subset(MK_C_Mono_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21) & TimePoint %in% c("Pre", "C1"))


saveRDS(MK_C_Mono_Cells_seurat_obj, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_C_Mono_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")

