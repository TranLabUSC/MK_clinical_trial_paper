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


MK_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")
# Extract UMAP coordinates
umap_coords <- Embeddings(MK_Cells_seurat_obj, "umap")
# Extract cluster annotations
cluster_annotations <- Idents(MK_Cells_seurat_obj)
# Combine data into a data frame
export_data <- data.frame(Barcode = rownames(umap_coords), umap_coords, Cluster = as.integer(cluster_annotations))
# Write to CSV
write.csv(export_data, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/all_cells_loupe.csv", row.names = FALSE, quote = FALSE)


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
