################################################################################
# CELL TYPE ANNOTATION USING MARKER GENES
################################################################################
#
# PURPOSE:
#   This script plots marker gene expression across different clustering resolutions
#   to assist with cell type annotation. The marker gene expression plots are used
#   to manually annotate cell clusters based on known cell type-specific markers.
#
# MANUSCRIPT FIGURES:
#   - Figure 5a: All cells UMAP with cell type annotations
#   - Figure 6a: T cell subpopulation UMAP with cluster annotations
#
# WORKFLOW:
#   1. Define a function to plot marker genes for each cell type
#   2. Test multiple clustering resolutions (0.1, 0.3, 1, 3) to find optimal granularity
#   3. Extract and recluster T cells for detailed subpopulation analysis
#   4. Remove irrelevant/low-quality clusters
#   5. Export annotated data for visualization (Loupe browser)
#
# INPUT FILES:
#   - markers.csv: List of marker genes for each cell type
#   - Seurat objects at different clustering resolutions
#
# OUTPUT FILES:
#   - Feature plots showing marker gene expression for each cell type
#   - Annotated UMAP coordinates (CSV format for Loupe browser)
#   - Final cleaned Seurat objects (T cells and all cells)
#
################################################################################

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

################################################################################
# FUNCTION: plotMarkersSeurat
################################################################################
# Generates feature plots for marker genes grouped by cell type
#
# PARAMETERS:
#   markers_csv    - Path to CSV file containing marker genes and cell types
#   output_folder  - Directory where plots will be saved
#   seurat_obj     - Seurat object containing clustered single-cell data
#
# OUTPUT:
#   - raw.png: UMAP plot with cluster labels
#   - {Cell_Type}.png: Feature plots showing expression of markers for each cell type
#
################################################################################
plotMarkersSeurat <- function(markers_csv, output_folder, seurat_obj) {
  print(output_folder)
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Generate and save base UMAP plot with cluster labels
  p <- DimPlot(seurat_obj, label = TRUE, raster = FALSE)
  ggsave(file.path(output_folder, "raw.png"), plot = p, device = "png", width = 10, height = 8)
  
  # Read the markers CSV file
  # Format: columns should include 'Name' (gene name) and 'Cell type'
  markers <- read.csv(markers_csv, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  markers$'Cell type' <- trimws(markers$'Cell type')
  
  # Loop over unique cell types and generate feature plots
  for (cell_type in unique(markers$'Cell type')) {
    
    # Error handling: skip cell types that cause plotting errors
    tryCatch({
      print(paste("Plotting for cell type:", cell_type))
      
      # Get all marker genes for the current cell type
      marker_genes <- markers$Name[markers$'Cell type' == cell_type]
      
      # Plot FeaturePlot showing expression of all markers for this cell type
      # Color scale: lightblue (low expression) to red (high expression)
      p <- FeaturePlot(seurat_obj, features = marker_genes, label = TRUE, 
                      cols = c("lightblue","red"), raster = FALSE)
      
      # Save the plot with cell type name (spaces replaced with underscores)
      png_filename <- file.path(output_folder, paste0(gsub(" ", "_", cell_type), ".png"))
      ggsave(png_filename, plot = p, device = "png", width = 20, height = 16)
      
    }, error = function(e) {
      # Print error message but continue with next cell type
      print(paste("Error plotting for cell type:", cell_type, "Error:", e$message))
    })
  }
}

################################################################################
# SECTION 1: TEST MULTIPLE CLUSTERING RESOLUTIONS
################################################################################
# Marker genes file containing known cell type-specific genes
markers_csv = "/project/dtran642_927/SonLe/USC_Source/source/Single_Cell/markers.csv"

# Resolution 0.1 - Coarse clustering (fewer clusters)
resolution_0.1 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.1/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.1/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_0.1$seurat_obj)

# Resolution 0.3
resolution_0.3 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.3/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_0.3/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_0.3$seurat_obj)

# Resolution 1 - CHOSEN RESOLUTION for downstream analysis
resolution_1 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_1$seurat_obj)

# Resolution 3 - Fine clustering (more clusters)
resolution_3 = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_3/seurat_obj_after_clustering.RDS")
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_3/plots/t_cell_annotation_plots"
plotMarkersSeurat(markers_csv, output_folder, resolution_3$seurat_obj)


################################################################################
# SECTION 2: EXPORT ALL CELLS DATA FOR LOUPE BROWSER VISUALIZATION
################################################################################
# Load the full dataset at resolution 1
MK_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")

# Extract UMAP coordinates for visualization
umap_coords <- Embeddings(MK_Cells_seurat_obj, "umap")

# Extract cluster annotations
cluster_annotations <- Idents(MK_Cells_seurat_obj)

# Combine barcode, UMAP coordinates, and cluster ID into a data frame
export_data <- data.frame(Barcode = rownames(umap_coords), 
                         umap_coords, 
                         Cluster = as.integer(cluster_annotations))

# Export to CSV for loading into 10x Genomics Loupe Browser
write.csv(export_data, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/all_cells_loupe.csv", row.names = FALSE, quote = FALSE)


################################################################################
# SECTION 3: ISOLATE AND RECLUSTER T CELLS
################################################################################
# T Cell clusters identified from resolution 1: (5,8,16,23,24,25,27,28,32,46)
# These were identified based on marker gene expression (CD3D, CD3E, CD3G, etc.)

# Subset T cell clusters from the full dataset
MK_T_Cells_seurat_obj <- subset(resolution_1$seurat_obj, 
                                subset = seurat_clusters %in% c(5,8,16,23,24,25,27,28,32,46))

################################################################################
# STANDARD SEURAT RECLUSTERING WORKFLOW FOR T CELLS
################################################################################

# Step 1: Identify highly variable features (genes)
# Using variance stabilizing transformation (vst) method
# Selecting top 2000 most variable genes
MK_T_Cells_seurat_obj <- FindVariableFeatures(MK_T_Cells_seurat_obj, 
                                              selection.method = "vst", 
                                              nfeatures = 2000)

# Step 2: Scale the data
# Z-score normalization across all cells
all.genes <- rownames(MK_T_Cells_seurat_obj)
MK_T_Cells_seurat_obj <- ScaleData(MK_T_Cells_seurat_obj, features = all.genes)

# Step 3: Perform linear dimensional reduction (PCA)
# Using the variable features identified above
MK_T_Cells_seurat_obj <- RunPCA(MK_T_Cells_seurat_obj, 
                                features = VariableFeatures(object = MK_T_Cells_seurat_obj))

# Step 4: Visualize PCA results (QC step)
# Print top genes contributing to first 5 PCs
print(MK_T_Cells_seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MK_T_Cells_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(MK_T_Cells_seurat_obj, reduction = "pca")

# Elbow plot to determine optimal number of PCs to use
ElbowPlot(MK_T_Cells_seurat_obj)

# Step 5: Determine dimensionality
# Based on ElbowPlot, using first 10 PCs
pcs_to_use <- 1:10

# Step 6: Construct k-nearest neighbor (KNN) graph
# Uses PCA space to find nearest neighbors
MK_T_Cells_seurat_obj <- FindNeighbors(MK_T_Cells_seurat_obj, dims = pcs_to_use)

# Step 7: Cluster the cells
# Using Louvain algorithm at resolution 1
MK_T_Cells_seurat_obj <- FindClusters(MK_T_Cells_seurat_obj, resolution = 1)

# Step 8: Run UMAP for visualization
# Non-linear dimensional reduction to 2D
MK_T_Cells_seurat_obj <- RunUMAP(MK_T_Cells_seurat_obj, dims = pcs_to_use)

# Save the reclustered T cell object
saveRDS(MK_T_Cells_seurat_obj, 
       file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")


################################################################################
# SECTION 4: GENERATE MARKER PLOTS FOR T CELL SUBPOPULATIONS
################################################################################
# Plot marker genes to identify T cell subpopulations
# (e.g., CD4+ T cells, CD8+ T cells, naive, memory, effector, exhausted, etc.)
output_folder = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_sub_population_annotation_helper_plots/res_1"
plotMarkersSeurat(markers_csv, output_folder, MK_T_Cells_seurat_obj)

# Export T cell UMAP coordinates for Loupe Browser
umap_coords <- Embeddings(MK_T_Cells_seurat_obj, "umap")
cluster_annotations <- Idents(MK_T_Cells_seurat_obj)
export_data <- data.frame(Barcode = rownames(umap_coords), 
                         umap_coords, 
                         Cluster = as.integer(cluster_annotations))
write.csv(export_data, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/plots/t_cell_sub_population_annotation_helper_plots/res_1/loupe.csv", row.names = FALSE, quote = FALSE)


################################################################################
# SECTION 5: QUALITY CONTROL - REMOVE IRRELEVANT CLUSTERS
################################################################################

# Load the T cell object
MK_T_Cells_seurat_obj = readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# Remove low-quality or contaminating clusters from T cells
# These clusters were identified as doublets, low-quality cells, or non-T cells
# based on marker gene expression and quality metrics
irrelevant_clusters <- c(11,15,19,20)
relevant_clusters <- unique(MK_T_Cells_seurat_obj@meta.data$seurat_clusters)[
  !unique(MK_T_Cells_seurat_obj@meta.data$seurat_clusters) %in% irrelevant_clusters
]
MK_T_Cells_seurat_obj <- subset(MK_T_Cells_seurat_obj, 
                                subset = seurat_clusters %in% relevant_clusters)

# Save the cleaned T cell object for downstream analysis
saveRDS(MK_T_Cells_seurat_obj, 
       file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")

# Remove low-quality or irrelevant clusters from ALL cells
# These include doublets, low-quality cells, and contaminating cell types
irrelevant_clusters <- c(19,20,27, 15, 17, 18, 21, 41, 42)
relevant_clusters <- unique(resolution_1$seurat_obj@meta.data$seurat_clusters)[
  !unique(resolution_1$seurat_obj@meta.data$seurat_clusters) %in% irrelevant_clusters
]
MK_Cells_seurat_obj <- subset(resolution_1$seurat_obj, 
                              subset = seurat_clusters %in% relevant_clusters)

# Save the cleaned all-cells object for downstream analysis
saveRDS(MK_Cells_seurat_obj, 
       file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_Cells_seurat_obj.RDS")

################################################################################
# END OF SCRIPT
################################################################################
