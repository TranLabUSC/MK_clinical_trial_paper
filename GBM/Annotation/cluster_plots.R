library(Seurat)
library(ggplot2)
library(dplyr)

plotTopMarkersSeurat <- function(seurat_obj, n_top, output_dir) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the DimPlot
  p <- DimPlot(seurat_obj, label = TRUE, raster = FALSE)
  ggsave(file.path(output_dir, "raw.png"), plot = p, device = "png", width = 10, height = 8)
  
  # Identify all markers using FindAllMarkers
  markers <- FindAllMarkers(seurat_obj)
  
  # Loop over each cluster
  for (cluster in unique(markers$cluster)) {
    tryCatch({
      print(paste("Plotting for cluster:", cluster))
      
      # Get top n markers for the cluster based on avg_log2FC
      cluster_markers <- markers[markers$cluster == cluster, ]
      top_markers <- cluster_markers %>% top_n(n = n_top, wt = avg_log2FC)
      marker_genes <- top_markers$gene
      
      # Plot FeaturePlots for the top markers
      p <- FeaturePlot(seurat_obj, features = marker_genes, label = TRUE, cols = c("lightblue", "red"), raster = FALSE)
      
      # Save the plots in a PDF file
      pdf_filename <- file.path(output_dir, paste0("cluster_", gsub(" ", "_", cluster), ".pdf"))
      ggsave(pdf_filename, plot = p, device = "pdf", width = 20, height = 16)
      
    }, error = function(e) {
      # Print an error message and continue to the next cluster
      print(paste("Error plotting for cluster:", cluster, "Error:", e$message))
    })
  }
}

MK_NC_Monocyte_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")
n = 20
output_dir <- paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/NC_Monocyte_res_1_Marker_Plots/top_",n)

plotTopMarkersSeurat(seurat_obj = MK_NC_Monocyte_Cells_seurat_obj, n_top = n, output_dir = output_dir)


######################################################################################################################################
writeTopMarkersSeurat <- function(seurat_obj, n_top, output_dir) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Identify all markers using FindAllMarkers
  markers <- FindAllMarkers(seurat_obj)
  
  # Loop over each cluster
  for (cluster in unique(markers$cluster)) {
    tryCatch({
      print(paste("Writing for cluster:", cluster))
      
      # Get top n markers for the cluster based on avg_log2FC
      cluster_markers <- markers[markers$cluster == cluster, ]
      top_markers <- cluster_markers %>% top_n(n = n_top, wt = avg_log2FC)
      marker_genes <- top_markers$gene
      
      # Write the marker genes to a CSV file
      csv_filename <- file.path(output_dir, paste0("cluster_", gsub(" ", "_", cluster), ".csv"))
      write.table(marker_genes, file = csv_filename, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
      
    }, error = function(e) {
      # Print an error message and continue to the next cluster
      print(paste("Error writing for cluster:", cluster, "Error:", e$message))
    })
  }
}

# Load your Seurat object
MK_NC_Monocyte_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")

# Set the number of top markers and output directory
n = 100
output_dir <- paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/NC_Monocyte_res_1_Marker_Plots/top_", n)

# Run the function
writeTopMarkersSeurat(MK_NC_Monocyte_Cells_seurat_obj, n, output_dir)

