# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

# Define the function
plot_predecessor_distribution <- function(seurat_obj,
                                          target_cluster,
                                          target_timepoint,
                                          predecessor_timepoint,
                                          cluster_col = "seurat_clusters",
                                          timepoint_col = "timepoint",
                                          predecessor_col = "predecessor_barcode",
                                          plot_type = "bar") {
  # Check if required columns exist
  required_cols <- c(cluster_col, timepoint_col, predecessor_col)
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("The following required metadata columns are missing in the Seurat object:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Subset target cells based on cluster and timepoint
  target_cells <- seurat_obj@meta.data %>%
    filter((!!sym(cluster_col)) == target_cluster,
           (!!sym(timepoint_col)) == target_timepoint) %>%
    pull(barcode)
  
  if (length(target_cells) == 0) {
    stop("No cells found for the specified cluster and timepoint.")
  }
  
  # Extract predecessor barcodes for target cells
  predecessor_barcodes <- seurat_obj@meta.data[target_cells, predecessor_col]
  
  # Remove NA predecessors
  valid_indices <- !is.na(predecessor_barcodes)
  predecessor_barcodes <- predecessor_barcodes[valid_indices]
  target_cells_valid <- target_cells[valid_indices]
  
  if (length(predecessor_barcodes) == 0) {
    stop("No valid predecessors found for the specified target cells.")
  }
  
  # Ensure predecessor barcodes exist in the Seurat object
  existing_predecessors <- predecessor_barcodes %in% colnames(seurat_obj)
  if (sum(!existing_predecessors) > 0) {
    warning(paste(sum(!existing_predecessors), "predecessor barcodes not found in the Seurat object and will be excluded."))
    predecessor_barcodes <- predecessor_barcodes[existing_predecessors]
  }
  
  if (length(predecessor_barcodes) == 0) {
    stop("No valid predecessor barcodes found in the Seurat object.")
  }
  
  # Subset predecessor cells based on predecessor_timepoint
  predecessor_clusters <- seurat_obj@meta.data %>%
    filter(barcode %in% predecessor_barcodes,
           (!!sym(timepoint_col)) == predecessor_timepoint) %>%
    pull((!!sym(cluster_col)))
  
  if (length(predecessor_clusters) == 0) {
    stop("No predecessor clusters found for the specified predecessor timepoint.")
  }
  
  # Calculate distribution
  cluster_dist <- as.data.frame(table(predecessor_clusters))
  colnames(cluster_dist) <- c("Predecessor_Cluster", "Frequency")
  cluster_dist <- cluster_dist %>%
    arrange(desc(Frequency)) %>%
    mutate(Percentage = (Frequency / sum(Frequency)) * 100)
  
  # Plotting
  if (plot_type == "bar") {
    p <- ggplot(cluster_dist, aes(x = reorder(Predecessor_Cluster, -Percentage), y = Percentage)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = paste("Predecessor Distribution for Cluster", target_cluster),
           x = "Predecessor Cluster",
           y = "Percentage (%)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else if (plot_type == "pie") {
    cluster_dist <- cluster_dist %>%
      arrange(desc(Predecessor_Cluster)) %>%
      mutate(ypos = cumsum(Percentage) - 0.5 * Percentage)
    
    p <- ggplot(cluster_dist, aes(x = "", y = Percentage, fill = Predecessor_Cluster)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      labs(title = paste("Predecessor Distribution for Cluster", target_cluster)) +
      theme(legend.title = element_blank())
  } else {
    stop("Invalid plot_type. Choose either 'bar' or 'pie'.")
  }
  
  return(p)
}

seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")
# Specify parameters
target_cluster <- "4"            # Cluster of interest
target_timepoint <- "C1"              # Timepoint of target cells
predecessor_timepoint <- "Pre"         # Timepoint of predecessor cells
cluster_col <- "seurat_clusters"   # Name of cluster column in metadata
timepoint_col <- "TimePoint"       # Name of timepoint column in metadata
predecessor_col <- "predecessor_barcode" # Name of predecessor barcode column

# Generate bar plot
predecessor_plot <- plot_predecessor_distribution(seurat_obj = seurat_obj,
                                                  target_cluster = target_cluster,
                                                  target_timepoint = target_timepoint,
                                                  predecessor_timepoint = predecessor_timepoint,
                                                  cluster_col = cluster_col,
                                                  timepoint_col = timepoint_col,
                                                  predecessor_col = predecessor_col,
                                                  plot_type = "bar")

# Display the plot
print(predecessor_plot)

# Generate pie chart (optional)
predecessor_pie <- plot_predecessor_distribution(seurat_obj = seurat_obj,
                                                 target_cluster = target_cluster,
                                                 target_timepoint = target_timepoint,
                                                 predecessor_timepoint = predecessor_timepoint,
                                                 cluster_col = cluster_col,
                                                 timepoint_col = timepoint_col,
                                                 predecessor_col = predecessor_col,
                                                 plot_type = "pie")

# Display the pie chart
print(predecessor_pie)


###########################################################################################################################################################
# Load necessary libraries
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

# Define the function
perform_de_between_cluster_and_predecessors <- function(seurat_obj,
                                                        target_cluster,
                                                        target_timepoint,
                                                        predecessor_timepoint,
                                                        cluster_col = "seurat_clusters",
                                                        timepoint_col = "TimePoint",
                                                        predecessor_col = "predecessor_barcode",
                                                        min.pct = 0,
                                                        logfc.threshold = 0.25) {
  
  # Check if required columns exist
  required_cols <- c(cluster_col, timepoint_col, predecessor_col)
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("The following required metadata columns are missing in the Seurat object:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Subset target cells based on cluster and timepoint
  target_cells <- seurat_obj@meta.data %>%
    filter((!!sym(cluster_col)) == target_cluster,
           (!!sym(timepoint_col)) == target_timepoint) %>%
    pull(barcode)
  
  if (length(target_cells) == 0) {
    stop("No cells found for the specified cluster and timepoint.")
  }
  
  # Extract predecessor barcodes for target cells
  predecessor_barcodes <- seurat_obj@meta.data[target_cells, predecessor_col]
  
  # Remove NA predecessors
  valid_indices <- !is.na(predecessor_barcodes)
  predecessor_barcodes <- predecessor_barcodes[valid_indices]
  target_cells_valid <- target_cells[valid_indices]
  
  if (length(predecessor_barcodes) == 0) {
    stop("No valid predecessors found for the specified target cells.")
  }
  
  # Ensure predecessor barcodes exist in the Seurat object
  existing_predecessors <- predecessor_barcodes %in% colnames(seurat_obj)
  if (sum(!existing_predecessors) > 0) {
    warning(paste(sum(!existing_predecessors), 
                  "predecessor barcodes not found in the Seurat object and will be excluded."))
    predecessor_barcodes <- predecessor_barcodes[existing_predecessors]
  }
  
  if (length(predecessor_barcodes) == 0) {
    stop("No valid predecessor barcodes found in the Seurat object.")
  }
  
  # Identify predecessor cells at the specified predecessor timepoint
  predecessor_cells <- predecessor_barcodes[seurat_obj@meta.data[predecessor_barcodes, timepoint_col] == predecessor_timepoint]
  
  if (length(predecessor_cells) == 0) {
    stop("No predecessor cells found at the specified predecessor timepoint.")
  }
  
  # Create a new identity class: 'Target' vs 'Predecessor'
  cells_for_de <- c(target_cells_valid, predecessor_cells)
  seurat_subset <- subset(seurat_obj, cells = cells_for_de)
  
  # Define identity classes
  seurat_subset$de_group <- ifelse(colnames(seurat_subset) %in% target_cells_valid, "Target", "Predecessor")
  
  # Set the identity for DE analysis
  Idents(seurat_subset) <- "de_group"
  
  # Perform differential expression analysis
  de_results <- FindMarkers(seurat_subset, 
                            ident.1 = "Target", 
                            ident.2 = "Predecessor", 
                            min.pct = min.pct, 
                            logfc.threshold = logfc.threshold)
  
  # Add gene names as a column
  de_results <- de_results %>%
    rownames_to_column(var = "gene")
  
  # Define the output filename
  output_filename <- paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/cluster_", target_cluster, "_de_result.csv")
  
  # Save DE results to CSV
  write.csv(de_results, file = output_filename, row.names = FALSE)
  
  message(paste("Differential expression results saved to", output_filename))
  
  # Return DE results and cell lists
  return(list(
    de_results = de_results,
    target_cells_valid = target_cells_valid,
    predecessor_cells = predecessor_cells
  ))
}

# Example Usage

# Load your updated Seurat object
# Replace 'updated_seurat_object.rds' with your actual file path
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")

# Specify parameters
target_cluster <- "1"                  # Cluster of interest (as a string)
target_timepoint <- "C1"                    # Timepoint of target cells
predecessor_timepoint <- "Pre"               # Timepoint of predecessor cells
cluster_col <- "seurat_clusters"         # Name of cluster column in metadata
timepoint_col <- "TimePoint"             # Name of timepoint column in metadata
predecessor_col <- "predecessor_barcode" # Name of predecessor barcode column

# Perform differential expression analysis and capture results
results <- perform_de_between_cluster_and_predecessors(
  seurat_obj = seurat_obj,
  target_cluster = target_cluster,
  target_timepoint = target_timepoint,
  predecessor_timepoint = predecessor_timepoint,
  cluster_col = cluster_col,
  timepoint_col = timepoint_col,
  predecessor_col = predecessor_col,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Access DE results
de_results <- results$de_results

# Access target and predecessor cells
target_cells_valid <- results$target_cells_valid
predecessor_cells <- results$predecessor_cells

# Proceed with expression data extraction
top_genes <- de_results %>%
  arrange(p_val_adj) %>%
  slice(1:20) %>%
  pull(gene)

# Extract expression data for these genes
expr_data <- FetchData(seurat_obj, vars = top_genes, cells = c(target_cells_valid, predecessor_cells))

# Create annotation for DE groups
annotation <- data.frame(Group = ifelse(colnames(expr_data) %in% target_cells_valid, "Target", "Predecessor"))
rownames(annotation) <- colnames(expr_data)

# Generate heatmap
library(pheatmap)

pheatmap(as.matrix(expr_data), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = annotation,
         show_colnames = FALSE,
         show_rownames = TRUE,
         scale = "row",
         main = paste("Top 20 DE Genes for Cluster", target_cluster))

ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = paste("Volcano Plot for Cluster", target_cluster),
       x = "Average Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "red")




############################################################################################################################################################################################
# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)

setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis")

# 1. Load the Seurat object
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")

# Check the metadata
metadata <- seurat_obj@meta.data
head(metadata)

# Ensure that 'TimePoint' and 'Patient' columns exist
if(!all(c("TimePoint", "Patient") %in% colnames(metadata))) {
  stop("The Seurat object does not contain 'TimePoint' and/or 'Patient' columns in metadata.")
}

# 2. Prepare pseudo-bulk RNA-seq counts
# Aggregate counts by Patient and TimePoint

# Extract the raw counts matrix
# Commonly in the "RNA" assay
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Add Patient and TimePoint information to the metadata
metadata <- seurat_obj@meta.data %>%
  select(Patient, TimePoint)

# **Modified Step:** Create a unique identifier for each pseudo-bulk sample with "patient_" prefix
metadata$Sample <- paste("patient", metadata$Patient, metadata$TimePoint, sep = "_")
# Example: "patient_1_Pre", "patient_1_C1", "patient_2_Pre", etc.

# Sum the counts for cells belonging to the same Sample
# Transpose counts to have cells as rows and genes as columns for aggregation
counts_df <- as.data.frame(t(as.matrix(counts)))
counts_df$Sample <- metadata$Sample

# Aggregate counts by Sample
pseudo_bulk_counts <- counts_df %>%
  group_by(Sample) %>%
  summarise_all(sum)

# **Corrected Step:** Convert back to a matrix with genes as rows and samples as columns
# Since 'pseudo_bulk_counts' has samples as rows and genes as columns,
# we need to transpose it to have genes as rows and samples as columns.

# Remove the 'Sample' column and transpose
pseudo_bulk_counts_mat <- t(as.matrix(pseudo_bulk_counts[,-1]))

# Assign gene names as rownames
rownames(pseudo_bulk_counts_mat) <- rownames(counts)  # Correct: Gene names

# Assign sample names as column names
colnames(pseudo_bulk_counts_mat) <- pseudo_bulk_counts$Sample

# Verify dimensions
dim(pseudo_bulk_counts_mat)  # Should be genes x samples (e.g., 19343 x 24)

write.csv(pseudo_bulk_counts_mat, file = "pseudo_bulk_RNASeq.csv")

# 3. Prepare the design matrix for edgeR
# Extract Patient and TimePoint information for each Sample
sample_info <- data.frame(
  Sample = colnames(pseudo_bulk_counts_mat)
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_")

# **Modified Step:** Reconstruct the Patient identifier with the "patient_" prefix
sample_info$Patient <- paste(sample_info$Patient_Prefix, sample_info$Patient_ID, sep = "_")
# Example: "patient_1", "patient_2", etc.

# Ensure factors are correctly set
sample_info$Patient <- factor(sample_info$Patient)
sample_info$TimePoint <- factor(sample_info$TimePoint, levels = c("Pre", "C1"))

# 4. Create DGEList object
dge <- DGEList(counts = pseudo_bulk_counts_mat)

# 5. Normalize the data
dge <- calcNormFactors(dge)

# 6. Create the design matrix
# Model: ~ Patient + TimePoint
design <- model.matrix(~ Patient + TimePoint, data = sample_info)

# 7. Estimate dispersion
dge <- estimateDisp(dge, design)

# 8. Fit the model and perform differential expression
fit <- glmFit(dge, design)
# Assuming "TimePointC1" is the coefficient for C1 vs Pre
lrt <- glmLRT(fit, coef = "TimePointC1")

# 9. Extract differential expression results
de_results <- topTags(lrt, n = Inf)$table
de_results <- as.data.frame(de_results)
# Add gene names as a column
de_results$Gene <- rownames(de_results)

# 10. Save the differential expression results
write.csv(de_results, file = "differential_expression_results.csv", row.names = FALSE)

# 11. Create input file for GSEA
# Typically, GSEA requires a ranked list of genes. One common ranking metric is logFC multiplied by the sign of the p-value.

# Here, we'll use the signed -log10(p-value) multiplied by the sign of logFC
# Alternatively, you can choose other ranking metrics based on your preference

# To handle p-values of 0, add a small pseudocount
de_results$PValue[de_results$PValue == 0] <- 1e-300

de_results$RankingMetric <- -log10(de_results$PValue) * sign(de_results$logFC)

# Order the genes by the ranking metric
ranked_genes <- de_results %>%
  arrange(desc(RankingMetric)) %>%
  select(Gene, RankingMetric)

# Save the ranked list for GSEA
# **Option 1:** GCT format
# GCT format requires specific headers. Here's how to structure it:

# Create a GCT header
gct_header <- data.frame(
  NAME = c("#1.2"),
  DESCRIPTION = c("RankingMetric"),
  stringsAsFactors = FALSE
)

# Combine header and data
gct_data <- data.frame(
  NAME = ranked_genes$Gene,
  DESCRIPTION = ranked_genes$Gene,  # GSEA often expects gene symbols or descriptions here
  RankingMetric = ranked_genes$RankingMetric,
  stringsAsFactors = FALSE
)

# Write the GCT file
write.table(gct_header, file = "GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gct_data, file = "GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE)

# **Option 2:** RNK format
# RNK format is simpler and often preferred for ranked lists

write.table(ranked_genes, file = "GSEA_ranked_genes.rnk", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Inform the user that the process is complete
cat("Pseudo-bulk RNA-seq aggregation, differential expression analysis, and GSEA input file creation are complete.\n")
cat("Differential expression results saved to 'differential_expression_results.csv'.\n")
cat("GSEA ranked gene list saved to 'GSEA_ranked_genes.gct' and 'GSEA_ranked_genes.rnk'.\n")





######################################################################################
library(dplyr)
library(tidyr)
# 1. Read the Pseudo-Bulk Counts Matrix
# -------------------------------------

# Specify the path to your pseudo_bulk_counts_mat.csv
counts_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/pseudo_bulk_RNASeq.csv"

# Read the CSV file
# Assuming the first column contains gene names and is set as rownames
counts_df <- read.csv(counts_csv_path, row.names = 1, check.names = FALSE)

# Verify the dimensions and a snippet of the data
cat("Dimensions of the counts matrix:", dim(counts_df), "\n")
head(counts_df)


# 4. Normalize Counts Using "Pre" as Baseline
# -------------------------------------------

# Subset the counts matrix to include only selected genes
selected_counts <- counts_df

# Parse sample names to extract Patient and TimePoint information
# Assuming sample names are in the format "patient_X_TimePoint", e.g., "patient_1_Pre"
sample_names <- colnames(selected_counts)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = sample_names,
  stringsAsFactors = FALSE
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_", remove = FALSE) %>%
  mutate(Patient = paste(Patient_Prefix, Patient_ID, sep = "_")) %>%
  select(Sample, Patient, TimePoint)

# View sample information
print(sample_info)

# Get a list of unique patients
patients <- unique(sample_info$Patient)
cat("Number of unique patients:", length(patients), "\n")

# Initialize a list to store log2 fold changes
log2fc_list <- list()

# Loop through each patient to compute log2FC
for (patient in patients) {
  # Identify "Pre" and "C1" samples for the patient
  pre_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "Pre") %>%
    pull(Sample)
  
  c1_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "C1") %>%
    pull(Sample)
  
  # Check if both samples are present
  if (length(pre_sample) == 1 & length(c1_sample) == 1) {
    # Extract counts for "Pre" and "C1" samples
    pre_counts <- selected_counts[, pre_sample]
    c1_counts <- selected_counts[, c1_sample]
    
    # Compute log2 fold change with a pseudocount of 1 to avoid division by zero
    log2fc <- log2((c1_counts + 1) / (pre_counts + 1))
    
    # Ensure that log2fc is a named vector with gene names
    names(log2fc) <- rownames(selected_counts)
    
    # Store in the list
    log2fc_list[[patient]] <- log2fc
  } else {
    warning(paste("Missing 'Pre' or 'C1' sample for patient:", patient))
  }
}

# Combine the list into a normalized matrix
# Rows: Genes, Columns: Patients (log2FC)
normalized_mat <- do.call(cbind, log2fc_list)

# Rename columns to indicate log2FC
colnames(normalized_mat) <- paste0(colnames(normalized_mat), "_log2FC")

# Inspect the normalized matrix
dim(normalized_mat)
head(normalized_mat)

# 5. Scale the Normalized Data
# ----------------------------

# Scaling can be done per gene to standardize the log2FC across patients
# This is useful for visualization purposes

# Transpose, scale, and transpose back
scaled_mat <- t(scale(t(normalized_mat)))

# Replace any NA values resulting from scaling with 0
scaled_mat[is.na(scaled_mat)] <- 0

# Inspect the scaled matrix
dim(scaled_mat)
head(scaled_mat)

# 6. Save the Normalized and Scaled Data as CSV
# ---------------------------------------------

# Prepare the normalized data for saving
normalized_df <- as.data.frame(normalized_mat)
normalized_df$Gene <- rownames(normalized_df)

# Reorder columns to have 'Gene' as the first column
normalized_df <- normalized_df %>%
  select(Gene, everything())

# Save the normalized log2FC data
write.csv(normalized_df, file = "normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Normalized log2 fold change data saved to 'normalized_log2FC_C1_vs_Pre.csv'.\n")

# Prepare the scaled data for saving
scaled_df <- as.data.frame(scaled_mat)
scaled_df$Gene <- rownames(scaled_df)

# Reorder columns to have 'Gene' as the first column
scaled_df <- scaled_df %>%
  select(Gene, everything())

# Save the scaled normalized data
write.csv(scaled_df, file = "scaled_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Scaled normalized data saved to 'scaled_normalized_log2FC_C1_vs_Pre.csv'.\n")



#########################################################################################################################################################
# Plotting heatmap

# Install necessary packages if not installed
# install.packages("pheatmap")

library(pheatmap)

# 1. Read the file
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/scaled_normalized_log2FC_C1_vs_Pre.csv", header = TRUE, stringsAsFactors = FALSE)

# The columns are expected to be:
# Gene, patient_10_log2FC, patient_12_log2FC, patient_13_log2FC, patient_14_log2FC,
# patient_18_log2FC, patient_19_log2FC, patient_20_log2FC, patient_21_log2FC,
# patient_2_log2FC, patient_3_log2FC, patient_5_log2FC, patient_7_log2FC

# 2. Define the gene list
gene_list <- c("ENSG00000269242", "IGLV1-40", "COL19A1", "GLB1L", "ENSG00000285542", "PKP2", "OIP5", 
               "MATR3", "VSIG4", "GAREM1", "CD163", "RNF227", "ALOX12", "TRAT1", "CD5", "CXCL3", "PKHD1L1", 
               "TAL1", "HBA1", "CXCL2", "PPARG", "EGR1", "SOCS3", "TEX14", "ANO5", "SIGLEC11", "CD163L1", 
               "ABTB2", "STOML1", "GPR153", "KLHDC7B", "LIME1", "H1-5", "RNF144A", "PARPBP", "TRAC", 
               "MAL", "ZFPM2", "IGF2BP1", "TPX2", "FGFBP2", "ANK2", "WDR62", "LHFPL5", "SDK1", "H3C2", 
               "CDC14B", "TRBC2", "STIL", "PRKG1", "H2AC11", "DYNC2H1", "ZNF532", "CENPF", "ENAH", "TAF9B", 
               "LYPD2", "ADGRB3", "SLIT2", "RAB27B", "BACH2", "ITGB7", "KLRC1", "MDGA2", "RGS7", "FCER2", 
               "H2BC15", "ENSG00000258311", "ELOVL7", "MAP7", "C1orf198", "MCM6", "ASAP2", "ST8SIA6", 
               "GPR174", "OXCT1", "HCAR3", "TIGIT", "C1QC", "LDHD", "VMO1", "PAX8", "TSHZ2", "FOS", "GPR183", 
               "NR4A2", "JUN", "PLEKHF1", "FHL1", "BMP6", "H2BC8", "SH3BGRL2", "SCN1B", "PCSK6", "CABP5", 
               "CXCL5", "BEND2", "ZNF417", "MFAP3L", "S1PR5", "HOMER2", "P2RY12", "NT5M", "CTTN", "PRKAR1B", 
               "THBS1", "SNCA", "PROS1", "ESAM", "LGALSL", "CTDSPL", "PDGFA", "PTCRA", "ITGA2B", "H2BC9", 
               "DAB2", "TSPAN33", "TRIM58", "PPBP", "GP9", "ITGB3", "C19orf33", "PF4", "SPARC", "F13A1", 
               "CAVIN2", "ACRBP", "MYL9", "CLU", "TSC22D1", "STON2", "SLC40A1", "LTBP1", "H2BC11", "FIBCD1", 
               "ABLIM3", "PF4V1", "TREML1", "TMEM40", "MMD", "GNG11", "TUBB1", "CMTM5", "MPIG6B", "RHOBTB1", 
               "IL21R", "DMTN", "PTGS1", "SAV1", "PANX1", "H4C8", "PLA2G12A", "PDLIM1", "H3C10", "BEX3", 
               "SLC12A4", "BACE1", "PLTP", "GATD3A", "MTG2", "LONRF3", "UBXN8", "LMNA", "ZNF358", "ACBD4", 
               "ENSG00000271741", "TNFRSF10D", "SIK1B", "NUDT6", "TNFAIP1", "CCL4", "ENSG00000205045", "NUP35",
               "CAMKK1", "ENSG00000285827", "FLT3", "MARVELD1", "RGL1", "PER1", "FMN1", "IL1R2", "ZBTB16", 
               "DDIT4", "FKBP5", "DGKH", "NAA80", "BIN1", "LEF1", "KIF5C", "CATSPER1", "HOXB2", "SH2D1A", 
               "C2orf88", "YPEL1", "ZNF707", "NUDCD1", "RPH3A", "CACNB4", "SACS", "PPOX", "ADCY3", "AMOT", 
               "NBEA", "CACNA2D3", "GGACT", "RNLS", "KCTD18", "RAB13", "OPN3", "CHD7", "OTUD3", "HSPA1B", 
               "NOXA1", "ZBTB49", "PTGDR2", "TIGD1", "RETREG1", "EPHB2", "CTTNBP2", "SLC35C1", "HEMK1", 
               "SLC2A9", "ATP8B4", "CLPB", "MAFG", "PIP4P2", "MRPS17", "CYP27A1", "MROH6", "MEI1", "C17orf80",
               "PDRG1", "ZNF432", "CNNM2", "CRISPLD2", "FOLR3", "BLZF1", "PTK2", "MVB12B", "SS18L1", "FASTKD1",
               "ENSG00000260729", "CLECL1", "NRG1", "ACTR3B", "TMEM150A", "CACNA1A", "ASGR2", "SLC38A7", 
               "GPRC5C", "RAB39A", "GLT1D1", "CES1", "CALCRL", "HOMER3", "SESN2", "CCR2", "MAF", "ZP3", 
               "ZNF354A", "TP53I11", "TRAF1", "CD247", "PAQR7", "FCRL6", "ARL4A", "SAP30", "PRKCQ", "SPON2", 
               "CPM", "FPR2", "ENSG00000275464", "ZNF331", "PADI4", "TAF4B", "ODAD3", "ZNF410", "PTMS", 
               "MMP17", "CLEC10A", "DYSF", "NOL12", "CAMK4", "LRRC57", "BCL11B", "COQ8A", "MAP2K6", "RXRB", 
               "ICAM1", "PRR5", "PPIF", "MICAL2", "ENSG00000260272", "EVA1B", "AP1S1", "SATB1", "RETN", 
               "ATP6V0A1", "MT1X", "FAM20C", "RFX2", "ATF5", "CPED1", "RELB", "COIL", "TOR1B", "OSBPL5", 
               "METTL8", "KIF13A", "MSR1", "PLCB1", "RNASE6", "IL13RA1", "GNA15", "SIRPA", "PID1", "GSTM1", 
               "CSNK2B", "SNRPD3", "EEF1G", "RNASE2", "TFEC", "CD93", "VGLL4", "MS4A6A", "RBP7", "TBC1D2",
               "ALDH2", "LGALS2", "ASGR1", "ZNF467", "ITGAM", "FOLR2", "ADAM15", "DOCK4", "IL1RN", "C1QB", 
               "C1QA", "NCF1", "TCN2", "FUOM", "HIP1", "PLBD1", "MARCO", "FES", "NAIP", "TMEM91", "SLC35D2", 
               "IRS2", "MLXIP", "HBEGF", "SLAMF7", "OAS3", "EPSTI1", "GBP5", "MEN1", "EHD4", "AMDHD2", 
               "FCF1", "RILP", "TRIM69", "TMEM106A", "FCGR1A", "ENSG00000124593", "AGAP3", "STAB1", "PTAFR", 
               "PXN", "ZNF101", "SNX25", "TMEM128", "LYSMD3", "RRP15", "RMND1", "POLR1G", "TMEM42", "PDE7B", 
               "C3", "DHX58", "ZNF224", "SPIRE1", "PPM1N", "DTD1", "PCGF1", "ICA1L", "ST7L", "SAMD3", "ETS1", 
               "DNHD1", "PLD4", "INKA2", "HSH2D", "TAGAP", "CTSL", "CYP4F22", "TRPT1", "PTPN4", "KIAA0586", 
               "RAB37", "NKG7", "LBH", "UPK3A", "PCBP3", "PRKN", "ZC3H10", "ABCA7", "CBR4", "CASP3", "NELFA", 
               "TMEM104", "TXN", "RARS2", "EI24", "TMEM80", "GLT8D1", "CBLB", "ZRANB3", "INSR", "DUSP1", 
               "NFKBIA", "ZFP36", "ZMYND11", "AHR", "ISG20", "RHBDD2", "CYSTM1", "PCED1B", "AQP9", "AKT3", 
               "PAPSS2", "TTC39B", "RPL41", "IFITM1", "TRMT61A", "CD101", "TCEAL4", "ST3GAL6", "SNED1", "NID1",
               "MLLT1", "ASB6", "SZT2", "CAPN10", "TRAF5", "NRGN", "PRKAR2B", "SKAP1", "KIF2A", "SYNE1", 
               "MAP3K7CL", "C10orf143", "H2BC12", "PDCL", "C5orf22", "NUDT18", "LSS", "SLC35E3", "GFOD1", 
               "MTPAP", "PEX5", "SMG9", "CRIPT", "SMARCD3", "U2AF1L4", "ZBED5", "NPIPB4", "ANKAR", "DDB2", 
               "CISD1", "ANKZF1", "H2BC4", "RPS6KA2", "DLG4", "PLEKHA7", "MGST1", "HPF1", "ACTN1", "LTB", 
               "TMEM102", "TSTD1", "ZDHHC2", "SH2D1B", "PGRMC1", "H2AC6", "SPHK2", "ASPH", "CAMTA1", "NBPF12",
               "ISG15", "HEG1", "SLC25A20", "ARID5B", "KLHL22", "SLFN5", "IFI6", "ZC2HC1A", "TBC1D15", 
               "L3MBTL3", "ITPKB", "ERCC5", "H1-0", "KANSL1L", "DUS2", "NAA15", "CYSLTR1", "APP", "SACM1L", 
               "MT-ATP8", "FLT3LG", "TMEM154", "TNFRSF8", "ACAP1", "ARL6IP6", "ANKRD36", "PLCL1", "DALRD3", 
               "RBM43", "TSEN15", "PLAGL1", "SAMSN1", "ERGIC2", "LMBRD1", "RHOF", "RPGR", "RNF146", "GPR160", 
               "NUDT15", "CEP85", "TTC9", "FAM200B", "ZMAT3", "C16orf54", "SYAP1", "SNRNP35", "ELL2", "MYBBP1A",
               "SHTN1", "DBNDD2", "GPAT4", "ABHD16A", "CHST2", "DYRK2", "PTP4A3", "PRR5L", "CX3CR1", "DUSP5", 
               "IVD", "SETBP1", "CSF3R", "MKKS", "TRAPPC5", "PLD3", "RPS17", "S100A12", "S100A8", "MLEC", "IPO5",
               "CARD19", "MTMR3", "CD14", "SASH1", "VCAN", "FRMD4B", "CRTAP", "FCN1", "LYZ", "CFP", "ARHGEF10L",
               "S100A9", "RBM47", "GPX1", "C1RL", "MERTK", "DCAF12", "SLC35F6", "ETHE1", "JARID2", "ARHGAP31",
               "PCBD1", "ARHGAP24", "CPNE8", "SPRED1", "CLIC4", "TAF2", "ATP2A2", "NUP210", "NR6A1", "PEX16",
               "BRD3", "ZMYM4", "SNTB1", "PARVB", "UBASH3B", "CREB5", "TBC1D9", "CBX6", "TNK2", "GOLIM4", 
               "TPM1", "CD36", "TRPS1", "LIMK2", "TNS3", "MCFD2", "SIRT6", "FBN2", "ZNF33A", "TST", "MXD1",
               "SLC2A3", "DGKG", "MTMR11", "TACC3", "IL6ST", "SMIM7", "TBCD", "EIF4G3", "XPO6", "HLCS", 
               "PITPNB", "NFU1", "HAGH", "SNAPC3", "RNF24", "BICRAL", "RASA1", "ITGAE", "SH3GLB2", "PLA2G4A", 
               "ENTPD1", "ANPEP", "CALML4", "MTMR10", "PMM2", "CAMTA2", "CD244", "S100Z", "SMAD1", "APOBEC3C", 
               "RBM22", "NFIL3", "AARSD1", "ACSL1", "TNFAIP3", "NAMPT", "KLF10", "B4GALT1", "EMILIN2", "LY86",
               "NME2", "HLA-DMB", "VAMP8", "HLA-DMA", "C3AR1", "C18orf32", "CTSB", "TNFRSF14", "ODF3B", "PIGQ",
               "MSL1", "TRIM22", "BLVRB", "LITAF", "LRP10", "SHISA5", "ADIPOR1", "ERAP1", "ATP1B3", "HSP90B1",
               "UBC", "SERPINA1", "SARAF", "CIRBP", "ARL6IP5", "RPL29", "EEF1B2", "CSTA", "HLA-DRB5", "AP1S2",
               "CCDC12", "UQCRH", "RPLP2", "RPLP1", "RPL13", "RPL18", "RPS24", "RPL30", "RPS8", "RPS18", 
               "RPS14", "RPL19", "RPL5", "RPL9", "RPL12", "RPS13", "RPS12", "RPL32", "RPL10", "RPL7A", "RPL18A",
               "RPL6", "RPS5", "RPS3A", "RPL24", "RPLP0", "RPL15", "RPS6", "GAPDH", "RPL23", "ETS2", "COPG2", 
               "SESN1", "FMNL2", "XPR1", "CEBPD", "JDP2", "FOXO3", "GSN", "CROCC2", "POLM", "FUCA1", "TFEB", 
               "STK25", "SNX20", "NBPF26", "TMBIM4", "SARNP", "RAB3D", "PCSK7", "ADAT1", "ALDH16A1", "TAPBPL",
               "SPSB3", "SLC66A3", "RBM4", "SNX18", "OAS2", "TMEM176B", "IFIT2", "KNDC1", "ADIPOR2", "PACS2",
               "EHBP1", "PAICS", "GLO1", "TTC19", "TRUB2", "RABEP1", "TIMM50", "PCK2", "AGPAT3", "ARHGAP45",
               "ITGAL", "LPP", "NCOA3", "RNASEK", "TNFRSF1B", "RHOC", "MYO1G", "LFNG", "CRIP1", "ARAP1",
               "PKN1", "INSIG1", "KLF4", "CCM2", "LNPEP", "GMIP", "CD58", "KIAA2026", "SH2D3C", "SPN", "PIK3CG",
               "PHTF2", "CD79B", "DRAP1", "USF3", "LYN", "SVIL", "RAP1GAP2", "ARAP2", "RAB29", "KLHL5", "CERK",
               "MXD3", "HSPA5", "SUN2", "RUNX3", "RBM38", "S1PR4", "CDKN1C", "VASP", "TUBA1A", "CEBPA", "MCM5",
               "GYPC", "HSPH1", "DNAJC10", "TRMT1", "PLAAT4", "MNDA", "CISD3", "CD52", "UCHL3", "METTL9",
               "NDUFB11", "YWHAG", "XPC", "HSPE1", "RPL22L1", "GLRX", "RPL39", "MRPL52", "RPS28", "TOMM5",
               "COMMD6", "IER3IP1", "C11orf21", "PFDN1", "ZDHHC1", "KLF12", "IPMK", "TBC1D10C", "APOBEC3G",
               "LIPA", "KLF3", "TRAF3IP3", "DDT", "ESYT1", "ZFP36L1", "PPP1R14B", "UBE2D1", "SPIDR", "MS4A4E",
               "MCL1", "PMVK", "DNAJB11", "NFE2L2", "EEA1", "SLC15A4", "FOXN3", "TMEM170B", "HSP90AA1", "QKI",
               "DYNLL1", "PISD", "FNDC3B", "WDFY3", "CMIP", "PSMB5", "YTHDF1", "IGF2R", "APEX1", "DDX21",
               "HEBP2", "HSPD1", "EPB41L3", "PTPRE", "COMT", "YBX3", "EMC6", "FERMT3", "AP2S1", "INPP4A",
               "GRK2", "ADK", "WDR43", "SSBP2", "MPC1", "PTPN1", "ZMYND8", "LEPROTL1", "ID2", "STIM2", "IRF7",
               "GBP4", "FGR", "TNFSF10", "GBP2", "LY6E", "WARS1", "MAD1L1", "SLC44A2", "PFKL", "TESC", "STIP1",
               "CCNDBP1", "OST4", "DHRS7", "FAM110A", "KLF2", "CFD", "SYNGR2", "FCGR3A", "RAC2", "UNC119",
               "CTSC", "OAS1", "CXCL16", "LAMTOR1", "SOD1", "TMC6", "UROS", "CDH23", "SETDB1", "ELOVL1", 
               "SF3B4", "PELI1", "CEBPB", "ARHGEF18", "DPEP2", "SLC31A2", "PML", "ZDHHC12", "SELPLG", "SSBP4",
               "HLA-E", "UNC93B1", "SIGLEC10", "WDR1", "IL16", "LEMD2", "DEDD2", "CCNL1", "NUTF2", "NOSIP",
               "LIMD2", "MAPKAPK3", "ADGRE2", "VSIR", "MTSS1", "C5AR2", "N4BP2L2", "CSTB", "TBCB", "PHF19", 
               "SP110", "SYTL1", "PECAM1", "CUX1", "HDGF", "NAP1L1", "PTPN6", "LYST", "GNG2", "FAM126A", 
               "RASAL3", "IPO7", "ZNF706", "STXBP2", "CANX", "HSP90AB1", "SPINT2", "ICAM3", "ADGRE3", "PPFIA1",
               "NUDC", "TMEM120A", "CD99", "ST3GAL1", "INPP5D", "CTBP2", "FRY", "TPST2", "LPCAT3", "RNF11",
               "ATP1A1", "DIP2A", "CYFIP2", "COQ2", "ANKRD44", "FGD3", "YWHAE", "MDFIC", "POLD4", "GPAA1",
               "TXNDC17", "GOLGA4", "PRKX", "LAT2", "HLA-F", "TRAPPC12", "LIMD1", "TPM4", "ARHGEF1", "UCP2",
               "RNF213", "PITPNC1", "LGALS9", "CD300LF", "ADGRE5", "UNC13D", "ITGAX", "MAP3K1", "CALHM2", 
               "EZR", "FLNA", "CSF1R", "GBP1", "JAK3", "LDB1", "TTYH3", "RFTN1", "TCIRG1", "SLC2A6", "VASH1",
               "PDIA4", "HEXD", "UBE2S", "SLC35E1", "METTL22", "IRAK4", "GRAMD1A", "VPS26B", "DNAJC2", "CHMP1A",
               "GALNT1", "TRAFD1", "SOX4", "RMC1", "ADA", "HES4", "CCDC115", "INPP5K", "CDIP1", "QPRT", "PRPF4",
               "SLC8B1", "NIP7", "GRK5", "KLF9", "DHPS", "NDST1", "PLSCR3", "BCAP29", "SLC39A1", "RIPK3", 
               "ENSG00000285589", "ENSG00000263620", "GATD3B", "PRDM1", "CD38", "GMEB2", "HPSE", "RNF25", 
               "CD22", "UTP4", "DNAJB9", "AARS1", "SYT17", "PTGDR", "PGLYRP2", "GMPR", "TCEAL3", "SVIP",
               "DNM3", "ARMCX5", "SMIM3", "CD9", "MTURN", "NKIRAS1", "GSTA4", "BRCA2", "SCARB1", "H4C5",
               "ZNF627", "TRIM3", "MTHFD2L", "GSKIP", "HAUS6", "BLMH", "HAL", "TMEM106C", "SIRT5", "FCRL3",
               "DHRS12", "STAT4", "B3GNT8", "RET", "NCR3", "CYP2U1", "BNIP3", "CISH", "ADAM8", "OPTN", 
               "NCALD", "TBX21", "P2RY10", "TOX", "TXK", "IL2RB", "PRF1", "GZMB", "SH2D2A", "CD96", "PTPRCAP", 
               "CD7", "NCKAP1", "GZMA", "CST7", "CTSW", "GNLY", "GZMM", "GZMH", "CLIC3", "AQP3", "KLRB1",
               "TRBC1", "FEN1", "TGFBR3", "FAM110B", "RGL4", "NME4", "CKB", "SFTPD", "RSKR", "TBC1D19", "HHAT",
               "FCMR", "SCML4", "GPC6", "MPP6", "S1PR1", "IL7R", "WDR27", "SPOCK2", "ZAP70", "GATA3", "MYBL1",
               "ABLIM1", "LCK", "IL32", "CD3E", "TCF7", "TMOD1", "LPL", "GPER1", "DEPTOR", "APOO", "IGLV3-19",
               "IGKV4-1", "TRBV4-2", "IL6R", "IL18", "IL27RA")

# 3. Subset the data for only the requested genes
subset_data <- data[data$Gene %in% gene_list,]

# Define patient order
patient_order <- c(18, 12, 14, 7, 10, 5, 20, 21, 19, 2, 3, 13)

# Create column names based on the patient order
ordered_columns <- paste0("patient_", patient_order, "_log2FC")

# Check that all ordered columns exist in the data
# If any column doesn't exist, you might need to verify the column naming in the CSV.
all(ordered_columns %in% colnames(subset_data))

# Reorder the data columns to match the desired patient order
subset_data <- subset_data[, c("Gene", ordered_columns)]

# Set rownames to Gene, remove Gene column for heatmap matrix
rownames(subset_data) <- subset_data$Gene
subset_data <- subset_data[,-1]  # remove the Gene column

# Convert to a numeric matrix if needed
mat <- as.matrix(subset_data)

# 4. Plot the heatmap
# You can adjust clustering, scaling, color palettes, etc. as needed.
# Run pheatmap and store the result in an object
res <- pheatmap(mat, 
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "row",
                show_rownames = TRUE, 
                show_colnames = TRUE,
                main = "Heatmap of Selected Genes",
                fontsize_row = 6,
                fontsize_col = 8)

# Extract the hierarchical clustering tree for rows
hc_rows <- res$tree_row

# Choose the number of clusters (e.g., k=5). This number can be adjusted based on 
# the structure you see in the dendrogram.
k <- 20

# Cut the dendrogram into k clusters
clusters <- cutree(hc_rows, k = k)

# Create a data frame with gene and cluster assignments
gene_clusters <- data.frame(
  Gene = rownames(mat),
  Cluster = clusters
)

# Write the cluster assignments to a CSV file so you can inspect them outside R
write.csv(gene_clusters, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/heatmap_gene_clusters.csv", row.names = FALSE)


# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Assuming 'mat' is already created from the previous steps:
# 'mat' is a numeric matrix with rows = genes and columns = patients.
# Columns should be named like "patient_18_log2FC", "patient_12_log2FC", etc.

# Define patient groups
short_term_survivors <- c(18, 12, 14, 7, 10)
long_term_survivors <- c(5, 20, 21, 19, 2, 3, 13)

# Extract the column names from mat
all_patients <- colnames(mat)

# Check that all patient columns are in the expected format: "patient_#_log2FC"
# If they are, we can parse them to find the patient numbers
# Let's extract patient numbers from column names:
# Example column name: "patient_18_log2FC"
patient_numbers <- sapply(strsplit(all_patients, "_"), function(x) x[2])
patient_numbers <- as.numeric(patient_numbers)

# Calculate the pathway activity (average expression of all genes in the pathway per patient)
# Since 'mat' rows are genes and columns are patients, we take the column means.
pathway_activity <- colMeans(mat, na.rm = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(
  patient = patient_numbers,
  activity = pathway_activity
)

# Add a group column
plot_data$group <- ifelse(plot_data$patient %in% short_term_survivors, "Short-term", "Long-term")

# Convert group to factor for plotting
plot_data$group <- factor(plot_data$group, levels = c("Short-term", "Long-term"))

# Perform a statistical test (Wilcoxon rank-sum test)
stat_test <- wilcox.test(activity ~ group, data = plot_data)

# Extract the p-value
p_value <- stat_test$p.value

# Create the box plot to match the template style
p <- ggplot(plot_data, aes(x = group, y = activity, fill = group)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 1.2), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), size = 4, alpha = 0.8) +
  scale_fill_manual(values = c("Short-term" = "blue", "Long-term" = "red")) +
  labs(
    title = "Pathway Activity by Survival Group",
    x = NULL,  # Remove x-axis label
    y = "Pathway Activity (Mean Expression)"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),  # Remove x-axis text (tick labels)
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(),
    legend.position = "none"  # Remove legend
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_data$activity, na.rm = TRUE),
    label = ifelse(is.na(p_value), "p-value = NA", paste("p-value =", round(p_value, 4))),
    size = 4,
    vjust = -0.5
  )

p




#########################################################################################################################################
# create 2 rank files using z-scores and t-value for GSEA

# Define your groups
short_term_survivors <- c(18, 12, 14, 7, 10)
long_term_survivors <- c(5, 20, 21, 19, 2, 3, 13)

# Construct the column names from the patient IDs
short_term_cols <- paste0("patient_", short_term_survivors, "_log2FC")
long_term_cols <- paste0("patient_", long_term_survivors, "_log2FC")

# 1. Read the file
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/scaled_normalized_log2FC_C1_vs_Pre.csv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

# Ensure that all required columns exist in the data
missing_st <- setdiff(short_term_cols, colnames(data))
missing_lt <- setdiff(long_term_cols, colnames(data))

if (length(missing_st) > 0 || length(missing_lt) > 0) {
  stop("Some patient columns are missing in the data: ",
       paste(c(missing_st, missing_lt), collapse = ", "))
}

# 2. Calculate z-score for each gene
z_scores <- apply(data, 1, function(row) {
  # Extract group values
  short_values <- as.numeric(row[short_term_cols])
  long_values <- as.numeric(row[long_term_cols])
  
  # Calculate means
  mean_short <- mean(short_values, na.rm = TRUE)
  mean_long <- mean(long_values, na.rm = TRUE)
  
  # Calculate standard deviations and sample sizes
  sd_short <- sd(short_values, na.rm = TRUE)
  sd_long <- sd(long_values, na.rm = TRUE)
  
  n_short <- sum(!is.na(short_values))
  n_long <- sum(!is.na(long_values))
  
  # Calculate the standard error of the difference
  se_diff <- sqrt((sd_short^2 / n_short) + (sd_long^2 / n_long))
  
  # If standard error is zero, handle gracefully
  if (se_diff == 0) {
    z_val <- NA
  } else {
    # Reverse sign so that Z is positive if higher in long-term survivors
    z_val <- (mean_long - mean_short) / se_diff
  }
  
  return(z_val)
})

data$Zscore <- z_scores

# 3. Calculate t-values using a two-sample t-test for each gene
# We'll perform a two-sample t-test comparing long_values vs short_values
# and extract the t-statistic. By default: t.test(long, short) will test mean(long)=mean(short).
# The t will be positive if mean(long) > mean(short).
t_values <- apply(data, 1, function(row) {
  short_values <- as.numeric(row[short_term_cols])
  long_values <- as.numeric(row[long_term_cols])
  
  # Perform two-sample t-test assuming unequal variances (Welch's t-test)
  # If you want a more traditional approach assuming equal variance, set var.equal=TRUE
  t_res <- t.test(long_values, short_values, var.equal = FALSE)
  
  # Extract the t-value
  t_val <- t_res$statistic
  return(t_val)
})

data$Tvalue <- t_values

# 4. Create two rank files:

# Remove NA values from Z-score file
rank_file_z <- data.frame(Gene = data$Gene, Score = data$Zscore, stringsAsFactors = FALSE)
# Filter out rows where Score is NA
rank_file_z <- rank_file_z[!is.na(rank_file_z$Score), ]
# Sort by Score descending
rank_file_z <- rank_file_z[order(rank_file_z$Score, decreasing = TRUE), ]
write.table(rank_file_z, 
            file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/GSEA/gene_zscore_rankfile.rnk", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# Rank file with T-values
# Remove NA values from T-value file
rank_file_t <- data.frame(Gene = data$Gene, Score = data$Tvalue, stringsAsFactors = FALSE)
# Filter out rows where Score is NA
rank_file_t <- rank_file_t[!is.na(rank_file_t$Score), ]
# Sort by Score descending
rank_file_t <- rank_file_t[order(rank_file_t$Score, decreasing = TRUE), ]
write.table(rank_file_t, 
            file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/GSEA/gene_tvalue_rankfile.rnk", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)






#################################################################################################################################################
# survival analysis with pathway activity 

# Load necessary libraries
library(tidyverse)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(Rcpp)
library(cowplot)
library(stringr)  # For string manipulation

# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Define the Cox proportional hazards model function
get_cox = function(cox_input_subset_df, time_column, event_column, covariates_list)
{
  df = cox_input_subset_df
  df =  df[, c(time_column, event_column, covariates_list)]
  
  res.cox <- NULL
  
  formula_string = paste0("Surv(", time_column, ", ", event_column, ") ~ ", paste(covariates_list, collapse = " + "))
  new_formula <- as.formula(formula_string)
  
  tryCatch(
    {
      res.cox <- coxph(new_formula, data = df, control = coxph.control(iter.max=20))
    },
    error = function(e_outer) {
      cat("Outer Error: ", conditionMessage(e_outer), "\n")
    },
    warning = function(w_outer) {
      if (grepl("Ran out of iterations and did not converge", conditionMessage(w_outer))) {
        cat("Caught the specific warning: Ran out of iterations and did not converge\n")
        tryCatch(
          {
            res.cox <- coxph(new_formula, data = df, control = coxph.control(iter.max=1000))
          },
          error = function(e_inner) {
            cat("Inner Error: ", conditionMessage(e_inner), "\n")
          },
          warning = function(w_inner) {
            warning(w_inner)
          }
        )
      } else {
        # Handle other warnings if needed
        warning(w_outer)
      }
    }
  )
  
  return(list(cox = res.cox))
}

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/scaled_normalized_log2FC_C1_vs_Pre.csv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

# gene_list_200 <- c(
#   "C1QC",
#   "C1QB",
#   "LYPD2",
#   "S100A12",
#   "VSIG4",
#   "VCAN",
#   "SLCO5A1",
#   "FOS",
#   "S100A8",
#   "C1QA",
#   "CD163",
#   "ADAMTS2",
#   "HFM1",
#   "VMO1",
#   "CEBPD",
#   "GPX1",
#   "RBP7",
#   "AREG",
#   "THBS1",
#   "S100A9",
#   "GNLY",
#   "PLD4",
#   "IL32",
#   "CYP4F22",
#   "GZMB",
#   "CD247",
#   "HES4",
#   "TPST1",
#   "MS4A6A",
#   "PPBP",
#   "G0S2",
#   "CD14",
#   "SCGB3A1",
#   "TEX14",
#   "ITGA2B",
#   "FMN1",
#   "FLT3",
#   "MYL9",
#   "IGHV1-46",
#   "NKG7",
#   "ACSL1",
#   "FKBP5",
#   "H3C10",
#   "NRGN",
#   "CLU",
#   "TFEC",
#   "TSPAN33",
#   "SASH1",
#   "CAVIN2",
#   "TREML1",
#   "RNASE2",
#   "LPL",
#   "SPON2",
#   "CKB",
#   "GZMH",
#   "PRF1",
#   "DDIT4",
#   "HSH2D",
#   "CDKN1C",
#   "HBB",
#   "SH2D1B",
#   "DOCK4",
#   "EGR1",
#   "F13A1",
#   "IGLV2-14",
#   "CTSW",
#   "RBFOX2",
#   "CST7",
#   "PLCB1",
#   "IGLV5-45",
#   "SLAMF7",
#   "LMNA",
#   "GZMA",
#   "SELPLG",
#   "DUSP5",
#   "SYTL1",
#   "FGFBP2",
#   "MPIG6B",
#   "PRKN",
#   "PF4",
#   "SHTN1",
#   "ARL4A",
#   "PRKAR2B",
#   "TSTD1",
#   "CMTM5",
#   "SPARC",
#   "CRIP1",
#   "LYZ",
#   "STAT4",
#   "ARID5B",
#   "SPN",
#   "IRS2",
#   "IL7R",
#   "NRG1",
#   "STAB1",
#   "APOO",
#   "ADGRB3",
#   "EIF4G3",
#   "BLVRB",
#   "MTURN",
#   "KLRB1",
#   "ASGR2",
#   "CHST2",
#   "BAIAP2L1",
#   "PGRMC1",
#   "IGKV2-30",
#   "MT-ATP8",
#   "ETS1",
#   "MSR1",
#   "OAS1",
#   "CSF3R",
#   "SLC44A2",
#   "C2orf88",
#   "BEX3",
#   "TRBC1",
#   "NFKBIA",
#   "SETBP1",
#   "HBG1",
#   "SLC2A6",
#   "RNASE6",
#   "CUX1",
#   "IL13RA1",
#   "SAMSN1",
#   "IGLV3-27",
#   "SMOX",
#   "FRMD4B",
#   "ETS2",
#   "CD36",
#   "CD99",
#   "IGLV1-51",
#   "TRIM58",
#   "PID1",
#   "SH3BGRL2",
#   "DUSP1",
#   "CD7",
#   "PIK3CG",
#   "HLA-DMB",
#   "H2AC6",
#   "GBP2",
#   "LGALS2",
#   "CPED1",
#   "TRBC2",
#   "FGR",
#   "MCEMP1",
#   "SIGLEC10",
#   "IGKV1-27",
#   "IL1R2",
#   "TMEM91",
#   "CTSC",
#   "SCART1",
#   "MGST1",
#   "GPR183",
#   "ITGAM",
#   "CCL4",
#   "ENTPD1",
#   "IL21R",
#   "TRBV29-1",
#   "MLXIP",
#   "KLF12",
#   "RHOC",
#   "EPB41L3",
#   "IGHV4-59",
#   "LINGO2",
#   "CCR7",
#   "SMIM7",
#   "LAT2",
#   "TNFAIP3",
#   "P2RY10",
#   "SSBP4",
#   "MMD",
#   "EPHX3",
#   "HCAR3",
#   "FAM110A",
#   "TMEM40",
#   "NAMPT",
#   "ENSG00000285589",
#   "TRGV9",
#   "CD101",
#   "HBEGF",
#   "CTSL",
#   "C1orf167",
#   "PER1",
#   "LTB",
#   "ZNF467",
#   "SAP30",
#   "RAB39B",
#   "CD9",
#   "CBLB",
#   "CDH23",
#   "AHR",
#   "PKN1",
#   "KLF10",
#   "SYAP1",
#   "JUN",
#   "PIGQ",
#   "VAMP8",
#   "GNG11",
#   "S1PR4",
#   "ADAM8",
#   "PTPRCAP",
#   "RETN",
#   "CREB5",
#   "WARS1",
#   "GP9",
#   "SH2D3C",
#   "TUBB1",
#   "AKT3",
#   "CALHM2"
# )


# Define the gene list for the pathway
# gene_list_1100 <- c("ENSG00000269242", "IGLV1-40", "COL19A1", "GLB1L", "ENSG00000285542", "PKP2", "OIP5",
#                "MATR3", "VSIG4", "GAREM1", "CD163", "RNF227", "ALOX12", "TRAT1", "CD5", "CXCL3", "PKHD1L1",
#                "TAL1", "HBA1", "CXCL2", "PPARG", "EGR1", "SOCS3", "TEX14", "ANO5", "SIGLEC11", "CD163L1",
#                "ABTB2", "STOML1", "GPR153", "KLHDC7B", "LIME1", "H1-5", "RNF144A", "PARPBP", "TRAC",
#                "MAL", "ZFPM2", "IGF2BP1", "TPX2", "FGFBP2", "ANK2", "WDR62", "LHFPL5", "SDK1", "H3C2",
#                "CDC14B", "TRBC2", "STIL", "PRKG1", "H2AC11", "DYNC2H1", "ZNF532", "CENPF", "ENAH", "TAF9B",
#                "LYPD2", "ADGRB3", "SLIT2", "RAB27B", "BACH2", "ITGB7", "KLRC1", "MDGA2", "RGS7", "FCER2",
#                "H2BC15", "ENSG00000258311", "ELOVL7", "MAP7", "C1orf198", "MCM6", "ASAP2", "ST8SIA6",
#                "GPR174", "OXCT1", "HCAR3", "TIGIT", "C1QC", "LDHD", "VMO1", "PAX8", "TSHZ2", "FOS", "GPR183",
#                "NR4A2", "JUN", "PLEKHF1", "FHL1", "BMP6", "H2BC8", "SH3BGRL2", "SCN1B", "PCSK6", "CABP5",
#                "CXCL5", "BEND2", "ZNF417", "MFAP3L", "S1PR5", "HOMER2", "P2RY12", "NT5M", "CTTN", "PRKAR1B",
#                "THBS1", "SNCA", "PROS1", "ESAM", "LGALSL", "CTDSPL", "PDGFA", "PTCRA", "ITGA2B", "H2BC9",
#                "DAB2", "TSPAN33", "TRIM58", "PPBP", "GP9", "ITGB3", "C19orf33", "PF4", "SPARC", "F13A1",
#                "CAVIN2", "ACRBP", "MYL9", "CLU", "TSC22D1", "STON2", "SLC40A1", "LTBP1", "H2BC11", "FIBCD1",
#                "ABLIM3", "PF4V1", "TREML1", "TMEM40", "MMD", "GNG11", "TUBB1", "CMTM5", "MPIG6B", "RHOBTB1",
#                "IL21R", "DMTN", "PTGS1", "SAV1", "PANX1", "H4C8", "PLA2G12A", "PDLIM1", "H3C10", "BEX3",
#                "SLC12A4", "BACE1", "PLTP", "GATD3A", "MTG2", "LONRF3", "UBXN8", "LMNA", "ZNF358", "ACBD4",
#                "ENSG00000271741", "TNFRSF10D", "SIK1B", "NUDT6", "TNFAIP1", "CCL4", "ENSG00000205045", "NUP35",
#                "CAMKK1", "ENSG00000285827", "FLT3", "MARVELD1", "RGL1", "PER1", "FMN1", "IL1R2", "ZBTB16",
#                "DDIT4", "FKBP5", "DGKH", "NAA80", "BIN1", "LEF1", "KIF5C", "CATSPER1", "HOXB2", "SH2D1A",
#                "C2orf88", "YPEL1", "ZNF707", "NUDCD1", "RPH3A", "CACNB4", "SACS", "PPOX", "ADCY3", "AMOT",
#                "NBEA", "CACNA2D3", "GGACT", "RNLS", "KCTD18", "RAB13", "OPN3", "CHD7", "OTUD3", "HSPA1B",
#                "NOXA1", "ZBTB49", "PTGDR2", "TIGD1", "RETREG1", "EPHB2", "CTTNBP2", "SLC35C1", "HEMK1",
#                "SLC2A9", "ATP8B4", "CLPB", "MAFG", "PIP4P2", "MRPS17", "CYP27A1", "MROH6", "MEI1", "C17orf80",
#                "PDRG1", "ZNF432", "CNNM2", "CRISPLD2", "FOLR3", "BLZF1", "PTK2", "MVB12B", "SS18L1", "FASTKD1",
#                "ENSG00000260729", "CLECL1", "NRG1", "ACTR3B", "TMEM150A", "CACNA1A", "ASGR2", "SLC38A7",
#                "GPRC5C", "RAB39A", "GLT1D1", "CES1", "CALCRL", "HOMER3", "SESN2", "CCR2", "MAF", "ZP3",
#                "ZNF354A", "TP53I11", "TRAF1", "CD247", "PAQR7", "FCRL6", "ARL4A", "SAP30", "PRKCQ", "SPON2",
#                "CPM", "FPR2", "ENSG00000275464", "ZNF331", "PADI4", "TAF4B", "ODAD3", "ZNF410", "PTMS",
#                "MMP17", "CLEC10A", "DYSF", "NOL12", "CAMK4", "LRRC57", "BCL11B", "COQ8A", "MAP2K6", "RXRB",
#                "ICAM1", "PRR5", "PPIF", "MICAL2", "ENSG00000260272", "EVA1B", "AP1S1", "SATB1", "RETN",
#                "ATP6V0A1", "MT1X", "FAM20C", "RFX2", "ATF5", "CPED1", "RELB", "COIL", "TOR1B", "OSBPL5",
#                "METTL8", "KIF13A", "MSR1", "PLCB1", "RNASE6", "IL13RA1", "GNA15", "SIRPA", "PID1", "GSTM1",
#                "CSNK2B", "SNRPD3", "EEF1G", "RNASE2", "TFEC", "CD93", "VGLL4", "MS4A6A", "RBP7", "TBC1D2",
#                "ALDH2", "LGALS2", "ASGR1", "ZNF467", "ITGAM", "FOLR2", "ADAM15", "DOCK4", "IL1RN", "C1QB",
#                "C1QA", "NCF1", "TCN2", "FUOM", "HIP1", "PLBD1", "MARCO", "FES", "NAIP", "TMEM91", "SLC35D2",
#                "IRS2", "MLXIP", "HBEGF", "SLAMF7", "OAS3", "EPSTI1", "GBP5", "MEN1", "EHD4", "AMDHD2",
#                "FCF1", "RILP", "TRIM69", "TMEM106A", "FCGR1A", "ENSG00000124593", "AGAP3", "STAB1", "PTAFR",
#                "PXN", "ZNF101", "SNX25", "TMEM128", "LYSMD3", "RRP15", "RMND1", "POLR1G", "TMEM42", "PDE7B",
#                "C3", "DHX58", "ZNF224", "SPIRE1", "PPM1N", "DTD1", "PCGF1", "ICA1L", "ST7L", "SAMD3", "ETS1",
#                "DNHD1", "PLD4", "INKA2", "HSH2D", "TAGAP", "CTSL", "CYP4F22", "TRPT1", "PTPN4", "KIAA0586",
#                "RAB37", "NKG7", "LBH", "UPK3A", "PCBP3", "PRKN", "ZC3H10", "ABCA7", "CBR4", "CASP3", "NELFA",
#                "TMEM104", "TXN", "RARS2", "EI24", "TMEM80", "GLT8D1", "CBLB", "ZRANB3", "INSR", "DUSP1",
#                "NFKBIA", "ZFP36", "ZMYND11", "AHR", "ISG20", "RHBDD2", "CYSTM1", "PCED1B", "AQP9", "AKT3",
#                "PAPSS2", "TTC39B", "RPL41", "IFITM1", "TRMT61A", "CD101", "TCEAL4", "ST3GAL6", "SNED1", "NID1",
#                "MLLT1", "ASB6", "SZT2", "CAPN10", "TRAF5", "NRGN", "PRKAR2B", "SKAP1", "KIF2A", "SYNE1",
#                "MAP3K7CL", "C10orf143", "H2BC12", "PDCL", "C5orf22", "NUDT18", "LSS", "SLC35E3", "GFOD1",
#                "MTPAP", "PEX5", "SMG9", "CRIPT", "SMARCD3", "U2AF1L4", "ZBED5", "NPIPB4", "ANKAR", "DDB2",
#                "CISD1", "ANKZF1", "H2BC4", "RPS6KA2", "DLG4", "PLEKHA7", "MGST1", "HPF1", "ACTN1", "LTB",
#                "TMEM102", "TSTD1", "ZDHHC2", "SH2D1B", "PGRMC1", "H2AC6", "SPHK2", "ASPH", "CAMTA1", "NBPF12",
#                "ISG15", "HEG1", "SLC25A20", "ARID5B", "KLHL22", "SLFN5", "IFI6", "ZC2HC1A", "TBC1D15",
#                "L3MBTL3", "ITPKB", "ERCC5", "H1-0", "KANSL1L", "DUS2", "NAA15", "CYSLTR1", "APP", "SACM1L",
#                "MT-ATP8", "FLT3LG", "TMEM154", "TNFRSF8", "ACAP1", "ARL6IP6", "ANKRD36", "PLCL1", "DALRD3",
#                "RBM43", "TSEN15", "PLAGL1", "SAMSN1", "ERGIC2", "LMBRD1", "RHOF", "RPGR", "RNF146", "GPR160",
#                "NUDT15", "CEP85", "TTC9", "FAM200B", "ZMAT3", "C16orf54", "SYAP1", "SNRNP35", "ELL2", "MYBBP1A",
#                "SHTN1", "DBNDD2", "GPAT4", "ABHD16A", "CHST2", "DYRK2", "PTP4A3", "PRR5L", "CX3CR1", "DUSP5",
#                "IVD", "SETBP1", "CSF3R", "MKKS", "TRAPPC5", "PLD3", "RPS17", "S100A12", "S100A8", "MLEC", "IPO5",
#                "CARD19", "MTMR3", "CD14", "SASH1", "VCAN", "FRMD4B", "CRTAP", "FCN1", "LYZ", "CFP", "ARHGEF10L",
#                "S100A9", "RBM47", "GPX1", "C1RL", "MERTK", "DCAF12", "SLC35F6", "ETHE1", "JARID2", "ARHGAP31",
#                "PCBD1", "ARHGAP24", "CPNE8", "SPRED1", "CLIC4", "TAF2", "ATP2A2", "NUP210", "NR6A1", "PEX16",
#                "BRD3", "ZMYM4", "SNTB1", "PARVB", "UBASH3B", "CREB5", "TBC1D9", "CBX6", "TNK2", "GOLIM4",
#                "TPM1", "CD36", "TRPS1", "LIMK2", "TNS3", "MCFD2", "SIRT6", "FBN2", "ZNF33A", "TST", "MXD1",
#                "SLC2A3", "DGKG", "MTMR11", "TACC3", "IL6ST", "SMIM7", "TBCD", "EIF4G3", "XPO6", "HLCS",
#                "PITPNB", "NFU1", "HAGH", "SNAPC3", "RNF24", "BICRAL", "RASA1", "ITGAE", "SH3GLB2", "PLA2G4A",
#                "ENTPD1", "ANPEP", "CALML4", "MTMR10", "PMM2", "CAMTA2", "CD244", "S100Z", "SMAD1", "APOBEC3C",
#                "RBM22", "NFIL3", "AARSD1", "ACSL1", "TNFAIP3", "NAMPT", "KLF10", "B4GALT1", "EMILIN2", "LY86",
#                "NME2", "HLA-DMB", "VAMP8", "HLA-DMA", "C3AR1", "C18orf32", "CTSB", "TNFRSF14", "ODF3B", "PIGQ",
#                "MSL1", "TRIM22", "BLVRB", "LITAF", "LRP10", "SHISA5", "ADIPOR1", "ERAP1", "ATP1B3", "HSP90B1",
#                "UBC", "SERPINA1", "SARAF", "CIRBP", "ARL6IP5", "RPL29", "EEF1B2", "CSTA", "HLA-DRB5", "AP1S2",
#                "CCDC12", "UQCRH", "RPLP2", "RPLP1", "RPL13", "RPL18", "RPS24", "RPL30", "RPS8", "RPS18",
#                "RPS14", "RPL19", "RPL5", "RPL9", "RPL12", "RPL32", "RPL10", "RPL7A", "RPL18A",
#                "RPL6", "RPS5", "RPS3A", "RPL24", "RPLP0", "RPL15", "RPS6", "GAPDH", "RPL23", "ETS2", "COPG2",
#                "SESN1", "FMNL2", "XPR1", "CEBPD", "JDP2", "FOXO3", "GSN", "CROCC2", "POLM", "FUCA1", "TFEB",
#                "STK25", "SNX20", "NBPF26", "TMBIM4", "SARNP", "RAB3D", "PCSK7", "ADAT1", "ALDH16A1", "TAPBPL",
#                "SPSB3", "SLC66A3", "RBM4", "SNX18", "OAS2", "TMEM176B", "IFIT2", "KNDC1", "ADIPOR2", "PACS2",
#                "EHBP1", "PAICS", "GLO1", "TTC19", "TRUB2", "RABEP1", "TIMM50", "PCK2", "AGPAT3", "ARHGAP45",
#                "ITGAL", "LPP", "NCOA3", "RNASEK", "TNFRSF1B", "RHOC", "MYO1G", "LFNG", "CRIP1", "ARAP1",
#                "PKN1", "INSIG1", "KLF4", "CCM2", "LNPEP", "GMIP", "CD58", "KIAA2026", "SH2D3C", "SPN", "PIK3CG",
#                "PHTF2", "CD79B", "DRAP1", "USF3", "LYN", "SVIL", "RAP1GAP2", "ARAP2", "RAB29", "KLHL5", "CERK",
#                "MXD3", "HSPA5", "SUN2", "RUNX3", "RBM38", "S1PR4", "CDKN1C", "VASP", "TUBA1A", "CEBPA", "MCM5",
#                "GYPC", "HSPH1", "DNAJC10", "TRMT1", "PLAAT4", "MNDA", "CISD3", "CD52", "UCHL3", "METTL9",
#                "NDUFB11", "YWHAG", "XPC", "HSPE1", "RPL22L1", "GLRX", "RPL39", "MRPL52", "RPS28", "TOMM5",
#                "COMMD6", "IER3IP1", "C11orf21", "PFDN1", "ZDHHC1", "KLF12", "IPMK", "TBC1D10C", "APOBEC3G",
#                "LIPA", "KLF3", "TRAF3IP3", "DDT", "ESYT1", "ZFP36L1", "PPP1R14B", "UBE2D1", "SPIDR", "MS4A4E",
#                "MCL1", "PMVK", "DNAJB11", "NFE2L2", "EEA1", "SLC15A4", "FOXN3", "TMEM170B", "HSP90AA1", "QKI",
#                "DYNLL1", "PISD", "FNDC3B", "WDFY3", "CMIP", "PSMB5", "YTHDF1", "IGF2R", "APEX1", "DDX21",
#                "HEBP2", "HSPD1", "EPB41L3", "PTPRE", "COMT", "YBX3", "EMC6", "FERMT3", "AP2S1", "INPP4A",
#                "GRK2", "ADK", "WDR43", "SSBP2", "MPC1", "PTPN1", "ZMYND8", "LEPROTL1", "ID2", "STIM2", "IRF7",
#                "GBP4", "FGR", "TNFSF10", "GBP2", "LY6E", "WARS1", "MAD1L1", "SLC44A2", "PFKL", "TESC", "STIP1",
#                "CCNDBP1", "OST4", "DHRS7", "FAM110A", "KLF2", "CFD", "SYNGR2", "FCGR3A", "RAC2", "UNC119",
#                "CTSC", "OAS1", "CXCL16", "LAMTOR1", "SOD1", "TMC6", "UROS", "CDH23", "SETDB1", "ELOVL1",
#                "SF3B4", "PELI1", "CEBPB", "ARHGEF18", "DPEP2", "SLC31A2", "PML", "ZDHHC12", "SELPLG", "SSBP4",
#                "HLA-E", "UNC93B1", "SIGLEC10", "WDR1", "IL16", "LEMD2", "DEDD2", "CCNL1", "NUTF2", "NOSIP",
#                "LIMD2", "MAPKAPK3", "ADGRE2", "VSIR", "MTSS1", "C5AR2", "N4BP2L2", "CSTB", "TBCB", "PHF19",
#                "SP110", "SYTL1", "PECAM1", "CUX1", "HDGF", "NAP1L1", "PTPN6", "LYST", "GNG2", "FAM126A",
#                "RASAL3", "IPO7", "ZNF706", "STXBP2", "CANX", "HSP90AB1", "SPINT2", "ICAM3", "ADGRE3", "PPFIA1",
#                "NUDC", "TMEM120A", "CD99", "ST3GAL1", "INPP5D", "CTBP2", "FRY", "TPST2", "LPCAT3", "RNF11",
#                "ATP1A1", "DIP2A", "CYFIP2", "COQ2", "ANKRD44", "FGD3", "YWHAE", "MDFIC", "POLD4", "GPAA1",
#                "TXNDC17", "GOLGA4", "PRKX", "LAT2", "HLA-F", "TRAPPC12", "LIMD1", "TPM4", "ARHGEF1", "UCP2",
#                "RNF213", "PITPNC1", "LGALS9", "CD300LF", "ADGRE5", "UNC13D", "ITGAX", "MAP3K1", "CALHM2",
#                "EZR", "FLNA", "CSF1R", "GBP1", "JAK3", "LDB1", "TTYH3", "RFTN1", "TCIRG1", "SLC2A6", "VASH1",
#                "PDIA4", "HEXD", "UBE2S", "SLC35E1", "METTL22", "IRAK4", "GRAMD1A", "VPS26B", "DNAJC2", "CHMP1A",
#                "GALNT1", "TRAFD1", "SOX4", "RMC1", "ADA", "HES4", "CCDC115", "INPP5K", "CDIP1", "QPRT", "PRPF4",
#                "SLC8B1", "NIP7", "GRK5", "KLF9", "DHPS", "NDST1", "PLSCR3", "BCAP29", "SLC39A1", "RIPK3",
#                "ENSG00000285589", "ENSG00000263620", "GATD3B", "PRDM1", "CD38", "GMEB2", "HPSE", "RNF25",
#                "CD22", "UTP4", "DNAJB9", "AARS1", "SYT17", "PTGDR", "PGLYRP2", "GMPR", "TCEAL3", "SVIP",
#                "DNM3", "ARMCX5", "SMIM3", "CD9", "MTURN", "NKIRAS1", "GSTA4", "BRCA2", "SCARB1", "H4C5",
#                "ZNF627", "TRIM3", "MTHFD2L", "GSKIP", "HAUS6", "BLMH", "HAL", "TMEM106C", "SIRT5", "FCRL3",
#                "DHRS12", "STAT4", "B3GNT8", "RET", "NCR3", "CYP2U1", "BNIP3", "CISH", "ADAM8", "OPTN",
#                "NCALD", "TBX21", "P2RY10", "TOX", "TXK", "IL2RB", "PRF1", "GZMB", "SH2D2A", "CD96", "PTPRCAP",
#                "CD7", "NCKAP1", "GZMA", "CST7", "CTSW", "GNLY", "GZMM", "GZMH", "CLIC3", "AQP3", "KLRB1",
#                "TRBC1", "FEN1", "TGFBR3", "FAM110B", "RGL4", "NME4", "CKB", "SFTPD", "RSKR", "TBC1D19", "HHAT",
#                "FCMR", "SCML4", "GPC6", "MPP6", "S1PR1", "IL7R", "WDR27", "SPOCK2", "ZAP70", "GATA3", "MYBL1",
#                "ABLIM1", "LCK", "IL32", "CD3E", "TCF7", "TMOD1", "LPL", "GPER1", "DEPTOR", "APOO", "IGLV3-19",
#                "IGKV4-1", "TRBV4-2", "IL6R", "IL18", "IL27RA")


# 1) Character vector for the first gene list
gene_list_500 <- c(
  "C1QC","C1QB","LYPD2","S100A12","VSIG4","VCAN","SLCO5A1","FOS","S100A8","C1QA","CD163",
  "ADAMTS2","HFM1","VMO1","CEBPD","GPX1","RBP7","AREG","THBS1","S100A9","GNLY","PLD4",
  "IL32","CYP4F22","GZMB","CD247","HES4","TPST1","MS4A6A","PPBP","G0S2","CD14","SCGB3A1",
  "TEX14","ITGA2B","FMN1","FLT3","MYL9","IGHV1-46","NKG7","ACSL1","FKBP5","H3C10","NRGN",
  "CLU","TFEC","TSPAN33","SASH1","CAVIN2","TREML1","RNASE2","LPL","SPON2","CKB","GZMH",
  "PRF1","DDIT4","HSH2D","CDKN1C","HBB","SH2D1B","DOCK4","EGR1","F13A1","IGLV2-14",
  "CTSW","RBFOX2","CST7","PLCB1","IGLV5-45","SLAMF7","LMNA","GZMA","SELPLG","DUSP5",
  "SYTL1","FGFBP2","MPIG6B","PRKN","PF4","SHTN1","ARL4A","PRKAR2B","TSTD1","CMTM5",
  "SPARC","CRIP1","LYZ","STAT4","ARID5B","SPN","IRS2","IL7R","NRG1","STAB1","APOO",
  "ADGRB3","EIF4G3","BLVRB","MTURN","KLRB1","ASGR2","CHST2","BAIAP2L1","PGRMC1","IGKV2-30",
  "MT-ATP8","ETS1","MSR1","OAS1","CSF3R","SLC44A2","C2orf88","BEX3","TRBC1","NFKBIA",
  "SETBP1","HBG1","SLC2A6","RNASE6","CUX1","IL13RA1","SAMSN1","IGLV3-27","SMOX","FRMD4B",
  "ETS2","CD36","CD99","IGLV1-51","TRIM58","PID1","SH3BGRL2","DUSP1","CD7","PIK3CG",
  "HLA-DMB","H2AC6","GBP2","LGALS2","CPED1","TRBC2","FGR","MCEMP1","SIGLEC10","IGKV1-27",
  "IL1R2","TMEM91","CTSC","SCART1","MGST1","GPR183","ITGAM","CCL4","ENTPD1","IL21R",
  "TRBV29-1","MLXIP","KLF12","RHOC","EPB41L3","IGHV4-59","LINGO2","CCR7","SMIM7","LAT2",
  "TNFAIP3","P2RY10","SSBP4","MMD","EPHX3","HCAR3","FAM110A","TMEM40","NAMPT",
  "ENSG00000285589","TRGV9","CD101","HBEGF","CTSL","C1orf167","PER1","LTB","ZNF467","SAP30",
  "RAB39B","CD9","CBLB","CDH23","AHR","PKN1","KLF10","SYAP1","JUN","PIGQ","VAMP8","GNG11",
  "S1PR4","ADAM8","PTPRCAP","RETN","CREB5","WARS1","GP9","SH2D3C","TUBB1","AKT3","CALHM2",
  "LY6E","TSHZ2","MXD1","RGL4","ITGB3","PTK2","ENSG00000259529","APOBEC3G","H2BC11","ACRBP",
  "SIRPA","TMEM128","PRKCQ","SPSB3","UPK3A","KNDC1","TSC22D1","RPLP0","NXF3","ADORA2A",
  "IFITM1","NCF1","STON2","ZAP70","UNC119","EEF1G","PDGFA","SH3BP4","GFOD1","IFIT2",
  "SNED1","CEBPA","S100P","FN1","TMEM176B","LPCAT3","IGHV2-70D","ITGAL","TTC39B",
  "IGLV4-60","RPL41","CLIC3","SNCA","LCK","CD300LF","CCR2","HOMER3","LY9","CSF1R","LGALSL",
  "LAPTM4B","POLD4","PTGDR2","GATD3B","ANPEP","SYNGR2","SLIT2","AP1S2","NT5M","TRBV6-6",
  "PPM1N","ARAP2","LFNG","LMBRD1","GBP5","TTC9","DRAP1","TP53I11","ADGRE2","ZDHHC1","RHOF",
  "MAL","COMT","TCIRG1","TRBV12-3","MYBL1","IGKV3D-11","ALOX12","GPRC5C","CD3E","ARHGAP24",
  "PPIF","TMEM120A","PTGS1","LDHD","CX3CR1","MEN1","DHPS","ZC2HC1A","TBC1D10C","TRDC",
  "ZNF532","STXBP2","ACAP1","TESC","IGLV2-5","YWHAG","PTPRE","CD244","S1PR5","SLC2A9",
  "UCP2","RNASE1","ODAD3","APOE","MCM6","H3C2","PEX5","PADI4","PLA2G4A","CAMTA1",
  "PLA2G12A","CLIC4","PLEKHA7","MARVELD1","SAV1","OBSCN","SIGLEC11","SPIDR","SLC38A7",
  "IGKV4-1","C1RL","ERAP1","ARL6IP5","TBCD","DNM3","CAV2","SLC22A16","IFI6",
  "ENSG00000124593","GBP4","JDP2","CFD","IGF2BP1","RPLP1","PECAM1","PRKG1","ANO5","AMOT",
  "CRISPLD2","IGHV5-51","H4C8","SNX18","MAP3K1","ID4","MERTK","SLC15A4","MS4A4E","HOXB2",
  "AARS1","LAIR2","TRBV25-1","GYPC","NAIP","H2BC4","TSPAN12","RASAL3","MICAL2","GPR174",
  "PTPN6","NUDCD1","TRBV2","COL26A1","HLA-DMA","ACTN1","BICRAL","MFSD2B","ELOVL7","TFPI",
  "CAMK4","H2AC11","NR4A3","NME4","PAQR7","PHF19","FCER1A","CSTA","IGLV1-47","LEF1","MCFD2",
  "UNC93B1","DGKH","DDT","RETREG1","CLPB","IGHV4-31","OST4","TRIM3","CD93","SLC40A1",
  "ARMCX5","ANO4","NPIPB4","C10orf143","HECTD2","SAPCD2","NCOA3","PFDN1","IER3IP1",
  "SLC2A3","CES1","LIME1","CXCL16","SHCBP1","TRBV18","MAFIP","LYN","YBX3","ZNF627","CMIP",
  "ZNF331","TRBV7-6","RPS8","RNF11","NBPF26","GALNT16","FCN1","TRPS1","CCDC9B","PLBD1",
  "MAD1L1","PF4V1","ENSG00000280571","CDIP1","TRAF3IP3","CPNE8","SFTPD","FLT3LG","DYSF",
  "TNFSF10","SYT17","MYO7A","SHISA5","CASP3","POLM","PRPF4","GPC6","ABHD16A","DUS2",
  "IGHV1-3","CDC45","SMARCD3","CAMKK1","PCK2","HSP90B1","TMC6","PLSCR3","PLCL1","MTSS1",
  "PDLIM1","RPS13","SERPINE1","H3C4","LIMD1","KCNB2","MNDA","ZFP36L1","PRR5L","OAS2","SKAP1",
  "CRTAP","MEIS1","CARD19","CERK","CCNDBP1","BIN1","TRMT1","PITPNB","RAB9B","PTPN4","ENAH",
  "RUNX3","LEPROTL1","NR6A1","HBA2","RAB13","PTMS","ADGRE5","RNF24","ARMC12","MKI67",
  "EPB41L1","SUN2","DBNDD2","ENSG00000272442","TRBV5-5","GSTM1","DPEP2","TNS3","SCML4",
  "STK25","GZMM"
)

# 2) Character vector for the second (longer) gene list
# gene_list_1000 <- c(
#   "C1QC","C1QB","LYPD2","S100A12","VSIG4","VCAN","SLCO5A1","FOS","S100A8","C1QA","CD163",
#   "ADAMTS2","HFM1","VMO1","CEBPD","GPX1","RBP7","AREG","THBS1","S100A9","GNLY","PLD4",
#   "IL32","CYP4F22","GZMB","CD247","HES4","TPST1","MS4A6A","PPBP","G0S2","CD14","SCGB3A1",
#   "TEX14","ITGA2B","FMN1","FLT3","MYL9","IGHV1-46","NKG7","ACSL1","FKBP5","H3C10","NRGN",
#   "CLU","TFEC","TSPAN33","SASH1","CAVIN2","TREML1","RNASE2","LPL","SPON2","CKB","GZMH",
#   "PRF1","DDIT4","HSH2D","CDKN1C","HBB","SH2D1B","DOCK4","EGR1","F13A1","IGLV2-14",
#   "CTSW","RBFOX2","CST7","PLCB1","IGLV5-45","SLAMF7","LMNA","GZMA","SELPLG","DUSP5",
#   "SYTL1","FGFBP2","MPIG6B","PRKN","PF4","SHTN1","ARL4A","PRKAR2B","TSTD1","CMTM5",
#   "SPARC","CRIP1","LYZ","STAT4","ARID5B","SPN","IRS2","IL7R","NRG1","STAB1","APOO",
#   "ADGRB3","EIF4G3","BLVRB","MTURN","KLRB1","ASGR2","CHST2","BAIAP2L1","PGRMC1","IGKV2-30",
#   "MT-ATP8","ETS1","MSR1","OAS1","CSF3R","SLC44A2","C2orf88","BEX3","TRBC1","NFKBIA",
#   "SETBP1","HBG1","SLC2A6","RNASE6","CUX1","IL13RA1","SAMSN1","IGLV3-27","SMOX","FRMD4B",
#   "ETS2","CD36","CD99","IGLV1-51","TRIM58","PID1","SH3BGRL2","DUSP1","CD7","PIK3CG",
#   "HLA-DMB","H2AC6","GBP2","LGALS2","CPED1","TRBC2","FGR","MCEMP1","SIGLEC10","IGKV1-27",
#   "IL1R2","TMEM91","CTSC","SCART1","MGST1","GPR183","ITGAM","CCL4","ENTPD1","IL21R",
#   "TRBV29-1","MLXIP","KLF12","RHOC","EPB41L3","IGHV4-59","LINGO2","CCR7","SMIM7","LAT2",
#   "TNFAIP3","P2RY10","SSBP4","MMD","EPHX3","HCAR3","FAM110A","TMEM40","NAMPT",
#   "ENSG00000285589","TRGV9","CD101","HBEGF","CTSL","C1orf167","PER1","LTB","ZNF467","SAP30",
#   "RAB39B","CD9","CBLB","CDH23","AHR","PKN1","KLF10","SYAP1","JUN","PIGQ","VAMP8","GNG11",
#   "S1PR4","ADAM8","PTPRCAP","RETN","CREB5","WARS1","GP9","SH2D3C","TUBB1","AKT3","CALHM2",
#   "LY6E","TSHZ2","MXD1","RGL4","ITGB3","PTK2","ENSG00000259529","APOBEC3G","H2BC11","ACRBP",
#   "SIRPA","TMEM128","PRKCQ","SPSB3","UPK3A","KNDC1","TSC22D1","RPLP0","NXF3","ADORA2A",
#   "IFITM1","NCF1","STON2","ZAP70","UNC119","EEF1G","PDGFA","SH3BP4","GFOD1","IFIT2",
#   "SNED1","CEBPA","S100P","FN1","TMEM176B","LPCAT3","IGHV2-70D","ITGAL","TTC39B",
#   "IGLV4-60","RPL41","CLIC3","SNCA","LCK","CD300LF","CCR2","HOMER3","LY9","CSF1R","LGALSL",
#   "LAPTM4B","POLD4","PTGDR2","GATD3B","ANPEP","SYNGR2","SLIT2","AP1S2","NT5M","TRBV6-6",
#   "PPM1N","ARAP2","LFNG","LMBRD1","GBP5","TTC9","DRAP1","TP53I11","ADGRE2","ZDHHC1","RHOF",
#   "MAL","COMT","TCIRG1","TRBV12-3","MYBL1","IGKV3D-11","ALOX12","GPRC5C","CD3E","ARHGAP24",
#   "PPIF","TMEM120A","PTGS1","LDHD","CX3CR1","MEN1","DHPS","ZC2HC1A","TBC1D10C","TRDC",
#   "ZNF532","STXBP2","ACAP1","TESC","IGLV2-5","YWHAG","PTPRE","CD244","S1PR5","SLC2A9",
#   "UCP2","RNASE1","ODAD3","APOE","MCM6","H3C2","PEX5","PADI4","PLA2G4A","CAMTA1",
#   "PLA2G12A","CLIC4","PLEKHA7","MARVELD1","SAV1","OBSCN","SIGLEC11","SPIDR","SLC38A7",
#   "IGKV4-1","C1RL","ERAP1","ARL6IP5","TBCD","DNM3","CAV2","SLC22A16","IFI6",
#   "ENSG00000124593","GBP4","JDP2","CFD","IGF2BP1","RPLP1","PECAM1","PRKG1","ANO5","AMOT",
#   "CRISPLD2","IGHV5-51","H4C8","SNX18","MAP3K1","ID4","MERTK","SLC15A4","MS4A4E","HOXB2",
#   "AARS1","LAIR2","TRBV25-1","GYPC","NAIP","H2BC4","TSPAN12","RASAL3","MICAL2","GPR174",
#   "PTPN6","NUDCD1","TRBV2","COL26A1","HLA-DMA","ACTN1","BICRAL","MFSD2B","ELOVL7","TFPI",
#   "CAMK4","H2AC11","NR4A3","NME4","PAQR7","PHF19","FCER1A","CSTA","IGLV1-47","LEF1","MCFD2",
#   "UNC93B1","DGKH","DDT","RETREG1","CLPB","IGHV4-31","OST4","TRIM3","CD93","SLC40A1",
#   "ARMCX5","ANO4","NPIPB4","C10orf143","HECTD2","SAPCD2","NCOA3","PFDN1","IER3IP1",
#   "SLC2A3","CES1","LIME1","CXCL16","SHCBP1","TRBV18","MAFIP","LYN","YBX3","ZNF627","CMIP",
#   "ZNF331","TRBV7-6","RPS8","RNF11","NBPF26","GALNT16","FCN1","TRPS1","CCDC9B","PLBD1",
#   "MAD1L1","PF4V1","ENSG00000280571","CDIP1","TRAF3IP3","CPNE8","SFTPD","FLT3LG","DYSF",
#   "TNFSF10","SYT17","MYO7A","SHISA5","CASP3","POLM","PRPF4","GPC6","ABHD16A","DUS2",
#   "IGHV1-3","CDC45","SMARCD3","CAMKK1","PCK2","HSP90B1","TMC6","PLSCR3","PLCL1","MTSS1",
#   "PDLIM1","RPS13","SERPINE1","H3C4","LIMD1","KCNB2","MNDA","ZFP36L1","PRR5L","OAS2","SKAP1",
#   "CRTAP","MEIS1","CARD19","CERK","CCNDBP1","BIN1","TRMT1","PITPNB","RAB9B","PTPN4","ENAH",
#   "RUNX3","LEPROTL1","NR6A1","HBA2","RAB13","PTMS","ADGRE5","RNF24","ARMC12","MKI67",
#   "EPB41L1","SUN2","DBNDD2","ENSG00000272442","TRBV5-5","GSTM1","DPEP2","TNS3","SCML4",
#   "STK25","GZMM","TRAV13-2","AGAP3","ECT2L","ADAM15","CHMP1A","MN1","HSD17B6","POPDC2",
#   "SPIRE1","ALDH2","MXD3","IGLV7-43","RPL18A","IGKV3-7","TRIM69","ZC3HAV1L","HAGH","PLAGL1",
#   "LYPD3","DNAJB9","BLZF1","C10orf95","ITGAX","TACC3","UNC13D","NUDT18","OSBPL5","BNIP3",
#   "C16orf54","HOMER1","STIM2","CENPF","RAB27B","IPMK","RNASE3","TBC1D9","MYEOV","SLC39A1",
#   "TUB","RET","ELL2","PCSK6","TCEAL4","CCNA2","TRIM22","ARL6IP6","FCRL6","HLA-G","EFR3B",
#   "ZNF224","DDB2","HEBP2","TUBB2B","KLHDC7B","SPC24","SIRT5","B3GNT8","POLR1G","PARPBP",
#   "TRIM72","TNFRSF1B","MTMR11","IGHV3-21","ESYT1","RNASEK","ENSG00000269242","TBX21",
#   "NME1-NME2","RFTN1","FAM189A1","HSPH1","SCN1B","TBC1D2","NGFR","ATP8B4","CD96","MYO1G",
#   "LPP","HEMK1","CXCL3","MYBBP1A","WASF3","MAPKAPK3","NFE2L2","HHIPL1","C2orf16","WDFY3",
#   "LTBP1","APP","ENSG00000263620","FAM13C","MAP3K7CL","CALML4","ABCA13","RBM38","TUBA1A",
#   "RPL39","TAPBPL","NAA80","MROH6","TCN2","RPS18","ZNF101","PTAFR","RAP1GAP2","FUOM","PTPN1",
#   "PAPSS2","CYP27A1","INSR","RPS12","TINCR","HEG1","THAP10","ADCY3","GPAA1","P2RY12","VASP",
#   "SCARB1","SOCS3","U2AF1L4","UQCRH","TOR1B","FAM126A","MTHFD2L","ENSG00000255641","ICAM1",
#   "IQCA1","ZNF311","ITPKB","NOXA1","SYDE2","IGLV3-16","H2BC9","MFAP3L","C17orf80","HLA-F",
#   "SLC4A3","ZIK1","TMEM106A","DNHD1","PGLYRP2","CACNB4","CSTB","RPGR","TGFBR3","C3",
#   "TNFRSF10D","PRKAR1B","HSP90AA1","SLC31A2","VSIG1","PLEKHG6","HOMER2","ZDHHC12","RHOBTB1",
#   "ENSG00000204422","GRK2","NOL12","HTRA1","TPM4","GMIP","LRP2","ARHGAP31","CISH",
#   "ENSG00000285827","BACE1","H2BC12","LSS","KLF3","PCSK7","TAL1","H1-0","AARSD1","DEDD2",
#   "ZNF432","TMBIM4","PLAAT1","APBB1","NAA11","SLC25A27","SAMD3","NCR3","BRCA2","KLF9",
#   "NCAPH","N4BP2L2","RPH3A","SELENBP1","ASPH","SOCS6","SMIM3","ADGRE3","ISG20","ETHE1",
#   "PELI1","AQP9","RXRB","ANK2","TIAF1","KIF2A","PML","MRPL52","AP1S1","CATSPER1","H4C5",
#   "EHD4","LRRC57","PFKL","NCKAP1","DYNC2H1","ARHGEF1","DYNLL1","DGKG","TRAC","EZR",
#   "ENSG00000260272","IRF7","PTPRU","FOXO3","FAM200B","KLRK1","TAS2R46","CCM2","ESAM","FCER2",
#   "TCF7","RELB","RBM4","PREX2","KIAA0586","CD1E","GNG2","ADIPOR2","C5AR2","GPR156","LIMK2",
#   "RNF25","SEPTIN4","STUM","POU4F1","FCGR3A","COPG2","SACS","IL2RB","MYBPC2","EPHB2",
#   "ASAP2","DLG4","CNGA1","RPL15","SERPINB10","NDST1","FCRL3","MTG2","CD79B","SNTB1",
#   "AKR1E2","FMNL2","DMTN","MTPAP","FES","IL6ST","TRAF1","ICAM3","PLLP","IGDCC3","MLEC",
#   "GATD3A","MUSK","MAFG","ALDH16A1","GOLGA4","RPL29","RARS2","SVIP","CAPN10","KCTD4",
#   "PTP4A3","SMAD1","TXN","SERPINA1","CASR","PMVK","B4GALT1","LGMN","CDC14B","IGHV2-26",
#   "LY86","S100Z","MCM5","RUNX1T1","IGLV1-40","SPOCK2","TAF9B","CACNA2D3","ZNF33A","AGPAT3",
#   "TMEM106C","ANKZF1","NPAS2","ADK","GPAT4","CD22","PCED1B","INPP4A","RAI14","IVD","SLC35D2",
#   "TST","ERCC5","TMEM52B","RNF175","RPS6KA2","NUDT15","NCALD","LYSMD3","HPSE","ZMYND10",
#   "TRAV38-2DV8","IL1RN","PARVB","DNAJC2","SMTNL1","PPARG","CISD1","MARCO","STOML1","CANX",
#   "EMC6","CABP5","TRAPPC12","SLC25A20","NKIRAS1","ZNF358","CD5","NFIL3","TBCB","PPP1R3G",
#   "CXCR1","SVIL","RBM22","CLEC10A","NBPF12","GPER1","ZNF706","NR4A2","CTBP2","HSPD1",
#   "CATSPER3","F8A3","RPL32","OPTN","TRAV8-1","KIF5C","KLF2","METTL24","TBKBP1","PTCRA",
#   "GNA15","ENSG00000258311","ATP1A1","INKA2","NUP210","ADIPOR1","TPM1","FBN2","HSPE1",
#   "SLC35C1","MMP23B","IL16","TXK","GOLGA7B","HSPA12B","ID2","CYP3A7-CYP3A51P","ZFHX2",
#   "PXN","ATP2A2","ARHGEF18","C8orf37","PYGO1","FRY","CISD3","CSNK2B","ENSG00000276087",
#   "EPHX2","EI24","OPN3","ZBED5","ZNF707","PLD3","VASH1","EHBP1","FGFR2","RABEP1","C11orf21",
#   "RBPMS2","RBPMS","OXCT1","IGHV3-49","TXNDC17","PCGF1","SLC66A3","COQ2","SNAPC3","CCDC183",
#   "ASPM","HSPA5","LIN28B","GATA3","NPTX1","FAM110B","FBLN7","MLLT1","INPP5K","LMNTD2",
#   "TOMM5","DYRK2","C1orf198","EGFLAM","DCAF12","MKKS","IGHV3-33","XPO6","ZBTB16","MEI1",
#   "ENSG00000205045","SOD1","TRBV7-9","SERF1A","EEF1B2","TRIM40","RAC2","PPFIA1","SMG9",
#   "CES3","NUDC","BEND2","KIF1A","HBA1","APEX1","PEX16","SLFN5","RIPK3","CACNA1A","DOCK6",
#   "STIP1","GPR160","SDK1","TIMM50","ISG15","MRPS17","IQCM","ODF3B","BEND4","PRDM1","ST3GAL6",
#   "GABRD","EPHA7","NEURL2","TNFRSF8","TRAPPC5","RPS3A","ZNF354A","RNF146","PHTF2","SNRPD3",
#   "LBH","SH2D2A","ALPK3","SDC2","FCGR1A","RAB39A","TNFAIP1","RPL6","COQ8A","C19orf81",
#   "ZFP36","NIP7","H1-5","INSIG1","GLT8D1","CCDC115","TMEM154","HLA-DRB5","ANKRD36","KLF4",
#   "PPP1R14B","CYP7B1","CFH","PDIA4","BCAP29","GLB1L","HSP90AB1","IGLV8-61","PROS1","TRPT1",
#   "XPR1","MAF"
# )


# Extract patient IDs from the log2FC data column headers
# Assuming the first column is 'gene_id', adjust if different
logfc_columns <- colnames(data)[-1]  # Exclude the gene identifier column
patient_ids_logfc <- str_extract(logfc_columns, "(?<=patient_)\\d+(?=_log2FC)")

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros

# Filter the log2FC data for genes in the pathway
# Assuming the first column is 'gene_id'; change if different
pathway_data <- data %>%
  filter(Gene %in% gene_list_1000)

# Check if any genes from the gene_list are missing in the data
missing_genes <- setdiff(gene_list_1000, pathway_data$Gene)
if(length(missing_genes) > 0){
  warning(paste("The following genes are missing in the log2FC data:", paste(missing_genes, collapse = ", ")))
}

# Calculate the mean log2FC for the pathway genes for each patient
# Transpose the data to have patients as rows and genes as columns
pathway_activity <- pathway_data %>%
  select(-Gene) %>%  # Remove gene identifier column
  pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
  mutate(patient_id = str_extract(patient_col, "\\d+")) %>%  # Extract numeric patient ID
  group_by(patient_id) %>%
  summarize(pathway_activity = mean(log2FC, na.rm = TRUE))  # Calculate mean log2FC

# Merge the pathway_activity with cox_input_df
cox_input_df <- cox_input_df %>%
  left_join(pathway_activity, by = "patient_id")


# Convert pathway_activity to a binary variable based on median
median_activity <- median(cox_input_df$pathway_activity, na.rm = TRUE)

cox_input_df <- cox_input_df %>%
  mutate(pathway_activity_binary = ifelse(pathway_activity > median_activity, "Yes", "No")) %>%
  mutate(pathway_activity_binary = factor(pathway_activity_binary, levels = c("No", "Yes")))

# Optionally, remove the continuous pathway_activity column if no longer needed
# cox_input_df <- cox_input_df %>% select(-pathway_activity)

# Now, 'cox_input_df' contains the 'pathway_activity_binary' column

# Perform survival analysis using the binary pathway activity as a covariate
cox_ph_result <- get_cox(cox_input_df, "OS.months.", "Dead", c("Age", "pathway_activity_binary"))

# Display the summary of the Cox proportional hazards model
summary(cox_ph_result$cox)
# Now, 'cox_input_df' contains the 'pathway_activity' column

# # Perform survival analysis using Arm as a covariate
# cox_ph_result <- get_cox(cox_input_df, "OS.months.", "Dead", c("Age", "Sex", "pathway_activity"))
# 
# # Display the summary of the Cox proportional hazards model
# summary(cox_ph_result$cox)







##############################################################################################################
# trying to build a relation between innate and adaptive immunity
# Load necessary libraries
library(tidyverse)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(Rcpp)
library(cowplot)
library(stringr)  # For string manipulation
library(Seurat)

# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/scaled_normalized_log2FC_C1_vs_Pre.csv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

gene_list_500 <- c(
  "C1QC","C1QB","LYPD2","S100A12","VSIG4","VCAN","SLCO5A1","FOS","S100A8","C1QA","CD163",
  "ADAMTS2","HFM1","VMO1","CEBPD","GPX1","RBP7","AREG","THBS1","S100A9","GNLY","PLD4",
  "IL32","CYP4F22","GZMB","CD247","HES4","TPST1","MS4A6A","PPBP","G0S2","CD14","SCGB3A1",
  "TEX14","ITGA2B","FMN1","FLT3","MYL9","IGHV1-46","NKG7","ACSL1","FKBP5","H3C10","NRGN",
  "CLU","TFEC","TSPAN33","SASH1","CAVIN2","TREML1","RNASE2","LPL","SPON2","CKB","GZMH",
  "PRF1","DDIT4","HSH2D","CDKN1C","HBB","SH2D1B","DOCK4","EGR1","F13A1","IGLV2-14",
  "CTSW","RBFOX2","CST7","PLCB1","IGLV5-45","SLAMF7","LMNA","GZMA","SELPLG","DUSP5",
  "SYTL1","FGFBP2","MPIG6B","PRKN","PF4","SHTN1","ARL4A","PRKAR2B","TSTD1","CMTM5",
  "SPARC","CRIP1","LYZ","STAT4","ARID5B","SPN","IRS2","IL7R","NRG1","STAB1","APOO",
  "ADGRB3","EIF4G3","BLVRB","MTURN","KLRB1","ASGR2","CHST2","BAIAP2L1","PGRMC1","IGKV2-30",
  "MT-ATP8","ETS1","MSR1","OAS1","CSF3R","SLC44A2","C2orf88","BEX3","TRBC1","NFKBIA",
  "SETBP1","HBG1","SLC2A6","RNASE6","CUX1","IL13RA1","SAMSN1","IGLV3-27","SMOX","FRMD4B",
  "ETS2","CD36","CD99","IGLV1-51","TRIM58","PID1","SH3BGRL2","DUSP1","CD7","PIK3CG",
  "HLA-DMB","H2AC6","GBP2","LGALS2","CPED1","TRBC2","FGR","MCEMP1","SIGLEC10","IGKV1-27",
  "IL1R2","TMEM91","CTSC","SCART1","MGST1","GPR183","ITGAM","CCL4","ENTPD1","IL21R",
  "TRBV29-1","MLXIP","KLF12","RHOC","EPB41L3","IGHV4-59","LINGO2","CCR7","SMIM7","LAT2",
  "TNFAIP3","P2RY10","SSBP4","MMD","EPHX3","HCAR3","FAM110A","TMEM40","NAMPT",
  "ENSG00000285589","TRGV9","CD101","HBEGF","CTSL","C1orf167","PER1","LTB","ZNF467","SAP30",
  "RAB39B","CD9","CBLB","CDH23","AHR","PKN1","KLF10","SYAP1","JUN","PIGQ","VAMP8","GNG11",
  "S1PR4","ADAM8","PTPRCAP","RETN","CREB5","WARS1","GP9","SH2D3C","TUBB1","AKT3","CALHM2",
  "LY6E","TSHZ2","MXD1","RGL4","ITGB3","PTK2","ENSG00000259529","APOBEC3G","H2BC11","ACRBP",
  "SIRPA","TMEM128","PRKCQ","SPSB3","UPK3A","KNDC1","TSC22D1","RPLP0","NXF3","ADORA2A",
  "IFITM1","NCF1","STON2","ZAP70","UNC119","EEF1G","PDGFA","SH3BP4","GFOD1","IFIT2",
  "SNED1","CEBPA","S100P","FN1","TMEM176B","LPCAT3","IGHV2-70D","ITGAL","TTC39B",
  "IGLV4-60","RPL41","CLIC3","SNCA","LCK","CD300LF","CCR2","HOMER3","LY9","CSF1R","LGALSL",
  "LAPTM4B","POLD4","PTGDR2","GATD3B","ANPEP","SYNGR2","SLIT2","AP1S2","NT5M","TRBV6-6",
  "PPM1N","ARAP2","LFNG","LMBRD1","GBP5","TTC9","DRAP1","TP53I11","ADGRE2","ZDHHC1","RHOF",
  "MAL","COMT","TCIRG1","TRBV12-3","MYBL1","IGKV3D-11","ALOX12","GPRC5C","CD3E","ARHGAP24",
  "PPIF","TMEM120A","PTGS1","LDHD","CX3CR1","MEN1","DHPS","ZC2HC1A","TBC1D10C","TRDC",
  "ZNF532","STXBP2","ACAP1","TESC","IGLV2-5","YWHAG","PTPRE","CD244","S1PR5","SLC2A9",
  "UCP2","RNASE1","ODAD3","APOE","MCM6","H3C2","PEX5","PADI4","PLA2G4A","CAMTA1",
  "PLA2G12A","CLIC4","PLEKHA7","MARVELD1","SAV1","OBSCN","SIGLEC11","SPIDR","SLC38A7",
  "IGKV4-1","C1RL","ERAP1","ARL6IP5","TBCD","DNM3","CAV2","SLC22A16","IFI6",
  "ENSG00000124593","GBP4","JDP2","CFD","IGF2BP1","RPLP1","PECAM1","PRKG1","ANO5","AMOT",
  "CRISPLD2","IGHV5-51","H4C8","SNX18","MAP3K1","ID4","MERTK","SLC15A4","MS4A4E","HOXB2",
  "AARS1","LAIR2","TRBV25-1","GYPC","NAIP","H2BC4","TSPAN12","RASAL3","MICAL2","GPR174",
  "PTPN6","NUDCD1","TRBV2","COL26A1","HLA-DMA","ACTN1","BICRAL","MFSD2B","ELOVL7","TFPI",
  "CAMK4","H2AC11","NR4A3","NME4","PAQR7","PHF19","FCER1A","CSTA","IGLV1-47","LEF1","MCFD2",
  "UNC93B1","DGKH","DDT","RETREG1","CLPB","IGHV4-31","OST4","TRIM3","CD93","SLC40A1",
  "ARMCX5","ANO4","NPIPB4","C10orf143","HECTD2","SAPCD2","NCOA3","PFDN1","IER3IP1",
  "SLC2A3","CES1","LIME1","CXCL16","SHCBP1","TRBV18","MAFIP","LYN","YBX3","ZNF627","CMIP",
  "ZNF331","TRBV7-6","RPS8","RNF11","NBPF26","GALNT16","FCN1","TRPS1","CCDC9B","PLBD1",
  "MAD1L1","PF4V1","ENSG00000280571","CDIP1","TRAF3IP3","CPNE8","SFTPD","FLT3LG","DYSF",
  "TNFSF10","SYT17","MYO7A","SHISA5","CASP3","POLM","PRPF4","GPC6","ABHD16A","DUS2",
  "IGHV1-3","CDC45","SMARCD3","CAMKK1","PCK2","HSP90B1","TMC6","PLSCR3","PLCL1","MTSS1",
  "PDLIM1","RPS13","SERPINE1","H3C4","LIMD1","KCNB2","MNDA","ZFP36L1","PRR5L","OAS2","SKAP1",
  "CRTAP","MEIS1","CARD19","CERK","CCNDBP1","BIN1","TRMT1","PITPNB","RAB9B","PTPN4","ENAH",
  "RUNX3","LEPROTL1","NR6A1","HBA2","RAB13","PTMS","ADGRE5","RNF24","ARMC12","MKI67",
  "EPB41L1","SUN2","DBNDD2","ENSG00000272442","TRBV5-5","GSTM1","DPEP2","TNS3","SCML4",
  "STK25","GZMM"
)

logfc_columns <- colnames(data)[-1]  # Exclude the gene identifier column
patient_ids_logfc <- str_extract(logfc_columns, "(?<=patient_)\\d+(?=_log2FC)")

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros

# Filter the log2FC data for genes in the pathway
# Assuming the first column is 'gene_id'; change if different
pathway_data <- data %>%
  filter(Gene %in% gene_list_500)

# Check if any genes from the gene_list are missing in the data
missing_genes <- setdiff(gene_list_500, pathway_data$Gene)
if(length(missing_genes) > 0){
  warning(paste("The following genes are missing in the log2FC data:", paste(missing_genes, collapse = ", ")))
}

# Calculate the mean log2FC for the pathway genes for each patient
# Transpose the data to have patients as rows and genes as columns
pathway_activity <- pathway_data %>%
  select(-Gene) %>%  # Remove gene identifier column
  pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
  mutate(patient_id = str_extract(patient_col, "\\d+")) %>%  # Extract numeric patient ID
  group_by(patient_id) %>%
  summarize(pathway_activity = mean(log2FC, na.rm = TRUE))  # Calculate mean log2FC

# Merge the pathway_activity with cox_input_df
cox_input_df <- cox_input_df %>%
  left_join(pathway_activity, by = "patient_id")

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

create_survival_data <- function(gmt_file, pathway_name, seurat_obj, survival_data) {
  # Read GMT file and extract genes for the specified pathway
  gmt <- readLines(gmt_file)
  pathway_genes <- NULL
  for (line in gmt) {
    split_line <- strsplit(line, "\t")[[1]]
    if (split_line[1] == pathway_name) {
      pathway_genes <- split_line[-c(1,2)]  # Assuming the first two columns are pathway name and description
      break
    }
  }
  
  # Check if pathway was found
  if (is.null(pathway_genes)) {
    stop("Pathway not found in GMT file.")
  }
  
  # Read survival data and sort patients by survival
  survival_df <- read.csv(survival_data)
  sorted_patients <- survival_df[order(-survival_df$OS.months.), 'patient_id']
  print(sorted_patients)
  
  # Create a grid of feature plots
  timepoints <- c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36")
  # Initialize columns for mean expression at each timepoint in survival_df
  for(tp in timepoints) {
    survival_df[[paste0("Mean_Expr_", tp)]] <- NA_real_
  }
  
  for (patient in sorted_patients) {
    # print(patient)
    patient_specific_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$donor == patient,])
    if (length(patient_specific_cells) == 0) {
      next
    }
    patient_seurat_obj <- subset(seurat_obj, cells = patient_specific_cells)
    # Get the genes present in the Seurat object
    seurat_genes <- rownames(patient_seurat_obj@assays$RNA@counts)
    
    # Find intersection of pathway genes and Seurat object genes
    common_genes <- intersect(pathway_genes, seurat_genes)
    
    # Check if there are any common genes
    if (length(common_genes) == 0) {
      stop("None of the pathway genes are found in the Seurat object.")
    }
    
    # Calculate average expression of pathway genes per cell
    # Note: This assumes that the data is already normalized
    pathway_avg_expression <- tryCatch({
      apply(GetAssayData(patient_seurat_obj, assay = "RNA", slot = "data")[common_genes, ], 2, mean, na.rm = TRUE)
    }, error = function(e) {
      for (tp in timepoints) {
        # Update survival_df with mean expression for this patient and timepoint
        survival_df[survival_df$patient_id == patient, paste0("Mean_Expr_", tp)] <- NA
      }
      return(NULL)
    })
    
    # Skip the rest of the loop if an error occurred
    if (is.null(pathway_avg_expression)) {
      next
    }
    
    # # Check for cells at timepoint N
    # baseline_timepoint <- intersect(timepoints, unique(patient_seurat_obj@meta.data$TimePoint))[1]
    # timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == baseline_timepoint, ])
    # 
    # if (length(timepoint_N_barcodes) > 0) {
    #   # Calculate mean and standard deviation for cells at timepoint N
    #   mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
    #   sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])
    #   print(sd_N)
    #   if (sd_N == 0 | is.na(sd_N)){
    #     scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    #   } else {
    #     # Scale based on timepoint N
    #     scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
    #   }
    # } else {
    #   # Use scale function if no cells at timepoint N
    #   scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    # }
    # 
    # pathway_avg_expression <- scaled_pathway_avg_expression
    names(pathway_avg_expression) <- colnames(patient_seurat_obj)
    
    # Add this as a metadata column
    patient_seurat_obj[["pathway_avg_expression"]] <- pathway_avg_expression
    
    for (tp in timepoints) {
      patient_tp_cells <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == tp, ])
      if (length(patient_tp_cells) > 0) {
        # Calculate mean expression for this timepoint
        mean_expression <- mean(pathway_avg_expression[patient_tp_cells], na.rm = TRUE)
      } else {
        mean_expression <- NA
      }
      # Update survival_df with mean expression for this patient and timepoint
      survival_df[survival_df$patient_id == patient, paste0("Mean_Expr_", tp)] <- mean_expression
    }
  }
  # Return the updated survival dataframe
  return(survival_df)
}


seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF" & survival_df$IDH != "POS",]

patient_col <- "Patient"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

include_cluster_proportion <- FALSE
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")
pathway <- "GO_ADAPTIVE_IMMUNE_RESPONSE"
method = "ratio"

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
celltype <- "Proliferating_Effector"
cluster_list <- as.vector(celltype_to_cluster[celltype_to_cluster$celltype == celltype, "cluster"])
seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cluster_list)
seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
# seurat_object_t_cell_subset <- seurat_object_t_cells
tcell_cluster_list <- cluster_list
seurat_metadata <- seurat_metadata_t_cells

updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
updated_survival_df <- updated_survival_df[updated_survival_df$IDH != "POS",]

temp_survival_df <- updated_survival_df
timepoint_a <- "C1"
timepoint_b <- "C2"

if (include_cluster_proportion) {
  timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
  timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
  temp_survival_df <- merge(temp_survival_df, timepoint_a_cluster_proportion_df, by = "patient_id")
  temp_survival_df <- merge(temp_survival_df, timepoint_b_cluster_proportion_df, by = "patient_id")
  temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_a)]
  temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)] * temp_survival_df[, paste0("Cluster_Proportion_", timepoint_b)]
  if (method == "difference") {
    temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
  } else if (method == "ratio") {
    temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
  }
} else {
  temp_survival_df[paste0("Pathway_Signal_", timepoint_a)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_a)]
  temp_survival_df[paste0("Pathway_Signal_", timepoint_b)] <- temp_survival_df[, paste0("Mean_Expr_", timepoint_b)]
  if (method == "difference") {
    temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) - (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
  } else if (method == "ratio") {
    temp_survival_df[, "signal_change"] <- (temp_survival_df[paste0("Pathway_Signal_", timepoint_b)]) / (temp_survival_df[paste0("Pathway_Signal_", timepoint_a)])
  }
}

temp_survival_df <- temp_survival_df[temp_survival_df$site == "UF" & temp_survival_df$Arm == "MK-3475 + MLA", ]

merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity", "OS.months.")], temp_survival_df[, c("patient_id", "signal_change")], by = "patient_id")
# Calculate Pearson correlation
correlation <- cor(merged_df$pathway_activity, merged_df$signal_change, method = "pearson")

# Print the correlation
print(correlation)


# couldn't find correlation between logFC and T cell activation or Adaptive immune response
######################################################################################################################################################
# now checking to see if there is correlation between logFC and CD8 clonal expansion







###########################################################################################################
# Load necessary libraries
library(tidyverse)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(Rcpp)
library(cowplot)
library(stringr)  # For string manipulation
library(Seurat)

# 1. Read the two Seurat objects from RDS files
seurat_obj1 <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_control_Pre_C1.RDS")
seurat_obj2 <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")

# 2. Merge the two Seurat objects into one
#    The 'add.cell.ids' parameter is useful for distinguishing cells from each original object.
#    The 'project' parameter can be set to a meaningful name for your combined dataset.
combined_seurat <- merge(
  x = seurat_obj1, 
  y = seurat_obj2, 
  add.cell.ids = c("Sample1", "Sample2"), 
  project = "CombinedProject"
)

# 3. (Optional) Check the combined metadata
# head(combined_seurat@meta.data)

# 4. Proceed with standard workflows on the combined object (QC, normalization, etc.)

# Check the metadata
metadata <- combined_seurat@meta.data
head(metadata)

# Ensure that 'TimePoint' and 'Patient' columns exist
if(!all(c("TimePoint", "Patient") %in% colnames(metadata))) {
  stop("The Seurat object does not contain 'TimePoint' and/or 'Patient' columns in metadata.")
}

# 2. Prepare pseudo-bulk RNA-seq counts
# Aggregate counts by Patient and TimePoint

# Extract the raw counts matrix
# Commonly in the "RNA" assay
counts <- GetAssayData(combined_seurat, assay = "RNA", slot = "counts")

# Add Patient and TimePoint information to the metadata
metadata <- combined_seurat@meta.data %>%
  select(Patient, TimePoint)

# **Modified Step:** Create a unique identifier for each pseudo-bulk sample with "patient_" prefix
metadata$Sample <- paste("patient", metadata$Patient, metadata$TimePoint, sep = "_")
# Example: "patient_1_Pre", "patient_1_C1", "patient_2_Pre", etc.

# Sum the counts for cells belonging to the same Sample
# Transpose counts to have cells as rows and genes as columns for aggregation
counts_df <- as.data.frame(t(as.matrix(counts)))
counts_df$Sample <- metadata$Sample

# Aggregate counts by Sample
pseudo_bulk_counts <- counts_df %>%
  group_by(Sample) %>%
  summarise_all(sum)

# **Corrected Step:** Convert back to a matrix with genes as rows and samples as columns
# Since 'pseudo_bulk_counts' has samples as rows and genes as columns,
# we need to transpose it to have genes as rows and samples as columns.

# Remove the 'Sample' column and transpose
pseudo_bulk_counts_mat <- t(as.matrix(pseudo_bulk_counts[,-1]))

# Assign gene names as rownames
rownames(pseudo_bulk_counts_mat) <- rownames(counts)  # Correct: Gene names

# Assign sample names as column names
colnames(pseudo_bulk_counts_mat) <- pseudo_bulk_counts$Sample

# Verify dimensions
dim(pseudo_bulk_counts_mat)  # Should be genes x samples (e.g., 19343 x 24)

# write.csv(pseudo_bulk_counts_mat, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_pre_C1_pseudo_bulk_RNASeq.csv")

# Specify the path to your pseudo_bulk_counts_mat.csv
counts_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_pre_C1_pseudo_bulk_RNASeq.csv"

# Read the CSV file
# Assuming the first column contains gene names and is set as rownames
counts_df <- read.csv(counts_csv_path, row.names = 1, check.names = FALSE)

# Verify the dimensions and a snippet of the data
cat("Dimensions of the counts matrix:", dim(counts_df), "\n")
head(counts_df)


# 4. Normalize Counts Using "Pre" as Baseline
# -------------------------------------------

# Subset the counts matrix to include only selected genes
selected_counts <- counts_df

# Parse sample names to extract Patient and TimePoint information
# Assuming sample names are in the format "patient_X_TimePoint", e.g., "patient_1_Pre"
sample_names <- colnames(selected_counts)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = sample_names,
  stringsAsFactors = FALSE
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_", remove = FALSE) %>%
  mutate(Patient = paste(Patient_Prefix, Patient_ID, sep = "_")) %>%
  select(Sample, Patient, TimePoint)

# View sample information
print(sample_info)

# Get a list of unique patients
patients <- unique(sample_info$Patient)
cat("Number of unique patients:", length(patients), "\n")

# Initialize a list to store log2 fold changes
log2fc_list <- list()

# Loop through each patient to compute log2FC
for (patient in patients) {
  # Identify "Pre" and "C1" samples for the patient
  pre_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "Pre") %>%
    pull(Sample)
  
  c1_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "C1") %>%
    pull(Sample)
  
  # Check if both samples are present
  if (length(pre_sample) == 1 & length(c1_sample) == 1) {
    # Extract counts for "Pre" and "C1" samples
    pre_counts <- selected_counts[, pre_sample]
    c1_counts <- selected_counts[, c1_sample]
    
    # Compute log2 fold change with a pseudocount of 1 to avoid division by zero
    log2fc <- log2((c1_counts + 1) / (pre_counts + 1))
    
    # Ensure that log2fc is a named vector with gene names
    names(log2fc) <- rownames(selected_counts)
    
    # Store in the list
    log2fc_list[[patient]] <- log2fc
  } else {
    warning(paste("Missing 'Pre' or 'C1' sample for patient:", patient))
  }
}

# Combine the list into a normalized matrix
# Rows: Genes, Columns: Patients (log2FC)
normalized_mat <- do.call(cbind, log2fc_list)

# Rename columns to indicate log2FC
colnames(normalized_mat) <- paste0(colnames(normalized_mat), "_log2FC")

# Inspect the normalized matrix
dim(normalized_mat)
head(normalized_mat)

# 5. Scale the Normalized Data
# ----------------------------

# Scaling can be done per gene to standardize the log2FC across patients
# This is useful for visualization purposes

# Transpose, scale, and transpose back
scaled_mat <- t(scale(t(normalized_mat)))

# Replace any NA values resulting from scaling with 0
scaled_mat[is.na(scaled_mat)] <- 0

# Inspect the scaled matrix
dim(scaled_mat)
head(scaled_mat)

# 6. Save the Normalized and Scaled Data as CSV
# ---------------------------------------------

# Prepare the normalized data for saving
normalized_df <- as.data.frame(normalized_mat)
normalized_df$Gene <- rownames(normalized_df)

# Reorder columns to have 'Gene' as the first column
normalized_df <- normalized_df %>%
  select(Gene, everything())

# Save the normalized log2FC data
write.csv(normalized_df, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Normalized log2 fold change data saved to 'normalized_log2FC_C1_vs_Pre.csv'.\n")

# Prepare the scaled data for saving
scaled_df <- as.data.frame(scaled_mat)
scaled_df$Gene <- rownames(scaled_df)

# Reorder columns to have 'Gene' as the first column
scaled_df <- scaled_df %>%
  select(Gene, everything())

# Save the scaled normalized data
write.csv(scaled_df, file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_scaled_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Scaled normalized data saved to 'scaled_normalized_log2FC_C1_vs_Pre.csv'.\n")







########################################################################################################################################################################
# trying to see if there is any correlation between innate and adaptive immunity
########################################################################################################################################################################

# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_scaled_normalized_log2FC_C1_vs_Pre.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# # Read the log2FC data
# data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_normalized_log2FC_C1_vs_Pre.csv",
#                  header = TRUE,
#                  stringsAsFactors = FALSE)

gene_list_500 <- c(
  "C1QC","C1QB","LYPD2","S100A12","VSIG4","VCAN","SLCO5A1","FOS","S100A8","C1QA","CD163",
  "ADAMTS2","HFM1","VMO1","CEBPD","GPX1","RBP7","AREG","THBS1","S100A9","GNLY","PLD4",
  "IL32","CYP4F22","GZMB","CD247","HES4","TPST1","MS4A6A","PPBP","G0S2","CD14","SCGB3A1",
  "TEX14","ITGA2B","FMN1","FLT3","MYL9","IGHV1-46","NKG7","ACSL1","FKBP5","H3C10","NRGN",
  "CLU","TFEC","TSPAN33","SASH1","CAVIN2","TREML1","RNASE2","LPL","SPON2","CKB","GZMH",
  "PRF1","DDIT4","HSH2D","CDKN1C","HBB","SH2D1B","DOCK4","EGR1","F13A1","IGLV2-14",
  "CTSW","RBFOX2","CST7","PLCB1","IGLV5-45","SLAMF7","LMNA","GZMA","SELPLG","DUSP5",
  "SYTL1","FGFBP2","MPIG6B","PRKN","PF4","SHTN1","ARL4A","PRKAR2B","TSTD1","CMTM5",
  "SPARC","CRIP1","LYZ","STAT4","ARID5B","SPN","IRS2","IL7R","NRG1","STAB1","APOO",
  "ADGRB3","EIF4G3","BLVRB","MTURN","KLRB1","ASGR2","CHST2","BAIAP2L1","PGRMC1","IGKV2-30",
  "MT-ATP8","ETS1","MSR1","OAS1","CSF3R","SLC44A2","C2orf88","BEX3","TRBC1","NFKBIA",
  "SETBP1","HBG1","SLC2A6","RNASE6","CUX1","IL13RA1","SAMSN1","IGLV3-27","SMOX","FRMD4B",
  "ETS2","CD36","CD99","IGLV1-51","TRIM58","PID1","SH3BGRL2","DUSP1","CD7","PIK3CG",
  "HLA-DMB","H2AC6","GBP2","LGALS2","CPED1","TRBC2","FGR","MCEMP1","SIGLEC10","IGKV1-27",
  "IL1R2","TMEM91","CTSC","SCART1","MGST1","GPR183","ITGAM","CCL4","ENTPD1","IL21R",
  "TRBV29-1","MLXIP","KLF12","RHOC","EPB41L3","IGHV4-59","LINGO2","CCR7","SMIM7","LAT2",
  "TNFAIP3","P2RY10","SSBP4","MMD","EPHX3","HCAR3","FAM110A","TMEM40","NAMPT",
  "ENSG00000285589","TRGV9","CD101","HBEGF","CTSL","C1orf167","PER1","LTB","ZNF467","SAP30",
  "RAB39B","CD9","CBLB","CDH23","AHR","PKN1","KLF10","SYAP1","JUN","PIGQ","VAMP8","GNG11",
  "S1PR4","ADAM8","PTPRCAP","RETN","CREB5","WARS1","GP9","SH2D3C","TUBB1","AKT3","CALHM2",
  "LY6E","TSHZ2","MXD1","RGL4","ITGB3","PTK2","ENSG00000259529","APOBEC3G","H2BC11","ACRBP",
  "SIRPA","TMEM128","PRKCQ","SPSB3","UPK3A","KNDC1","TSC22D1","RPLP0","NXF3","ADORA2A",
  "IFITM1","NCF1","STON2","ZAP70","UNC119","EEF1G","PDGFA","SH3BP4","GFOD1","IFIT2",
  "SNED1","CEBPA","S100P","FN1","TMEM176B","LPCAT3","IGHV2-70D","ITGAL","TTC39B",
  "IGLV4-60","RPL41","CLIC3","SNCA","LCK","CD300LF","CCR2","HOMER3","LY9","CSF1R","LGALSL",
  "LAPTM4B","POLD4","PTGDR2","GATD3B","ANPEP","SYNGR2","SLIT2","AP1S2","NT5M","TRBV6-6",
  "PPM1N","ARAP2","LFNG","LMBRD1","GBP5","TTC9","DRAP1","TP53I11","ADGRE2","ZDHHC1","RHOF",
  "MAL","COMT","TCIRG1","TRBV12-3","MYBL1","IGKV3D-11","ALOX12","GPRC5C","CD3E","ARHGAP24",
  "PPIF","TMEM120A","PTGS1","LDHD","CX3CR1","MEN1","DHPS","ZC2HC1A","TBC1D10C","TRDC",
  "ZNF532","STXBP2","ACAP1","TESC","IGLV2-5","YWHAG","PTPRE","CD244","S1PR5","SLC2A9",
  "UCP2","RNASE1","ODAD3","APOE","MCM6","H3C2","PEX5","PADI4","PLA2G4A","CAMTA1",
  "PLA2G12A","CLIC4","PLEKHA7","MARVELD1","SAV1","OBSCN","SIGLEC11","SPIDR","SLC38A7",
  "IGKV4-1","C1RL","ERAP1","ARL6IP5","TBCD","DNM3","CAV2","SLC22A16","IFI6",
  "ENSG00000124593","GBP4","JDP2","CFD","IGF2BP1","RPLP1","PECAM1","PRKG1","ANO5","AMOT",
  "CRISPLD2","IGHV5-51","H4C8","SNX18","MAP3K1","ID4","MERTK","SLC15A4","MS4A4E","HOXB2",
  "AARS1","LAIR2","TRBV25-1","GYPC","NAIP","H2BC4","TSPAN12","RASAL3","MICAL2","GPR174",
  "PTPN6","NUDCD1","TRBV2","COL26A1","HLA-DMA","ACTN1","BICRAL","MFSD2B","ELOVL7","TFPI",
  "CAMK4","H2AC11","NR4A3","NME4","PAQR7","PHF19","FCER1A","CSTA","IGLV1-47","LEF1","MCFD2",
  "UNC93B1","DGKH","DDT","RETREG1","CLPB","IGHV4-31","OST4","TRIM3","CD93","SLC40A1",
  "ARMCX5","ANO4","NPIPB4","C10orf143","HECTD2","SAPCD2","NCOA3","PFDN1","IER3IP1",
  "SLC2A3","CES1","LIME1","CXCL16","SHCBP1","TRBV18","MAFIP","LYN","YBX3","ZNF627","CMIP",
  "ZNF331","TRBV7-6","RPS8","RNF11","NBPF26","GALNT16","FCN1","TRPS1","CCDC9B","PLBD1",
  "MAD1L1","PF4V1","ENSG00000280571","CDIP1","TRAF3IP3","CPNE8","SFTPD","FLT3LG","DYSF",
  "TNFSF10","SYT17","MYO7A","SHISA5","CASP3","POLM","PRPF4","GPC6","ABHD16A","DUS2",
  "IGHV1-3","CDC45","SMARCD3","CAMKK1","PCK2","HSP90B1","TMC6","PLSCR3","PLCL1","MTSS1",
  "PDLIM1","RPS13","SERPINE1","H3C4","LIMD1","KCNB2","MNDA","ZFP36L1","PRR5L","OAS2","SKAP1",
  "CRTAP","MEIS1","CARD19","CERK","CCNDBP1","BIN1","TRMT1","PITPNB","RAB9B","PTPN4","ENAH",
  "RUNX3","LEPROTL1","NR6A1","HBA2","RAB13","PTMS","ADGRE5","RNF24","ARMC12","MKI67",
  "EPB41L1","SUN2","DBNDD2","ENSG00000272442","TRBV5-5","GSTM1","DPEP2","TNS3","SCML4",
  "STK25","GZMM"
)

logfc_columns <- colnames(data)[-1]  # Exclude the gene identifier column
patient_ids_logfc <- str_extract(logfc_columns, "(?<=patient_)\\d+(?=_log2FC)")

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF", ]
# cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros

# Filter the log2FC data for genes in the pathway
# Assuming the first column is 'gene_id'; change if different
pathway_data <- data %>%
  filter(Gene %in% gene_list_500)

# Check if any genes from the gene_list are missing in the data
missing_genes <- setdiff(gene_list_500, pathway_data$Gene)
if(length(missing_genes) > 0){
  warning(paste("The following genes are missing in the log2FC data:", paste(missing_genes, collapse = ", ")))
}

# Calculate the mean log2FC for the pathway genes for each patient
# Transpose the data to have patients as rows and genes as columns
pathway_activity <- pathway_data %>%
  select(-Gene) %>%  # Remove gene identifier column
  pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
  mutate(patient_id = str_extract(patient_col, "\\d+")) %>%  # Extract numeric patient ID
  group_by(patient_id) %>%
  summarize(pathway_activity = mean(log2FC, na.rm = TRUE))  # Calculate mean log2FC

# Merge the pathway_activity with cox_input_df
cox_input_df <- cox_input_df %>%
  left_join(pathway_activity, by = "patient_id")
# cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9), ]
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9, 1, 4), ]



clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_All_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion$patient_id <- rownames(clonal_expansion)
merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity", "Arm", "MGMT", "IDH", "Sex", "Age")], clonal_expansion, by = "patient_id")

merged_df <- merged_df[!is.na(merged_df$pathway_activity) & !is.na(merged_df$clonal_expansion), ]

# spearman_value <- cor(merged_df$pathway_activity,
#                       merged_df$clonal_expansion,
#                       method = "spearman")
# 
# spearman_value

merged_df[merged_df$IDH == "UNK", "IDH"] = "NEG"
merged_df

#############################################################################################################
# gemini correlation metrics
#############################################################################################################
# Load the mccr package
library(mccr)
# Initialize list to store results
results <- list()

# 2. Both columns numeric
results$'Pearson (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "pearson")
results$'Spearman (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "spearman")
results$'Kendall (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "kendall")

# 3. Numeric pathway_activity vs. Binary clonal_expansion (> 1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
# Check for variance in the binary column
if (length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- cor(merged_df$pathway_activity, merged_df$clonal_binary_gt1, method = "pearson")
} else {
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- NA
  warning("Binary clonal expansion column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 4. Binary pathway_activity (> median) vs. Numeric clonal_expansion
median_pathway <- median(merged_df$pathway_activity)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity > median_pathway)
# Check for variance in the binary column
if (length(unique(merged_df$pathway_binary_median)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson")
} else {
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- NA
  warning("Binary pathway activity column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 5. Both columns binary
# Check for variance in both binary columns
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Phi coefficient (equivalent to Pearson on 0/1 binary data)
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  
  # Matthews Correlation Coefficient (MCC) using the mccr package
  # Note: mccr expects factors or numeric {0, 1} or {-1, 1}. Our 0/1 format works.
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
  
} else {
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  warning("One or both binary columns lack variance. Cannot calculate Phi or MCC.")
}


# Find the metric with the most negative correlation
best_metric_name <- NULL
highest_inverse_corr <- 0 # Initialize to 0, looking for most negative

for (metric_name in names(results)) {
  value <- results[[metric_name]]
  if (!is.na(value) && value < highest_inverse_corr) {
    highest_inverse_corr <- value
    best_metric_name <- metric_name
  }
}

# Print the results
cat("Correlation Results:\n\n")
for (metric_name in names(results)) {
  cat(sprintf("%s: %.4f\n", metric_name, results[[metric_name]]))
}

cat("\n--------------------------------------------------\n")
cat("Best Inverse Correlation Found:\n")
cat(sprintf("Method: %s\n", best_metric_name))
cat(sprintf("Correlation Value: %.4f\n", highest_inverse_corr))
cat("--------------------------------------------------\n")

# Display the median used for pathway activity binarization
cat(sprintf("\nMedian used for 'pathway_activity' binarization: %.4f\n", median_pathway))





###############################################################################################################
# applying linear model
###############################################################################################################

# 1) Make sure Sex is treated as a factor (categorical) -----------------------
merged_df$Sex <- factor(merged_df$Sex)      # ‘F’ will be the reference level by default

# 2) Fit the model ------------------------------------------------------------
fit <- lm(clonal_expansion ~ pathway_activity + Sex + Age, data = merged_df)

# 3) Inspect the results ------------------------------------------------------
summary(fit)            # classic R output
#– or –#
library(broom)
tidy(fit)               # neat tibble with estimates, SEs and p-values
###############################################

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# 2. Add a grouping column ----------------------------------------------------------
merged_df$Group <- dplyr::case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                          ~ "Uncategorised"
)

# 3. Plot with ggplot2 --------------------------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = pathway_activity,
               y = clonal_expansion,
               colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Pathway activity vs. clonal expansion",
       subtitle = "Points coloured by survivor group",
       x = "Pathway activity",
       y = "Clonal expansion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor"  = "#4DBBD5",
    "Uncategorised"       = "grey60"
  ))




################################################################################################################################
# starting with the DC Analysis (all DC Cells)
################################################################################################################################

############################################################################################################################################################################################
# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)

setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis")

# 1. Load the Seurat object
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_DC_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")

# Check the metadata
metadata <- seurat_obj@meta.data
head(metadata)

# Ensure that 'TimePoint' and 'Patient' columns exist
if(!all(c("TimePoint", "Patient") %in% colnames(metadata))) {
  stop("The Seurat object does not contain 'TimePoint' and/or 'Patient' columns in metadata.")
}

# 2. Prepare pseudo-bulk RNA-seq counts
# Aggregate counts by Patient and TimePoint

# Extract the raw counts matrix
# Commonly in the "RNA" assay
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Add Patient and TimePoint information to the metadata
metadata <- seurat_obj@meta.data %>%
  select(Patient, TimePoint)

# **Modified Step:** Create a unique identifier for each pseudo-bulk sample with "patient_" prefix
metadata$Sample <- paste("patient", metadata$Patient, metadata$TimePoint, sep = "_")
# Example: "patient_1_Pre", "patient_1_C1", "patient_2_Pre", etc.

# Sum the counts for cells belonging to the same Sample
# Transpose counts to have cells as rows and genes as columns for aggregation
counts_df <- as.data.frame(t(as.matrix(counts)))
counts_df$Sample <- metadata$Sample

# Aggregate counts by Sample
pseudo_bulk_counts <- counts_df %>%
  group_by(Sample) %>%
  summarise_all(sum)

# **Corrected Step:** Convert back to a matrix with genes as rows and samples as columns
# Since 'pseudo_bulk_counts' has samples as rows and genes as columns,
# we need to transpose it to have genes as rows and samples as columns.

# Remove the 'Sample' column and transpose
pseudo_bulk_counts_mat <- t(as.matrix(pseudo_bulk_counts[,-1]))

# Assign gene names as rownames
rownames(pseudo_bulk_counts_mat) <- rownames(counts)  # Correct: Gene names

# Assign sample names as column names
colnames(pseudo_bulk_counts_mat) <- pseudo_bulk_counts$Sample

# Verify dimensions
dim(pseudo_bulk_counts_mat)  # Should be genes x samples (e.g., 19343 x 24)

write.csv(pseudo_bulk_counts_mat, file = "exp_DC_pseudo_bulk_RNASeq.csv")

# 3. Prepare the design matrix for edgeR
# Extract Patient and TimePoint information for each Sample
sample_info <- data.frame(
  Sample = colnames(pseudo_bulk_counts_mat)
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_")

# **Modified Step:** Reconstruct the Patient identifier with the "patient_" prefix
sample_info$Patient <- paste(sample_info$Patient_Prefix, sample_info$Patient_ID, sep = "_")
# Example: "patient_1", "patient_2", etc.

# Ensure factors are correctly set
sample_info$Patient <- factor(sample_info$Patient)
sample_info$TimePoint <- factor(sample_info$TimePoint, levels = c("Pre", "C1"))

# 4. Create DGEList object
dge <- DGEList(counts = pseudo_bulk_counts_mat)

# 5. Normalize the data
dge <- calcNormFactors(dge)

# 6. Create the design matrix
# Model: ~ Patient + TimePoint
design <- model.matrix(~ Patient + TimePoint, data = sample_info)

# 7. Estimate dispersion
dge <- estimateDisp(dge, design)

# 8. Fit the model and perform differential expression
fit <- glmFit(dge, design)
# Assuming "TimePointC1" is the coefficient for C1 vs Pre
lrt <- glmLRT(fit, coef = "TimePointC1")

# 9. Extract differential expression results
de_results <- topTags(lrt, n = Inf)$table
de_results <- as.data.frame(de_results)
# Add gene names as a column
de_results$Gene <- rownames(de_results)

# 10. Save the differential expression results
write.csv(de_results, file = "exp_DC_differential_expression_results.csv", row.names = FALSE)

# 11. Create input file for GSEA
# Typically, GSEA requires a ranked list of genes. One common ranking metric is logFC multiplied by the sign of the p-value.

# Here, we'll use the signed -log10(p-value) multiplied by the sign of logFC
# Alternatively, you can choose other ranking metrics based on your preference

# To handle p-values of 0, add a small pseudocount
de_results$PValue[de_results$PValue == 0] <- 1e-300

de_results$RankingMetric <- -log10(de_results$PValue) * sign(de_results$logFC)

# Order the genes by the ranking metric
ranked_genes <- de_results %>%
  arrange(desc(RankingMetric)) %>%
  select(Gene, RankingMetric)

# Save the ranked list for GSEA
# **Option 1:** GCT format
# GCT format requires specific headers. Here's how to structure it:

# Create a GCT header
gct_header <- data.frame(
  NAME = c("#1.2"),
  DESCRIPTION = c("RankingMetric"),
  stringsAsFactors = FALSE
)

# Combine header and data
gct_data <- data.frame(
  NAME = ranked_genes$Gene,
  DESCRIPTION = ranked_genes$Gene,  # GSEA often expects gene symbols or descriptions here
  RankingMetric = ranked_genes$RankingMetric,
  stringsAsFactors = FALSE
)

# Write the GCT file
write.table(gct_header, file = "exp_DC_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gct_data, file = "exp_DC_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE)

# **Option 2:** RNK format
# RNK format is simpler and often preferred for ranked lists

write.table(ranked_genes, file = "exp_DC_GSEA_ranked_genes.rnk", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Inform the user that the process is complete
cat("Pseudo-bulk RNA-seq aggregation, differential expression analysis, and GSEA input file creation are complete.\n")
cat("Differential expression results saved to 'differential_expression_results.csv'.\n")
cat("GSEA ranked gene list saved to 'GSEA_ranked_genes.gct' and 'GSEA_ranked_genes.rnk'.\n")





######################################################################################
library(dplyr)
library(tidyr)
# 1. Read the Pseudo-Bulk Counts Matrix
# -------------------------------------

# Specify the path to your pseudo_bulk_counts_mat.csv
counts_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_pseudo_bulk_RNASeq.csv"

# Read the CSV file
# Assuming the first column contains gene names and is set as rownames
counts_df <- read.csv(counts_csv_path, row.names = 1, check.names = FALSE)

# Verify the dimensions and a snippet of the data
cat("Dimensions of the counts matrix:", dim(counts_df), "\n")
head(counts_df)


# 4. Normalize Counts Using "Pre" as Baseline
# -------------------------------------------

# Subset the counts matrix to include only selected genes
selected_counts <- counts_df

# Parse sample names to extract Patient and TimePoint information
# Assuming sample names are in the format "patient_X_TimePoint", e.g., "patient_1_Pre"
sample_names <- colnames(selected_counts)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = sample_names,
  stringsAsFactors = FALSE
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_", remove = FALSE) %>%
  mutate(Patient = paste(Patient_Prefix, Patient_ID, sep = "_")) %>%
  select(Sample, Patient, TimePoint)

# View sample information
print(sample_info)

# Get a list of unique patients
patients <- unique(sample_info$Patient)
cat("Number of unique patients:", length(patients), "\n")

# Initialize a list to store log2 fold changes
log2fc_list <- list()

# Loop through each patient to compute log2FC
for (patient in patients) {
  # Identify "Pre" and "C1" samples for the patient
  pre_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "Pre") %>%
    pull(Sample)
  
  c1_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "C1") %>%
    pull(Sample)
  
  # Check if both samples are present
  if (length(pre_sample) == 1 & length(c1_sample) == 1) {
    # Extract counts for "Pre" and "C1" samples
    pre_counts <- selected_counts[, pre_sample]
    c1_counts <- selected_counts[, c1_sample]
    
    # Compute log2 fold change with a pseudocount of 1 to avoid division by zero
    log2fc <- log2((c1_counts + 1) / (pre_counts + 1))
    
    # Ensure that log2fc is a named vector with gene names
    names(log2fc) <- rownames(selected_counts)
    
    # Store in the list
    log2fc_list[[patient]] <- log2fc
  } else {
    warning(paste("Missing 'Pre' or 'C1' sample for patient:", patient))
  }
}

# Combine the list into a normalized matrix
# Rows: Genes, Columns: Patients (log2FC)
normalized_mat <- do.call(cbind, log2fc_list)

# Rename columns to indicate log2FC
colnames(normalized_mat) <- paste0(colnames(normalized_mat), "_log2FC")

# Inspect the normalized matrix
dim(normalized_mat)
head(normalized_mat)

# 5. Scale the Normalized Data
# ----------------------------

# Scaling can be done per gene to standardize the log2FC across patients
# This is useful for visualization purposes

# Transpose, scale, and transpose back
scaled_mat <- t(scale(t(normalized_mat)))

# Replace any NA values resulting from scaling with 0
scaled_mat[is.na(scaled_mat)] <- 0

# Inspect the scaled matrix
dim(scaled_mat)
head(scaled_mat)

# 6. Save the Normalized and Scaled Data as CSV
# ---------------------------------------------

# Prepare the normalized data for saving
normalized_df <- as.data.frame(normalized_mat)
normalized_df$Gene <- rownames(normalized_df)

# Reorder columns to have 'Gene' as the first column
normalized_df <- normalized_df %>%
  select(Gene, everything())

# Save the normalized log2FC data
write.csv(normalized_df, file = "exp_DC_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Normalized log2 fold change data saved to 'exp_DC_normalized_log2FC_C1_vs_Pre.csv'.\n")

# Prepare the scaled data for saving
scaled_df <- as.data.frame(scaled_mat)
scaled_df$Gene <- rownames(scaled_df)

# Reorder columns to have 'Gene' as the first column
scaled_df <- scaled_df %>%
  select(Gene, everything())

# Save the scaled normalized data
write.csv(scaled_df, file = "exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Scaled normalized data saved to 'exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv'.\n")



#########################################################################################################################################################
# Plotting heatmap

# Install necessary packages if not installed
# install.packages("pheatmap")

library(pheatmap)

# 1. Read the file
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv", header = TRUE, stringsAsFactors = FALSE)

# The columns are expected to be:
# Gene, patient_10_log2FC, patient_12_log2FC, patient_13_log2FC, patient_14_log2FC,
# patient_18_log2FC, patient_19_log2FC, patient_20_log2FC, patient_21_log2FC,
# patient_2_log2FC, patient_3_log2FC, patient_5_log2FC, patient_7_log2FC

# 2. Define the gene list
gene_list <- c(
  "EREG", "IGLV1-40", "AREG", "TAL1", "CAVIN2", "S100A12", "MYL9", "SLCO5A1", "IGKV3-11", "COL24A1",
  "LAMP5", "COL26A1", "SZT2", "S100A8", "CUX2", "TTC24", "CUEDC1", "NAALADL2", "ITGA2B", "IGLV2-11",
  "CLEC1B", "CRYM", "IRS2", "TRAF4", "NRGN", "PPBP", "THBS1", "CRYAA2", "COBLL1", "AKAP6",
  "SETBP1", "TAMALIN", "ADGRB3", "PMEPA1", "HHAT", "WNT10A", "H3C10", "NRP1", "KCNH8", "SCAMP5",
  "FAM193B", "TRBV5-1", "OSBPL5", "ST3GAL2", "DDIT4", "PHEX", "KCNK10", "CLIP2", "ARPP21", "ACRBP",
  "KRT5", "PADI4", "SH3D19", "NFE2", "EPDR1", "SAMSN1", "CMKLR1", "TUBB6", "FOS", "PYHIN1",
  "SLC1A1", "SLC7A11", "CYP46A1", "S100A9", "TUBB1", "FLNB", "CAMKK1", "BTNL9", "LHFPL2", "TGFBR3",
  "NAT8L", "ARL5C", "GNAZ", "KCTD5", "GP9", "IGKV1-5", "CBLB", "CREB5", "TBC1D4", "RBP7",
  "TRPC1", "PTK7", "NTM", "C19orf33", "ANKRD36", "HLCS", "PHLPP2", "ESAM", "DAPK2", "PTGDS",
  "PLVAP", "ALDH1A1", "CDA", "IL1B", "C12orf42", "ZNRF2", "IL3RA", "NFIL3", "ENHO", "CCL8",
  "CLEC4C", "AJAP1", "BAIAP2L2", "NAMPT", "ENSG00000281593", "BAIAP2L1", "PLXNA4", "SH2D4A", "RGS7", "TPM1",
  "ALOX12", "ACSL1", "E2F5", "HBD", "PAPLN", "APOO", "CCDC50", "MOXD1", "CH25H", "CLU",
  "IGHV3-30", "RNF165", "IRAK3", "SPDYE21", "SDC2", "AQP9", "UGCG", "CIP2A", "C1orf74", "L2HGDH",
  "NIBAN3", "PER1", "ST6GALNAC4", "ITGB3", "RUBCN", "DOK3", "FAM160A1", "MAP1A", "IFITM3", "GPM6B",
  "SLC35F3", "RPS6KA2", "ST6GALNAC3", "SPRY3", "LAIR1", "PRKAR2B", "RCBTB2", "GRAMD1B", "MCEMP1", "AFF1",
  "ABHD6", "MAPKAPK2", "MZB1", "HS3ST3B1", "TCF4", "IGHV1-46", "NPC1", "GREM2", "MSRB3", "BSPRY",
  "SEL1L3", "APP", "WDFY4", "VCAN", "MED12L", "TRBV5-6", "GSG1L", "PIWIL4", "P2RY6", "PARP10",
  "PTPRS", "ENPP4", "KIRREL3", "PIP5K1A", "SLC25A10", "PAFAH2", "LCN2", "CARD11", "TRBV9", "IGLV1-51",
  "RASGEF1B", "ETS2", "TREML1", "RHEX", "CD14", "HS3ST1", "FARP2", "CFAP44", "NET1", "VASH2",
  "NR4A3", "APOBEC3H", "TMEM40", "BCL11A", "KCNK1", "RMC1", "MPP1", "GLCCI1", "ARAP2", "KRT7",
  "JAM3", "ZSCAN5A", "DRD4", "NFX1", "FCER1A", "IGF2BP3", "SPIB", "GLMN", "MRFAP1", "PGRMC1",
  "FRS3", "CXCL8", "NFKBIA", "EPHB1", "FAM110A", "GLCE", "AMIGO2", "MCC", "SYCP2L", "LAPTM4B",
  "ZNF827", "CDH8", "PACSIN1", "AGAP5", "GDPGP1", "SIT1", "ITCH", "SLC15A4", "HACE1", "RUBCNL",
  "ZNF521", "SEC14L5", "NHSL1", "DNASE1L3", "CD247", "MPP6", "RREB1", "SYNRG", "MAP1LC3B2", "IGHV3-64D",
  "ERN1", "NECTIN1", "KLHL13", "MAOB", "USP24", "PPIL1", "LTK", "ZFP36L1", "ZNF500", "LTB",
  "PTGS2", "PLPP1", "SMC6", "GLT1D1", "TNFRSF21", "VEPH1", "LGMN", "GGA2", "SFMBT2", "SLAMF6",
  "LPL", "LIN52", "GABPB2", "LIN7A", "CCL3", "VEGFB", "TMEM91", "NEFM", "ENSG00000259753", "RNF170",
  "ZHX2", "IGHV3-7", "SPC24", "MPIG6B", "FMNL3", "CCDC138", "TENM4", "AGBL2", "TRDMT1", "ARRDC5",
  "GFI1B", "TSPAN33", "HOXA5", "HP", "P3H2", "CPED1", "POLE", "PPM1J", "SLC30A1", "PRF1",
  "CDH1", "TDRD1", "GLI3", "PLEKHD1", "ERO1B", "TSC22D3", "ODAD3", "AHI1", "CEBPD", "LYPD3",
  "H2BC11", "IGF2R", "AAR2", "SKAP1", "SLC9C2", "ADGRG5", "PDE4A", "SMARCD3", "CDR2", "CCDC159",
  "CCDC78", "SLC12A2", "RNASE2", "MS4A4E", "SBF1", "EDIL3", "NAA25", "IL18R1", "LYG1", "ELF3",
  "C18orf25", "FBH1", "DAAM1", "CCNB2", "C5orf63", "AVEN", "HOXA10", "ENSG00000124593", "OGT", "MITF",
  "GCC1", "ZNF124", "OPHN1", "LLGL2", "RAB38", "IGKV3D-20", "FARP1", "PCED1B", "HYAL1", "S100A3",
  "GPR39", "TMEM41B", "SLC24A4", "CXCL2", "IL13RA1", "ZFAT", "COL23A1", "THSD1", "IGKV3-15", "SMARCAD1",
  "RAB13", "NOTCH4", "SCARB2", "SIK2", "CTSA", "CRIP3", "SCML2", "CTSD", "DUSP1", "ANKRD20A1",
  "FOXO3B", "STAT4", "TUBA4A", "ENSG00000258311", "XKR6", "LARP4B", "ASGR1", "TBCEL", "C14orf93", "CXCL5",
  "NOTCH2NLC", "CHRM3", "MIB2", "C20orf194", "TLR5", "ABHD15", "IRF8", "ENSG00000284946", "OFD1", "NOD1",
  "MROH8", "SEPTIN11", "TP53TG5", "HSPA1A", "MCOLN2", "GOLGA8B", "MAFF", "TNFSF10", "FAM149A", "SAP30",
  "IGKV1D-8", "PROSER3", "EPB41L5", "MDFIC", "PLAC8", "NEK8", "IGFBP7", "GGTA1", "ZFYVE26", "HOXA3",
  "RETN", "MICAL3", "DUSP2", "NOTCH3", "DERL3", "RYR3", "DEAF1", "ENSG00000266086", "OXTR", "SH3BGRL2",
  "GPR52", "PAG1", "ZNF334", "MRNIP", "PBX3", "SGCZ", "KDM4B", "TRIM14", "CTDSPL", "CIB2",
  "ZMYM2", "FKBP14", "SPARC", "DHX57", "SEC61A2", "PMAIP1", "C12orf75", "CFAP299", "CLCN5", "TBC1D32",
  "TULP4", "CRB1", "TPM2", "KLRB1", "PAPPA", "SMPD3", "SLITRK5", "GAA", "RBFOX2", "HSD17B13",
  "SLC23A2", "PROC", "GUCY1B1", "CYSLTR1", "ENDOD1", "MXD1", "MPP3", "STIP1", "IGHV5-51", "MTMR1",
  "ADGRG1", "RGS2", "NAV2", "TUBE1", "ASGR2", "RAB27B", "SP2", "ETS1", "ESR2", "IQSEC3",
  "TCF19", "RNF115", "CALM3", "CASZ1", "NUDT12", "YBX3", "COL9A2", "SAMD14", "CCDC69", "TUNAR",
  "ZNF407", "MIS12", "SEMA5A", "CLIC5", "PF4V1", "ROR1", "ZNF17", "ZC2HC1C", "PTAFR", "NUDT17",
  "TLCD3A", "MDM2", "ADNP", "NR2F2", "SCAMP4", "ABCF3", "METTL8", "CENPE", "PRRT4", "AGAP4",
  "F13A1", "NOL3", "FPR1", "NPL", "ATP6V0A2", "SLC41A2", "PIBF1", "RAB11FIP5", "PF4", "ZBTB34"
)

# 3. Subset the data for only the requested genes
subset_data <- data[data$Gene %in% gene_list,]

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Define patient order
patient_order <- c(18, 12, 14, 7, 5, 20, 21, 19, 2, 3, 13)

# Create column names based on the patient order
ordered_columns <- paste0("patient_", patient_order, "_log2FC")

# Check that all ordered columns exist in the data
# If any column doesn't exist, you might need to verify the column naming in the CSV.
all(ordered_columns %in% colnames(subset_data))

# columns that exist
ordered_columns <- intersect(ordered_columns, colnames(subset_data))

# Reorder the data columns to match the desired patient order
subset_data <- subset_data[, c("Gene", ordered_columns)]

# Set rownames to Gene, remove Gene column for heatmap matrix
rownames(subset_data) <- subset_data$Gene
subset_data <- subset_data[,-1]  # remove the Gene column

# Convert to a numeric matrix if needed
mat <- as.matrix(subset_data)

patient_ids <- patient_order                              # same order you plotted
group_tag   <- ifelse(patient_ids %in% short_term_survivor_group,
                      "short term", "long term")
new_labels  <- paste0("patient_", patient_ids, " (", group_tag, ")")

colnames(mat) <- new_labels        # replace the matrix column names

# 4. Plot the heatmap
# You can adjust clustering, scaling, color palettes, etc. as needed.
# Run pheatmap and store the result in an object
res <- pheatmap(mat, 
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "row",
                show_rownames = TRUE, 
                show_colnames = TRUE,
                main = "Heatmap of Selected Genes",
                fontsize_row = 6,
                fontsize_col = 8)

# Extract the hierarchical clustering tree for rows
hc_rows <- res$tree_row

# Choose the number of clusters (e.g., k=5). This number can be adjusted based on 
# the structure you see in the dendrogram.
k <- 20

# Cut the dendrogram into k clusters
clusters <- cutree(hc_rows, k = k)

# Create a data frame with gene and cluster assignments
gene_clusters <- data.frame(
  Gene = rownames(mat),
  Cluster = clusters
)

# Write the cluster assignments to a CSV file so you can inspect them outside R
write.csv(gene_clusters, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_heatmap_gene_clusters.csv", row.names = FALSE)





########################################################################################################################################################################
# trying to see if there is any correlation between innate and adaptive immunity
########################################################################################################################################################################

# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# # Read the log2FC data
# data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_normalized_log2FC_C1_vs_Pre.csv",
#                  header = TRUE,
#                  stringsAsFactors = FALSE)

################################################################################
# 1. Get the 500 most-significant genes (by FDR) from the DE results
################################################################################
de_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_differential_expression_results.csv"

de_results <- read.csv(de_file, header = TRUE, stringsAsFactors = FALSE)

sig500 <- de_results %>% 
  arrange(FDR) %>%                   # lowest FDR first
  slice_head(n = 500)                # take top 500 rows

gene_list_500   <- sig500$Gene
gene_list_up    <- sig500 %>% filter(logFC  > 0) %>% pull(Gene)
gene_list_down  <- sig500 %>% filter(logFC  < 0) %>% pull(Gene)

################################################################################
# 2. Helper to collapse a gene set to a per-patient activity
################################################################################
get_activity <- function(expr_df, genes_vec, new_name){
  expr_df %>%
    filter(Gene %in% genes_vec) %>%            # keep only the gene set
    select(-Gene) %>%                          # drop gene column
    pivot_longer(cols = everything(),
                 names_to  = "patient_col",
                 values_to = "log2FC") %>%
    mutate(patient_id = str_extract(patient_col, "\\d+")) %>%
    group_by(patient_id) %>%
    summarize("{new_name}" := mean(log2FC, na.rm = TRUE),
              .groups = "drop")
}

pathway_all  <- get_activity(data, gene_list_500,  "pathway_activity_all")
pathway_up   <- get_activity(data, gene_list_up,   "pathway_activity_up")
pathway_down <- get_activity(data, gene_list_down, "pathway_activity_down")

################################################################################
# 3. Add the three pathway-activity variables to the Cox input
################################################################################
cox_input_df <- cox_input_df %>%
  left_join(pathway_all,  by = "patient_id") %>%
  left_join(pathway_up,   by = "patient_id") %>%
  left_join(pathway_down, by = "patient_id")

# cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9), ]
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9, 1, 4), ]



clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_All_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion$patient_id <- rownames(clonal_expansion)
merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity_all", "pathway_activity_up", "pathway_activity_down", "Arm", "MGMT", "IDH", "Sex", "Age")], clonal_expansion, by = "patient_id")

merged_df <- merged_df[!is.na(merged_df$pathway_activity_all) & !is.na(merged_df$clonal_expansion), ]

# spearman_value <- cor(merged_df$pathway_activity,
#                       merged_df$clonal_expansion,
#                       method = "spearman")
# 
# spearman_value

merged_df[merged_df$IDH == "UNK", "IDH"] = "NEG"
merged_df

#############################################################################################################
# gemini correlation metrics
#############################################################################################################
# Load the mccr package
library(mccr)
# Initialize list to store results
results <- list()

# 2. Both columns numeric
results$'Pearson (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "pearson")
results$'Spearman (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "spearman")
results$'Kendall (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "kendall")

# 3. Numeric pathway_activity vs. Binary clonal_expansion (> 1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
# Check for variance in the binary column
if (length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_binary_gt1, method = "pearson")
} else {
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- NA
  warning("Binary clonal expansion column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 4. Binary pathway_activity (> median) vs. Numeric clonal_expansion
median_pathway <- median(merged_df$pathway_activity_all)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity_all > median_pathway)
# Check for variance in the binary column
if (length(unique(merged_df$pathway_binary_median)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson")
} else {
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- NA
  warning("Binary pathway activity column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 5. Both columns binary
# Check for variance in both binary columns
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Phi coefficient (equivalent to Pearson on 0/1 binary data)
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  
  # Matthews Correlation Coefficient (MCC) using the mccr package
  # Note: mccr expects factors or numeric {0, 1} or {-1, 1}. Our 0/1 format works.
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
  
} else {
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  warning("One or both binary columns lack variance. Cannot calculate Phi or MCC.")
}


# Find the metric with the most negative correlation
best_metric_name <- NULL
highest_inverse_corr <- 0 # Initialize to 0, looking for most negative

for (metric_name in names(results)) {
  value <- results[[metric_name]]
  if (!is.na(value) && value < highest_inverse_corr) {
    highest_inverse_corr <- value
    best_metric_name <- metric_name
  }
}

# Print the results
cat("Correlation Results:\n\n")
for (metric_name in names(results)) {
  cat(sprintf("%s: %.4f\n", metric_name, results[[metric_name]]))
}

cat("\n--------------------------------------------------\n")
cat("Best Inverse Correlation Found:\n")
cat(sprintf("Method: %s\n", best_metric_name))
cat(sprintf("Correlation Value: %.4f\n", highest_inverse_corr))
cat("--------------------------------------------------\n")

# Display the median used for pathway activity binarization
cat(sprintf("\nMedian used for 'pathway_activity' binarization: %.4f\n", median_pathway))



###############################################################################################################
# applying linear model
###############################################################################################################

# 1) Make sure Sex is treated as a factor (categorical) -----------------------
merged_df$Sex <- factor(merged_df$Sex)      # ‘F’ will be the reference level by default

# 2) Fit the model ------------------------------------------------------------
fit <- lm(clonal_expansion ~ pathway_activity_down + Sex + Age, data = merged_df)

# 3) Inspect the results ------------------------------------------------------
summary(fit)            # classic R output
#– or –#
library(broom)
tidy(fit)               # neat tibble with estimates, SEs and p-values
###############################################

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# 2. Add a grouping column ----------------------------------------------------------
merged_df$Group <- dplyr::case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                          ~ "Uncategorised"
)

# 3. Plot with ggplot2 --------------------------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = pathway_activity_down,
                      y = clonal_expansion,
                      colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Pathway activity vs. clonal expansion",
       subtitle = "Points coloured by survivor group",
       x = "Pathway activity",
       y = "Clonal expansion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor"  = "#4DBBD5",
    "Uncategorised"       = "grey60"
  ))

















################################################################################################################################
# starting with the DC Analysis (CCR2 Positive DC Cells)
################################################################################################################################

############################################################################################################################################################################################
# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)

setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis")

# 1. Load the Seurat object
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_DC_CCR2_Pos_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")

# Check the metadata
metadata <- seurat_obj@meta.data
head(metadata)

# Ensure that 'TimePoint' and 'Patient' columns exist
if(!all(c("TimePoint", "Patient") %in% colnames(metadata))) {
  stop("The Seurat object does not contain 'TimePoint' and/or 'Patient' columns in metadata.")
}

# 2. Prepare pseudo-bulk RNA-seq counts
# Aggregate counts by Patient and TimePoint

# Extract the raw counts matrix
# Commonly in the "RNA" assay
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Add Patient and TimePoint information to the metadata
metadata <- seurat_obj@meta.data %>%
  select(Patient, TimePoint)

# **Modified Step:** Create a unique identifier for each pseudo-bulk sample with "patient_" prefix
metadata$Sample <- paste("patient", metadata$Patient, metadata$TimePoint, sep = "_")
# Example: "patient_1_Pre", "patient_1_C1", "patient_2_Pre", etc.

# Sum the counts for cells belonging to the same Sample
# Transpose counts to have cells as rows and genes as columns for aggregation
counts_df <- as.data.frame(t(as.matrix(counts)))
counts_df$Sample <- metadata$Sample

# Aggregate counts by Sample
pseudo_bulk_counts <- counts_df %>%
  group_by(Sample) %>%
  summarise_all(sum)

# **Corrected Step:** Convert back to a matrix with genes as rows and samples as columns
# Since 'pseudo_bulk_counts' has samples as rows and genes as columns,
# we need to transpose it to have genes as rows and samples as columns.

# Remove the 'Sample' column and transpose
pseudo_bulk_counts_mat <- t(as.matrix(pseudo_bulk_counts[,-1]))

# Assign gene names as rownames
rownames(pseudo_bulk_counts_mat) <- rownames(counts)  # Correct: Gene names

# Assign sample names as column names
colnames(pseudo_bulk_counts_mat) <- pseudo_bulk_counts$Sample

# Verify dimensions
dim(pseudo_bulk_counts_mat)  # Should be genes x samples (e.g., 19343 x 24)

write.csv(pseudo_bulk_counts_mat, file = "exp_DC_CCR2_Pos_Cells_pseudo_bulk_RNASeq.csv")

# 3. Prepare the design matrix for edgeR
# Extract Patient and TimePoint information for each Sample
sample_info <- data.frame(
  Sample = colnames(pseudo_bulk_counts_mat)
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_")

# **Modified Step:** Reconstruct the Patient identifier with the "patient_" prefix
sample_info$Patient <- paste(sample_info$Patient_Prefix, sample_info$Patient_ID, sep = "_")
# Example: "patient_1", "patient_2", etc.

# Ensure factors are correctly set
sample_info$Patient <- factor(sample_info$Patient)
sample_info$TimePoint <- factor(sample_info$TimePoint, levels = c("Pre", "C1"))

# 4. Create DGEList object
dge <- DGEList(counts = pseudo_bulk_counts_mat)

# 5. Normalize the data
dge <- calcNormFactors(dge)

# 6. Create the design matrix
# Model: ~ Patient + TimePoint
design <- model.matrix(~ Patient + TimePoint, data = sample_info)

# 7. Estimate dispersion
dge <- estimateDisp(dge, design)

# 8. Fit the model and perform differential expression
fit <- glmFit(dge, design)
# Assuming "TimePointC1" is the coefficient for C1 vs Pre
lrt <- glmLRT(fit, coef = "TimePointC1")

# 9. Extract differential expression results
de_results <- topTags(lrt, n = Inf)$table
de_results <- as.data.frame(de_results)
# Add gene names as a column
de_results$Gene <- rownames(de_results)

# 10. Save the differential expression results
write.csv(de_results, file = "exp_DC_CCR2_Pos_Cells_differential_expression_results.csv", row.names = FALSE)

# 11. Create input file for GSEA
# Typically, GSEA requires a ranked list of genes. One common ranking metric is logFC multiplied by the sign of the p-value.

# Here, we'll use the signed -log10(p-value) multiplied by the sign of logFC
# Alternatively, you can choose other ranking metrics based on your preference

# To handle p-values of 0, add a small pseudocount
de_results$PValue[de_results$PValue == 0] <- 1e-300

de_results$RankingMetric <- -log10(de_results$PValue) * sign(de_results$logFC)

# Order the genes by the ranking metric
ranked_genes <- de_results %>%
  arrange(desc(RankingMetric)) %>%
  select(Gene, RankingMetric)

# Save the ranked list for GSEA
# **Option 1:** GCT format
# GCT format requires specific headers. Here's how to structure it:

# Create a GCT header
gct_header <- data.frame(
  NAME = c("#1.2"),
  DESCRIPTION = c("RankingMetric"),
  stringsAsFactors = FALSE
)

# Combine header and data
gct_data <- data.frame(
  NAME = ranked_genes$Gene,
  DESCRIPTION = ranked_genes$Gene,  # GSEA often expects gene symbols or descriptions here
  RankingMetric = ranked_genes$RankingMetric,
  stringsAsFactors = FALSE
)

# Write the GCT file
write.table(gct_header, file = "exp_DC_CCR2_Pos_Cells_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gct_data, file = "exp_DC_CCR2_Pos_Cells_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE)

# **Option 2:** RNK format
# RNK format is simpler and often preferred for ranked lists

write.table(ranked_genes, file = "exp_DC_CCR2_Pos_Cells_GSEA_ranked_genes.rnk", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Inform the user that the process is complete
cat("Pseudo-bulk RNA-seq aggregation, differential expression analysis, and GSEA input file creation are complete.\n")
cat("Differential expression results saved to 'differential_expression_results.csv'.\n")
cat("GSEA ranked gene list saved to 'GSEA_ranked_genes.gct' and 'GSEA_ranked_genes.rnk'.\n")





######################################################################################
library(dplyr)
library(tidyr)
# 1. Read the Pseudo-Bulk Counts Matrix
# -------------------------------------

# Specify the path to your pseudo_bulk_counts_mat.csv
counts_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_pseudo_bulk_RNASeq.csv"

# Read the CSV file
# Assuming the first column contains gene names and is set as rownames
counts_df <- read.csv(counts_csv_path, row.names = 1, check.names = FALSE)

# Verify the dimensions and a snippet of the data
cat("Dimensions of the counts matrix:", dim(counts_df), "\n")
head(counts_df)


# 4. Normalize Counts Using "Pre" as Baseline
# -------------------------------------------

# Subset the counts matrix to include only selected genes
selected_counts <- counts_df

# Parse sample names to extract Patient and TimePoint information
# Assuming sample names are in the format "patient_X_TimePoint", e.g., "patient_1_Pre"
sample_names <- colnames(selected_counts)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = sample_names,
  stringsAsFactors = FALSE
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_", remove = FALSE) %>%
  mutate(Patient = paste(Patient_Prefix, Patient_ID, sep = "_")) %>%
  select(Sample, Patient, TimePoint)

# View sample information
print(sample_info)

# Get a list of unique patients
patients <- unique(sample_info$Patient)
cat("Number of unique patients:", length(patients), "\n")

# Initialize a list to store log2 fold changes
log2fc_list <- list()

# Loop through each patient to compute log2FC
for (patient in patients) {
  # Identify "Pre" and "C1" samples for the patient
  pre_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "Pre") %>%
    pull(Sample)
  
  c1_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "C1") %>%
    pull(Sample)
  
  # Check if both samples are present
  if (length(pre_sample) == 1 & length(c1_sample) == 1) {
    # Extract counts for "Pre" and "C1" samples
    pre_counts <- selected_counts[, pre_sample]
    c1_counts <- selected_counts[, c1_sample]
    
    # Compute log2 fold change with a pseudocount of 1 to avoid division by zero
    log2fc <- log2((c1_counts + 1) / (pre_counts + 1))
    
    # Ensure that log2fc is a named vector with gene names
    names(log2fc) <- rownames(selected_counts)
    
    # Store in the list
    log2fc_list[[patient]] <- log2fc
  } else {
    warning(paste("Missing 'Pre' or 'C1' sample for patient:", patient))
  }
}

# Combine the list into a normalized matrix
# Rows: Genes, Columns: Patients (log2FC)
normalized_mat <- do.call(cbind, log2fc_list)

# Rename columns to indicate log2FC
colnames(normalized_mat) <- paste0(colnames(normalized_mat), "_log2FC")

# Inspect the normalized matrix
dim(normalized_mat)
head(normalized_mat)

# 5. Scale the Normalized Data
# ----------------------------

# Scaling can be done per gene to standardize the log2FC across patients
# This is useful for visualization purposes

# Transpose, scale, and transpose back
scaled_mat <- t(scale(t(normalized_mat)))

# Replace any NA values resulting from scaling with 0
scaled_mat[is.na(scaled_mat)] <- 0

# Inspect the scaled matrix
dim(scaled_mat)
head(scaled_mat)

# 6. Save the Normalized and Scaled Data as CSV
# ---------------------------------------------

# Prepare the normalized data for saving
normalized_df <- as.data.frame(normalized_mat)
normalized_df$Gene <- rownames(normalized_df)

# Reorder columns to have 'Gene' as the first column
normalized_df <- normalized_df %>%
  select(Gene, everything())

# Save the normalized log2FC data
write.csv(normalized_df, file = "exp_DC_CCR2_Pos_Cells_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Normalized log2 fold change data saved to 'exp_DC_normalized_log2FC_C1_vs_Pre.csv'.\n")

# Prepare the scaled data for saving
scaled_df <- as.data.frame(scaled_mat)
scaled_df$Gene <- rownames(scaled_df)

# Reorder columns to have 'Gene' as the first column
scaled_df <- scaled_df %>%
  select(Gene, everything())

# Save the scaled normalized data
write.csv(scaled_df, file = "exp_DC_CCR2_Pos_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Scaled normalized data saved to 'exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv'.\n")



#########################################################################################################################################################
# Plotting heatmap

# Install necessary packages if not installed
# install.packages("pheatmap")

library(pheatmap)

# 1. Read the file
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv", header = TRUE, stringsAsFactors = FALSE)

# The columns are expected to be:
# Gene, patient_10_log2FC, patient_12_log2FC, patient_13_log2FC, patient_14_log2FC,
# patient_18_log2FC, patient_19_log2FC, patient_20_log2FC, patient_21_log2FC,
# patient_2_log2FC, patient_3_log2FC, patient_5_log2FC, patient_7_log2FC

# 2. Define the gene list
gene_list <- c(
  "IGKV3-11", "AREG", "HBB", "COL26A1", "CAVIN2", "NTM", "TUBB1", "PTGDS", "EPHB1", "S100A12",
  "RGS7", "PHEX", "H3C10", "SEMA5A", "CUX2", "SETBP1", "NRP1", "GNLY", "IGHV3-23", "S100A8",
  "DAPK2", "KCTD5", "CARD11", "MYL9", "SLCO5A1", "P3H2", "GP9", "PRF1", "CLEC1B", "TBC1D4",
  "IGKV4-1", "COBLL1", "FPR1", "EREG", "SAMSN1", "SLC35F3", "COL24A1", "PER1", "DDIT4", "FOS",
  "PLXNA4", "CMTM5", "TSC22D3", "PF4", "TBC1D32", "TMEM40", "STAT4", "TUBB6", "ALOX12", "SH2D4A",
  "NRGN", "KCNH8", "CD247", "IGHV4-34", "HFM1", "NIBAN3", "THBS1", "PPBP", "S100A9", "TNFRSF21",
  "IGLV1-51", "ARAP2", "FGFBP2", "DERL3", "ST3GAL2", "BCL11A", "ZHX2", "STXBP1", "EGR1", "AUTS2",
  "MAP1A", "TAL1", "CLIP2", "SLC24A4", "FLNB", "IGKV1-5", "VCAN", "PLVAP", "MCC", "MOXD1",
  "PRKAR2B", "ACRBP", "MMD", "LGMN", "L2HGDH", "NAMPT", "BAIAP2L2", "RETN", "KRT5", "PPP3CC",
  "ZFAT", "TAMALIN", "HHAT", "CD9", "GRAMD1B", "TREML1", "WNT2B", "TRAF4", "TCL1B", "MNDA",
  "ARID3B", "LHFPL2", "FAM160A1", "SIDT1", "CRYAA2", "ARHGAP5", "PIP5K1A", "ACSL1", "CLEC4C", "CCDC186",
  "JCHAIN", "IGKC", "PTPRS", "ATXN1", "NIBAN2", "TDRD1", "BORCS6", "PELI1", "ENSG00000281593", "SLC1A1",
  "STRBP", "AK5", "CTSA", "TRBV5-1", "KHDRBS2", "GNG11", "PFKFB2", "MZB1", "RMC1", "SPARC",
  "RHEX", "ST6GALNAC4", "ITGA2B", "CD3E", "IL3RA", "ERN1", "LRRC36", "APOO", "C5AR1", "DYSF",
  "ANKRD28", "GLDN", "PLXDC2", "ESAM", "EXT1", "TRBV4-2", "PYHIN1", "APBA2", "PLCB4", "USP24",
  "IFITM3", "PGAP1", "GFI1B", "LIMS1", "LAMP5", "H2AC6", "ETS2", "PPP2R2B", "POLD1", "PDXK",
  "MCEMP1", "UGCG", "CLU", "TLR3", "CREB5", "RREB1", "NR4A1", "KRT7", "RAB20", "SCAI",
  "PRSS23", "NFIL3", "NPC1", "GLUL", "ZNF778", "GLCE", "CD14", "MAP1LC3B2", "DNASE1L3", "TRIM37",
  "ZSCAN16", "FGFR1", "ADAM19", "IGF1R", "FHIT", "RUBCNL", "USP11", "FCGR3A", "GLCCI1", "NTNG2",
  "CDR2", "ENSG00000196826", "DACH1", "TRBV7-9", "GREM2", "C14orf93", "TPM2", "CAPS", "CCDC50", "GZMB",
  "TEP1", "AJAP1", "IL7R", "ADCY6", "CDKN1A", "KLHL13", "ITGB3", "FRMD3", "SLC44A5", "LYPD3",
  "PLOD2", "PBX3", "IL32", "KNOP1", "ENHO", "GUCY1B1", "RABGAP1L", "NAA25", "ABHD8", "CYP46A1",
  "PDE4B", "ATP11A", "MFNG", "AQP9", "RPL36A", "RCBTB2", "SERPINA1", "MAP4K5", "COG6", "TRBV4-1",
  "DUSP1", "MYBL2", "LRRC26", "RASGEF1B", "ERO1B", "NOL3", "CCDC91", "HEY1", "LCN2", "NEIL1",
  "GRAP2", "PMEPA1", "MOSMO", "ASPM", "WDFY4", "FBXO9", "GPR39", "H2BC4", "TMEM101", "KLRC1",
  "DDX21", "CCDC138", "P2RY14", "TTC39C", "VEGFB", "DNAJB5", "IRAK3", "KIF4A", "HYAL1", "PRKD3",
  "TP53TG5", "MCL1", "GATA2", "TP53I3", "IL1B", "TCL1A", "LEF1", "PACS2", "UPP1", "THAP8",
  "SNURF", "WWOX", "AHI1", "ETS1", "PGM3", "ZNF398", "IFITM2", "GADD45A", "PTPRM", "RNLS",
  "TSC22D1", "RNF165", "APPBP2", "CUEDC1", "CSTA", "KRBA1", "SOD2", "CLCN5", "KBTBD3", "SPDYE21",
  "SSX2IP", "AGAP6", "CCDC88C", "F13A1", "GAS2L1", "LIN7A", "GATM", "GPATCH8", "NUDT12", "ZNF610",
  "TRAV39", "C2orf88", "LINC02210-CRHR1", "TSPAN33", "DCAF5", "SERPINH1", "FAM214B", "ZDHHC17", "LAIR1", "IRS2",
  "SH3D19", "PIF1", "GZMM", "ZNF189", "RAE1", "ARNTL2", "PILRA", "ENSG00000124593", "SKAP1", "C9orf116",
  "RBFOX2", "ZNF444", "LILRA4", "APH1B", "ELP2", "PCED1B", "GRIP1", "MIB2", "C5orf63", "SPHK1",
  "CDK1", "TRBV24-1", "LTBP1", "ZNF521", "CTTN", "TLN2", "MYLK", "MAP3K3", "GSTP1", "PIK3IP1",
  "HBD", "LRG1", "F5", "LDHD", "KIF13A", "SCAMP4", "ABCF3", "CCND2", "RNF115", "SIGLEC6",
  "ABLIM1", "ZSCAN5A", "SAP30", "CHST2", "SYNE1", "APOBR", "IGHV3-7", "RAPGEFL1", "ST18", "GIN1",
  "NKG7", "NADK2", "CD96", "PDE4A", "PGRMC1", "METTL8", "AIFM2", "ANKRD36", "MRPL49", "FCRL6",
  "TPM1", "FCER1A", "SMS", "TMEM91", "RAB34", "EPDR1", "IFNGR2", "CYB561D1", "IGSF6", "CCNB2",
  "TRBV9", "PLXNA2", "EFCAB2", "SOX6", "NDUFB2", "NIPA1", "SYNM", "SULT1A1", "HLA-DRB5", "RERE",
  "RORA", "KIR3DL1", "GEM", "MRE11", "GAA", "SPIB", "NUDT1", "NPIPB2", "FYCO1", "LY86",
  "S100A6", "MBOAT2", "RARA", "SLC9A9", "H1-0", "ST6GALNAC3", "ATP5PF", "TMEM63B", "RBP7", "IGFBP3",
  "UGGT2", "SLC45A4", "KCMF1", "PLCB2", "PDE6G", "PPP1R14A", "GOLM1", "YWHAH", "CASP1", "SCML2",
  "NPL", "ERVK3-1", "TTC12", "ENSG00000272921", "GLS", "ARID3A", "UBE2D1", "MFSD9", "NOP16", "C11orf71",
  "CEP128", "PROSER3", "SLC7A11", "CDA", "ZFP62", "MPP7", "TSC22D2", "ARL5C", "ARHGEF26", "MRC2",
  "PEX7", "TRIM59", "GSDMB", "RIC3", "GATA1", "ADGRB3", "PARN", "SIGMAR1", "FARP1", "ADM",
  "NOXA1", "SLC14A2", "LIPN", "ENSG00000269237", "KLRB1", "TRMT10B", "JAZF1", "PHLPP2", "PLEKHG2", "TNS2",
  "FGGY", "GPM6B", "ZSCAN20", "BICD1", "FBXO34", "CAMK4", "LST1", "CD83", "SLC40A1", "TSTD1",
  "IGLV4-69", "UBAC1", "SDC2", "DIS3L2", "CCDC141", "CHPF2", "CTTNBP2NL", "PPARA", "IL2RB", "TLCD3A",
  "ENSG00000263620", "SBF1", "PIM3", "SLC30A1", "ENSG00000273748", "CD244", "ATP5F1D", "DAAM1", "DCAF13", "VPS26A",
  "PPM1A", "DNAAF2", "PTPN12", "ABHD6", "MESP1", "ENSG00000283782", "CTSD", "METTL3", "CRIM1", "CST7"
)

# 3. Subset the data for only the requested genes
subset_data <- data[data$Gene %in% gene_list,]

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Define patient order
patient_order <- c(14, 7, 5, 20, 21, 19, 2, 3, 13)

# Create column names based on the patient order
ordered_columns <- paste0("patient_", patient_order, "_log2FC")

# Check that all ordered columns exist in the data
# If any column doesn't exist, you might need to verify the column naming in the CSV.
all(ordered_columns %in% colnames(subset_data))

# columns that exist
ordered_columns <- intersect(ordered_columns, colnames(subset_data))

# Reorder the data columns to match the desired patient order
subset_data <- subset_data[, c("Gene", ordered_columns)]

# Set rownames to Gene, remove Gene column for heatmap matrix
rownames(subset_data) <- subset_data$Gene
subset_data <- subset_data[,-1]  # remove the Gene column

# Convert to a numeric matrix if needed
mat <- as.matrix(subset_data)

patient_ids <- patient_order                              # same order you plotted
group_tag   <- ifelse(patient_ids %in% short_term_survivor_group,
                      "short term", "long term")
new_labels  <- paste0("patient_", patient_ids, " (", group_tag, ")")

colnames(mat) <- new_labels        # replace the matrix column names

# 4. Plot the heatmap
# You can adjust clustering, scaling, color palettes, etc. as needed.
# Run pheatmap and store the result in an object
res <- pheatmap(mat, 
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "row",
                show_rownames = TRUE, 
                show_colnames = TRUE,
                main = "Heatmap of Selected Genes",
                fontsize_row = 6,
                fontsize_col = 8)

# Extract the hierarchical clustering tree for rows
hc_rows <- res$tree_row

# Choose the number of clusters (e.g., k=5). This number can be adjusted based on 
# the structure you see in the dendrogram.
k <- 20

# Cut the dendrogram into k clusters
clusters <- cutree(hc_rows, k = k)

# Create a data frame with gene and cluster assignments
gene_clusters <- data.frame(
  Gene = rownames(mat),
  Cluster = clusters
)

# Write the cluster assignments to a CSV file so you can inspect them outside R
write.csv(gene_clusters, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_heatmap_gene_clusters.csv", row.names = FALSE)





########################################################################################################################################################################
# trying to see if there is any correlation between innate and adaptive immunity
########################################################################################################################################################################

# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# # Read the log2FC data
# data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_normalized_log2FC_C1_vs_Pre.csv",
#                  header = TRUE,
#                  stringsAsFactors = FALSE)

################################################################################
# 1. Get the 500 most-significant genes (by FDR) from the DE results
################################################################################
de_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_differential_expression_results.csv"

de_results <- read.csv(de_file, header = TRUE, stringsAsFactors = FALSE)

sig500 <- de_results %>% 
  arrange(FDR) %>%                   # lowest FDR first
  slice_head(n = 500)                # take top 500 rows

gene_list_500   <- sig500$Gene
gene_list_up    <- sig500 %>% filter(logFC  > 0) %>% pull(Gene)
gene_list_down  <- sig500 %>% filter(logFC  < 0) %>% pull(Gene)

################################################################################
# 2. Helper to collapse a gene set to a per-patient activity
################################################################################
get_activity <- function(expr_df, genes_vec, new_name){
  expr_df %>%
    filter(Gene %in% genes_vec) %>%            # keep only the gene set
    select(-Gene) %>%                          # drop gene column
    pivot_longer(cols = everything(),
                 names_to  = "patient_col",
                 values_to = "log2FC") %>%
    mutate(patient_id = str_extract(patient_col, "\\d+")) %>%
    group_by(patient_id) %>%
    summarize("{new_name}" := mean(log2FC, na.rm = TRUE),
              .groups = "drop")
}

pathway_all  <- get_activity(data, gene_list_500,  "pathway_activity_all")
pathway_up   <- get_activity(data, gene_list_up,   "pathway_activity_up")
pathway_down <- get_activity(data, gene_list_down, "pathway_activity_down")

################################################################################
# 3. Add the three pathway-activity variables to the Cox input
################################################################################

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF", ]
# cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros


cox_input_df <- cox_input_df %>%
  left_join(pathway_all,  by = "patient_id") %>%
  left_join(pathway_up,   by = "patient_id") %>%
  left_join(pathway_down, by = "patient_id")

# cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9), ]
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9, 1, 4), ]



clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_All_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion$patient_id <- rownames(clonal_expansion)
merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity_all", "pathway_activity_up", "pathway_activity_down", "Arm", "MGMT", "IDH", "Sex", "Age")], clonal_expansion, by = "patient_id")

merged_df <- merged_df[!is.na(merged_df$pathway_activity_all) & !is.na(merged_df$clonal_expansion), ]

# spearman_value <- cor(merged_df$pathway_activity,
#                       merged_df$clonal_expansion,
#                       method = "spearman")
# 
# spearman_value

merged_df[merged_df$IDH == "UNK", "IDH"] = "NEG"
merged_df

#############################################################################################################
# gemini correlation metrics
#############################################################################################################
# Load the mccr package
library(mccr)
# Initialize list to store results
results <- list()

# 2. Both columns numeric
results$'Pearson (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "pearson")
results$'Spearman (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "spearman")
results$'Kendall (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "kendall")

# 3. Numeric pathway_activity vs. Binary clonal_expansion (> 1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
# Check for variance in the binary column
if (length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_binary_gt1, method = "pearson")
} else {
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- NA
  warning("Binary clonal expansion column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 4. Binary pathway_activity (> median) vs. Numeric clonal_expansion
median_pathway <- median(merged_df$pathway_activity_all)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity_all > median_pathway)
# Check for variance in the binary column
if (length(unique(merged_df$pathway_binary_median)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson")
} else {
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- NA
  warning("Binary pathway activity column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 5. Both columns binary
# Check for variance in both binary columns
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Phi coefficient (equivalent to Pearson on 0/1 binary data)
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  
  # Matthews Correlation Coefficient (MCC) using the mccr package
  # Note: mccr expects factors or numeric {0, 1} or {-1, 1}. Our 0/1 format works.
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
  
} else {
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  warning("One or both binary columns lack variance. Cannot calculate Phi or MCC.")
}


# Find the metric with the most negative correlation
best_metric_name <- NULL
highest_inverse_corr <- 0 # Initialize to 0, looking for most negative

for (metric_name in names(results)) {
  value <- results[[metric_name]]
  if (!is.na(value) && value < highest_inverse_corr) {
    highest_inverse_corr <- value
    best_metric_name <- metric_name
  }
}

# Print the results
cat("Correlation Results:\n\n")
for (metric_name in names(results)) {
  cat(sprintf("%s: %.4f\n", metric_name, results[[metric_name]]))
}

cat("\n--------------------------------------------------\n")
cat("Best Inverse Correlation Found:\n")
cat(sprintf("Method: %s\n", best_metric_name))
cat(sprintf("Correlation Value: %.4f\n", highest_inverse_corr))
cat("--------------------------------------------------\n")

# Display the median used for pathway activity binarization
cat(sprintf("\nMedian used for 'pathway_activity' binarization: %.4f\n", median_pathway))



###############################################################################################################
# applying linear model
###############################################################################################################

# 1) Make sure Sex is treated as a factor (categorical) -----------------------
merged_df$Sex <- factor(merged_df$Sex)      # ‘F’ will be the reference level by default

# 2) Fit the model ------------------------------------------------------------
fit <- lm(clonal_expansion ~ pathway_activity_down + Sex + Age, data = merged_df)

# 3) Inspect the results ------------------------------------------------------
summary(fit)            # classic R output
#– or –#
library(broom)
tidy(fit)               # neat tibble with estimates, SEs and p-values
###############################################

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# 2. Add a grouping column ----------------------------------------------------------
merged_df$Group <- dplyr::case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                          ~ "Uncategorised"
)

# 3. Plot with ggplot2 --------------------------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = pathway_activity_down,
                      y = clonal_expansion,
                      colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Pathway activity vs. clonal expansion",
       subtitle = "Points coloured by survivor group",
       x = "Pathway activity",
       y = "Clonal expansion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor"  = "#4DBBD5",
    "Uncategorised"       = "grey60"
  ))























########################################################################################################################################################################
# trying to see if there is any correlation between innate and adaptive immunity (Central Memory CD8 clonal expansion)
########################################################################################################################################################################
# Load necessary libraries
library(tidyverse)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(Rcpp)
library(cowplot)
library(stringr)  # For string manipulation
library(Seurat)


# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_scaled_normalized_log2FC_C1_vs_Pre.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# # Read the log2FC data
# data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_normalized_log2FC_C1_vs_Pre.csv",
#                  header = TRUE,
#                  stringsAsFactors = FALSE)

gene_list_500 <- c(
  "C1QC","C1QB","LYPD2","S100A12","VSIG4","VCAN","SLCO5A1","FOS","S100A8","C1QA","CD163",
  "ADAMTS2","HFM1","VMO1","CEBPD","GPX1","RBP7","AREG","THBS1","S100A9","GNLY","PLD4",
  "IL32","CYP4F22","GZMB","CD247","HES4","TPST1","MS4A6A","PPBP","G0S2","CD14","SCGB3A1",
  "TEX14","ITGA2B","FMN1","FLT3","MYL9","IGHV1-46","NKG7","ACSL1","FKBP5","H3C10","NRGN",
  "CLU","TFEC","TSPAN33","SASH1","CAVIN2","TREML1","RNASE2","LPL","SPON2","CKB","GZMH",
  "PRF1","DDIT4","HSH2D","CDKN1C","HBB","SH2D1B","DOCK4","EGR1","F13A1","IGLV2-14",
  "CTSW","RBFOX2","CST7","PLCB1","IGLV5-45","SLAMF7","LMNA","GZMA","SELPLG","DUSP5",
  "SYTL1","FGFBP2","MPIG6B","PRKN","PF4","SHTN1","ARL4A","PRKAR2B","TSTD1","CMTM5",
  "SPARC","CRIP1","LYZ","STAT4","ARID5B","SPN","IRS2","IL7R","NRG1","STAB1","APOO",
  "ADGRB3","EIF4G3","BLVRB","MTURN","KLRB1","ASGR2","CHST2","BAIAP2L1","PGRMC1","IGKV2-30",
  "MT-ATP8","ETS1","MSR1","OAS1","CSF3R","SLC44A2","C2orf88","BEX3","TRBC1","NFKBIA",
  "SETBP1","HBG1","SLC2A6","RNASE6","CUX1","IL13RA1","SAMSN1","IGLV3-27","SMOX","FRMD4B",
  "ETS2","CD36","CD99","IGLV1-51","TRIM58","PID1","SH3BGRL2","DUSP1","CD7","PIK3CG",
  "HLA-DMB","H2AC6","GBP2","LGALS2","CPED1","TRBC2","FGR","MCEMP1","SIGLEC10","IGKV1-27",
  "IL1R2","TMEM91","CTSC","SCART1","MGST1","GPR183","ITGAM","CCL4","ENTPD1","IL21R",
  "TRBV29-1","MLXIP","KLF12","RHOC","EPB41L3","IGHV4-59","LINGO2","CCR7","SMIM7","LAT2",
  "TNFAIP3","P2RY10","SSBP4","MMD","EPHX3","HCAR3","FAM110A","TMEM40","NAMPT",
  "ENSG00000285589","TRGV9","CD101","HBEGF","CTSL","C1orf167","PER1","LTB","ZNF467","SAP30",
  "RAB39B","CD9","CBLB","CDH23","AHR","PKN1","KLF10","SYAP1","JUN","PIGQ","VAMP8","GNG11",
  "S1PR4","ADAM8","PTPRCAP","RETN","CREB5","WARS1","GP9","SH2D3C","TUBB1","AKT3","CALHM2",
  "LY6E","TSHZ2","MXD1","RGL4","ITGB3","PTK2","ENSG00000259529","APOBEC3G","H2BC11","ACRBP",
  "SIRPA","TMEM128","PRKCQ","SPSB3","UPK3A","KNDC1","TSC22D1","RPLP0","NXF3","ADORA2A",
  "IFITM1","NCF1","STON2","ZAP70","UNC119","EEF1G","PDGFA","SH3BP4","GFOD1","IFIT2",
  "SNED1","CEBPA","S100P","FN1","TMEM176B","LPCAT3","IGHV2-70D","ITGAL","TTC39B",
  "IGLV4-60","RPL41","CLIC3","SNCA","LCK","CD300LF","CCR2","HOMER3","LY9","CSF1R","LGALSL",
  "LAPTM4B","POLD4","PTGDR2","GATD3B","ANPEP","SYNGR2","SLIT2","AP1S2","NT5M","TRBV6-6",
  "PPM1N","ARAP2","LFNG","LMBRD1","GBP5","TTC9","DRAP1","TP53I11","ADGRE2","ZDHHC1","RHOF",
  "MAL","COMT","TCIRG1","TRBV12-3","MYBL1","IGKV3D-11","ALOX12","GPRC5C","CD3E","ARHGAP24",
  "PPIF","TMEM120A","PTGS1","LDHD","CX3CR1","MEN1","DHPS","ZC2HC1A","TBC1D10C","TRDC",
  "ZNF532","STXBP2","ACAP1","TESC","IGLV2-5","YWHAG","PTPRE","CD244","S1PR5","SLC2A9",
  "UCP2","RNASE1","ODAD3","APOE","MCM6","H3C2","PEX5","PADI4","PLA2G4A","CAMTA1",
  "PLA2G12A","CLIC4","PLEKHA7","MARVELD1","SAV1","OBSCN","SIGLEC11","SPIDR","SLC38A7",
  "IGKV4-1","C1RL","ERAP1","ARL6IP5","TBCD","DNM3","CAV2","SLC22A16","IFI6",
  "ENSG00000124593","GBP4","JDP2","CFD","IGF2BP1","RPLP1","PECAM1","PRKG1","ANO5","AMOT",
  "CRISPLD2","IGHV5-51","H4C8","SNX18","MAP3K1","ID4","MERTK","SLC15A4","MS4A4E","HOXB2",
  "AARS1","LAIR2","TRBV25-1","GYPC","NAIP","H2BC4","TSPAN12","RASAL3","MICAL2","GPR174",
  "PTPN6","NUDCD1","TRBV2","COL26A1","HLA-DMA","ACTN1","BICRAL","MFSD2B","ELOVL7","TFPI",
  "CAMK4","H2AC11","NR4A3","NME4","PAQR7","PHF19","FCER1A","CSTA","IGLV1-47","LEF1","MCFD2",
  "UNC93B1","DGKH","DDT","RETREG1","CLPB","IGHV4-31","OST4","TRIM3","CD93","SLC40A1",
  "ARMCX5","ANO4","NPIPB4","C10orf143","HECTD2","SAPCD2","NCOA3","PFDN1","IER3IP1",
  "SLC2A3","CES1","LIME1","CXCL16","SHCBP1","TRBV18","MAFIP","LYN","YBX3","ZNF627","CMIP",
  "ZNF331","TRBV7-6","RPS8","RNF11","NBPF26","GALNT16","FCN1","TRPS1","CCDC9B","PLBD1",
  "MAD1L1","PF4V1","ENSG00000280571","CDIP1","TRAF3IP3","CPNE8","SFTPD","FLT3LG","DYSF",
  "TNFSF10","SYT17","MYO7A","SHISA5","CASP3","POLM","PRPF4","GPC6","ABHD16A","DUS2",
  "IGHV1-3","CDC45","SMARCD3","CAMKK1","PCK2","HSP90B1","TMC6","PLSCR3","PLCL1","MTSS1",
  "PDLIM1","RPS13","SERPINE1","H3C4","LIMD1","KCNB2","MNDA","ZFP36L1","PRR5L","OAS2","SKAP1",
  "CRTAP","MEIS1","CARD19","CERK","CCNDBP1","BIN1","TRMT1","PITPNB","RAB9B","PTPN4","ENAH",
  "RUNX3","LEPROTL1","NR6A1","HBA2","RAB13","PTMS","ADGRE5","RNF24","ARMC12","MKI67",
  "EPB41L1","SUN2","DBNDD2","ENSG00000272442","TRBV5-5","GSTM1","DPEP2","TNS3","SCML4",
  "STK25","GZMM"
)

logfc_columns <- colnames(data)[-1]  # Exclude the gene identifier column
patient_ids_logfc <- str_extract(logfc_columns, "(?<=patient_)\\d+(?=_log2FC)")

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF", ]
# cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros

# Filter the log2FC data for genes in the pathway
# Assuming the first column is 'gene_id'; change if different
pathway_data <- data %>%
  filter(Gene %in% gene_list_500)

# Check if any genes from the gene_list are missing in the data
missing_genes <- setdiff(gene_list_500, pathway_data$Gene)
if(length(missing_genes) > 0){
  warning(paste("The following genes are missing in the log2FC data:", paste(missing_genes, collapse = ", ")))
}

# Calculate the mean log2FC for the pathway genes for each patient
# Transpose the data to have patients as rows and genes as columns
pathway_activity <- pathway_data %>%
  select(-Gene) %>%  # Remove gene identifier column
  pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
  mutate(patient_id = str_extract(patient_col, "\\d+")) %>%  # Extract numeric patient ID
  group_by(patient_id) %>%
  summarize(pathway_activity = mean(log2FC, na.rm = TRUE))  # Calculate mean log2FC

# Merge the pathway_activity with cox_input_df
cox_input_df <- cox_input_df %>%
  left_join(pathway_activity, by = "patient_id")
# cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9), ]
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9, 1, 4), ]


# Central Memory CD8 clonal expansion
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_Memory_Precursor_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Memory_Precursor_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_Memory_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_Memory_Precursor_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C1/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C1.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion$patient_id <- rownames(clonal_expansion)
merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity", "Arm", "MGMT", "IDH", "Sex", "Age")], clonal_expansion, by = "patient_id")

merged_df <- merged_df[!is.na(merged_df$pathway_activity) & !is.na(merged_df$clonal_expansion), ]

# spearman_value <- cor(merged_df$pathway_activity,
#                       merged_df$clonal_expansion,
#                       method = "spearman")
# 
# spearman_value

merged_df[merged_df$IDH == "UNK", "IDH"] = "NEG"
merged_df
merged_df <- merged_df[!is.infinite(merged_df$clonal_expansion), ]
merged_df <- merged_df[merged_df$clonal_expansion != 0, ]

#############################################################################################################
# gemini correlation metrics
#############################################################################################################
# Load the mccr package
library(mccr)
# Initialize list to store results
results <- list()

# 2. Both columns numeric
results$'Pearson (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "pearson")
results$'Spearman (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "spearman")
results$'Kendall (Numeric)' <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "kendall")

# 3. Numeric pathway_activity vs. Binary clonal_expansion (> 1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
# Check for variance in the binary column
if (length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- cor(merged_df$pathway_activity, merged_df$clonal_binary_gt1, method = "pearson")
} else {
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- NA
  warning("Binary clonal expansion column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 4. Binary pathway_activity (> median) vs. Numeric clonal_expansion
median_pathway <- median(merged_df$pathway_activity)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity > median_pathway)
# Check for variance in the binary column
if (length(unique(merged_df$pathway_binary_median)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson")
} else {
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- NA
  warning("Binary pathway activity column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 5. Both columns binary
# Check for variance in both binary columns
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Phi coefficient (equivalent to Pearson on 0/1 binary data)
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  
  # Matthews Correlation Coefficient (MCC) using the mccr package
  # Note: mccr expects factors or numeric {0, 1} or {-1, 1}. Our 0/1 format works.
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
  
} else {
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  warning("One or both binary columns lack variance. Cannot calculate Phi or MCC.")
}


# Find the metric with the most negative correlation
best_metric_name <- NULL
highest_inverse_corr <- 0 # Initialize to 0, looking for most negative

for (metric_name in names(results)) {
  value <- results[[metric_name]]
  if (!is.na(value) && value < highest_inverse_corr) {
    highest_inverse_corr <- value
    best_metric_name <- metric_name
  }
}

# Print the results
cat("Correlation Results:\n\n")
for (metric_name in names(results)) {
  cat(sprintf("%s: %.4f\n", metric_name, results[[metric_name]]))
}

cat("\n--------------------------------------------------\n")
cat("Best Inverse Correlation Found:\n")
cat(sprintf("Method: %s\n", best_metric_name))
cat(sprintf("Correlation Value: %.4f\n", highest_inverse_corr))
cat("--------------------------------------------------\n")

# Display the median used for pathway activity binarization
cat(sprintf("\nMedian used for 'pathway_activity' binarization: %.4f\n", median_pathway))





###############################################################################################################
# applying linear model
###############################################################################################################

# 1) Make sure Sex is treated as a factor (categorical) -----------------------
merged_df$Sex <- factor(merged_df$Sex)      # ‘F’ will be the reference level by default

# 2) Fit the model ------------------------------------------------------------
fit <- lm(clonal_expansion ~ pathway_activity + Sex + Age, data = merged_df)

# 3) Inspect the results ------------------------------------------------------
summary(fit)            # classic R output
#– or –#
library(broom)
tidy(fit)               # neat tibble with estimates, SEs and p-values
###############################################

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# 2. Add a grouping column ----------------------------------------------------------
merged_df$Group <- dplyr::case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                          ~ "Uncategorised"
)

# 3. Plot with ggplot2 --------------------------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = pathway_activity,
                      y = clonal_expansion,
                      colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Pathway activity vs. clonal expansion",
       subtitle = "Points coloured by survivor group",
       x = "Pathway activity",
       y = "Clonal expansion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor"  = "#4DBBD5",
    "Uncategorised"       = "grey60"
  ))

























################################################################################################################################
# starting with the Classical Monocyte Analysis
################################################################################################################################

############################################################################################################################################################################################
# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)

setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis")

# 1. Load the Seurat object
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_C_Mono_Cells_res_1_seurat_obj_exp_Pre_C1.RDS")

# Check the metadata
metadata <- seurat_obj@meta.data
head(metadata)

# Ensure that 'TimePoint' and 'Patient' columns exist
if(!all(c("TimePoint", "Patient") %in% colnames(metadata))) {
  stop("The Seurat object does not contain 'TimePoint' and/or 'Patient' columns in metadata.")
}

# 2. Prepare pseudo-bulk RNA-seq counts
# Aggregate counts by Patient and TimePoint

# Extract the raw counts matrix
# Commonly in the "RNA" assay
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Add Patient and TimePoint information to the metadata
metadata <- seurat_obj@meta.data %>%
  select(Patient, TimePoint)

# **Modified Step:** Create a unique identifier for each pseudo-bulk sample with "patient_" prefix
metadata$Sample <- paste("patient", metadata$Patient, metadata$TimePoint, sep = "_")
# Example: "patient_1_Pre", "patient_1_C1", "patient_2_Pre", etc.

# Sum the counts for cells belonging to the same Sample
# Transpose counts to have cells as rows and genes as columns for aggregation
counts_df <- as.data.frame(t(as.matrix(counts)))
counts_df$Sample <- metadata$Sample

# Aggregate counts by Sample
pseudo_bulk_counts <- counts_df %>%
  group_by(Sample) %>%
  summarise_all(sum)

# **Corrected Step:** Convert back to a matrix with genes as rows and samples as columns
# Since 'pseudo_bulk_counts' has samples as rows and genes as columns,
# we need to transpose it to have genes as rows and samples as columns.

# Remove the 'Sample' column and transpose
pseudo_bulk_counts_mat <- t(as.matrix(pseudo_bulk_counts[,-1]))

# Assign gene names as rownames
rownames(pseudo_bulk_counts_mat) <- rownames(counts)  # Correct: Gene names

# Assign sample names as column names
colnames(pseudo_bulk_counts_mat) <- pseudo_bulk_counts$Sample

# Verify dimensions
dim(pseudo_bulk_counts_mat)  # Should be genes x samples (e.g., 19343 x 24)

write.csv(pseudo_bulk_counts_mat, file = "exp_C_Mono_Cells_pseudo_bulk_RNASeq.csv")

# 3. Prepare the design matrix for edgeR
# Extract Patient and TimePoint information for each Sample
sample_info <- data.frame(
  Sample = colnames(pseudo_bulk_counts_mat)
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_")

# **Modified Step:** Reconstruct the Patient identifier with the "patient_" prefix
sample_info$Patient <- paste(sample_info$Patient_Prefix, sample_info$Patient_ID, sep = "_")
# Example: "patient_1", "patient_2", etc.

# Ensure factors are correctly set
sample_info$Patient <- factor(sample_info$Patient)
sample_info$TimePoint <- factor(sample_info$TimePoint, levels = c("Pre", "C1"))

# 4. Create DGEList object
dge <- DGEList(counts = pseudo_bulk_counts_mat)

# 5. Normalize the data
dge <- calcNormFactors(dge)

# 6. Create the design matrix
# Model: ~ Patient + TimePoint
design <- model.matrix(~ Patient + TimePoint, data = sample_info)

# 7. Estimate dispersion
dge <- estimateDisp(dge, design)

# 8. Fit the model and perform differential expression
fit <- glmFit(dge, design)
# Assuming "TimePointC1" is the coefficient for C1 vs Pre
lrt <- glmLRT(fit, coef = "TimePointC1")

# 9. Extract differential expression results
de_results <- topTags(lrt, n = Inf)$table
de_results <- as.data.frame(de_results)
# Add gene names as a column
de_results$Gene <- rownames(de_results)

# 10. Save the differential expression results
write.csv(de_results, file = "exp_C_Mono_Cells_differential_expression_results.csv", row.names = FALSE)

# 11. Create input file for GSEA
# Typically, GSEA requires a ranked list of genes. One common ranking metric is logFC multiplied by the sign of the p-value.

# Here, we'll use the signed -log10(p-value) multiplied by the sign of logFC
# Alternatively, you can choose other ranking metrics based on your preference

# To handle p-values of 0, add a small pseudocount
de_results$PValue[de_results$PValue == 0] <- 1e-300

de_results$RankingMetric <- -log10(de_results$PValue) * sign(de_results$logFC)

# Order the genes by the ranking metric
ranked_genes <- de_results %>%
  arrange(desc(RankingMetric)) %>%
  select(Gene, RankingMetric)

# Save the ranked list for GSEA
# **Option 1:** GCT format
# GCT format requires specific headers. Here's how to structure it:

# Create a GCT header
gct_header <- data.frame(
  NAME = c("#1.2"),
  DESCRIPTION = c("RankingMetric"),
  stringsAsFactors = FALSE
)

# Combine header and data
gct_data <- data.frame(
  NAME = ranked_genes$Gene,
  DESCRIPTION = ranked_genes$Gene,  # GSEA often expects gene symbols or descriptions here
  RankingMetric = ranked_genes$RankingMetric,
  stringsAsFactors = FALSE
)

# Write the GCT file
write.table(gct_header, file = "exp_C_Mono_Cells_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gct_data, file = "exp_C_Mono_Cells_GSEA_ranked_genes.gct", sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE)

# **Option 2:** RNK format
# RNK format is simpler and often preferred for ranked lists

write.table(ranked_genes, file = "exp_C_Mono_Cells_GSEA_ranked_genes.rnk", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Inform the user that the process is complete
cat("Pseudo-bulk RNA-seq aggregation, differential expression analysis, and GSEA input file creation are complete.\n")
cat("Differential expression results saved to 'differential_expression_results.csv'.\n")
cat("GSEA ranked gene list saved to 'GSEA_ranked_genes.gct' and 'GSEA_ranked_genes.rnk'.\n")





######################################################################################
library(dplyr)
library(tidyr)
# 1. Read the Pseudo-Bulk Counts Matrix
# -------------------------------------

# Specify the path to your pseudo_bulk_counts_mat.csv
counts_csv_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_C_Mono_Cells_pseudo_bulk_RNASeq.csv"

# Read the CSV file
# Assuming the first column contains gene names and is set as rownames
counts_df <- read.csv(counts_csv_path, row.names = 1, check.names = FALSE)

# Verify the dimensions and a snippet of the data
cat("Dimensions of the counts matrix:", dim(counts_df), "\n")
head(counts_df)


# 4. Normalize Counts Using "Pre" as Baseline
# -------------------------------------------

# Subset the counts matrix to include only selected genes
selected_counts <- counts_df

# Parse sample names to extract Patient and TimePoint information
# Assuming sample names are in the format "patient_X_TimePoint", e.g., "patient_1_Pre"
sample_names <- colnames(selected_counts)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = sample_names,
  stringsAsFactors = FALSE
) %>%
  separate(Sample, into = c("Patient_Prefix", "Patient_ID", "TimePoint"), sep = "_", remove = FALSE) %>%
  mutate(Patient = paste(Patient_Prefix, Patient_ID, sep = "_")) %>%
  select(Sample, Patient, TimePoint)

# View sample information
print(sample_info)

# Get a list of unique patients
patients <- unique(sample_info$Patient)
cat("Number of unique patients:", length(patients), "\n")

# Initialize a list to store log2 fold changes
log2fc_list <- list()

# Loop through each patient to compute log2FC
for (patient in patients) {
  # Identify "Pre" and "C1" samples for the patient
  pre_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "Pre") %>%
    pull(Sample)
  
  c1_sample <- sample_info %>%
    filter(Patient == patient & TimePoint == "C1") %>%
    pull(Sample)
  
  # Check if both samples are present
  if (length(pre_sample) == 1 & length(c1_sample) == 1) {
    # Extract counts for "Pre" and "C1" samples
    pre_counts <- selected_counts[, pre_sample]
    c1_counts <- selected_counts[, c1_sample]
    
    # Compute log2 fold change with a pseudocount of 1 to avoid division by zero
    log2fc <- log2((c1_counts + 1) / (pre_counts + 1))
    
    # Ensure that log2fc is a named vector with gene names
    names(log2fc) <- rownames(selected_counts)
    
    # Store in the list
    log2fc_list[[patient]] <- log2fc
  } else {
    warning(paste("Missing 'Pre' or 'C1' sample for patient:", patient))
  }
}

# Combine the list into a normalized matrix
# Rows: Genes, Columns: Patients (log2FC)
normalized_mat <- do.call(cbind, log2fc_list)

# Rename columns to indicate log2FC
colnames(normalized_mat) <- paste0(colnames(normalized_mat), "_log2FC")

# Inspect the normalized matrix
dim(normalized_mat)
head(normalized_mat)

# 5. Scale the Normalized Data
# ----------------------------

# Scaling can be done per gene to standardize the log2FC across patients
# This is useful for visualization purposes

# Transpose, scale, and transpose back
scaled_mat <- t(scale(t(normalized_mat)))

# Replace any NA values resulting from scaling with 0
scaled_mat[is.na(scaled_mat)] <- 0

# Inspect the scaled matrix
dim(scaled_mat)
head(scaled_mat)

# 6. Save the Normalized and Scaled Data as CSV
# ---------------------------------------------

# Prepare the normalized data for saving
normalized_df <- as.data.frame(normalized_mat)
normalized_df$Gene <- rownames(normalized_df)

# Reorder columns to have 'Gene' as the first column
normalized_df <- normalized_df %>%
  select(Gene, everything())

# Save the normalized log2FC data
write.csv(normalized_df, file = "exp_C_Mono_Cells_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Normalized log2 fold change data saved to 'exp_C_Mono_normalized_log2FC_C1_vs_Pre.csv'.\n")

# Prepare the scaled data for saving
scaled_df <- as.data.frame(scaled_mat)
scaled_df$Gene <- rownames(scaled_df)

# Reorder columns to have 'Gene' as the first column
scaled_df <- scaled_df %>%
  select(Gene, everything())

# Save the scaled normalized data
write.csv(scaled_df, file = "exp_C_Mono_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv", row.names = FALSE)
cat("Scaled normalized data saved to 'exp_C_Mono_scaled_normalized_log2FC_C1_vs_Pre.csv'.\n")



#########################################################################################################################################################
# Plotting heatmap

# Install necessary packages if not installed
# install.packages("pheatmap")

library(pheatmap)

# 1. Read the file
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_C_Mono_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv", header = TRUE, stringsAsFactors = FALSE)

# The columns are expected to be:
# Gene, patient_10_log2FC, patient_12_log2FC, patient_13_log2FC, patient_14_log2FC,
# patient_18_log2FC, patient_19_log2FC, patient_20_log2FC, patient_21_log2FC,
# patient_2_log2FC, patient_3_log2FC, patient_5_log2FC, patient_7_log2FC

# 2. Define the gene list
gene_list <- c(
  "PCDH1", "KRT86", "CH25H", "RNASE1", "RASA4", "IGKV3D-11", "IGHV4-34", "IL1R2", "TMEM45A", "HLA-DQA1",
  "AREG", "FOXN4", "MYL4", "IGHV3-11", "THBS1", "TRANK1", "ANO9", "WFDC1", "IFIT1", "EBF4",
  "BTN3A1", "SALL3", "TRBV6-7", "ENSG00000284989", "CCR7", "IL2RG", "IFI44L", "PPARG", "EPSTI1", "ZG16B",
  "DEFA3", "NXF3", "HSH2D", "DLGAP2", "IGKV1-8", "CDK14", "TNNT1", "AGMO", "CTSV", "MCEMP1",
  "FN1", "HERC6", "DNAH17", "LTF", "H2AC13", "ADAMTS2", "TMEM89", "C1orf115", "IGLVI-70", "IGKV6-21",
  "ACSL1", "IGHV1-24", "GLYATL3", "CETP", "RFLNB", "NSMF", "RBPMS2", "ENSG00000283782", "AIM2", "IGKV3D-15",
  "IGHV3-23", "IGLV8-61", "ADAMTS9", "NOXA1", "PRKAG1", "IGHV3-64D", "FFAR3", "ITGA4", "IGKV1-9", "VASH1",
  "RAB13", "ALOX5AP", "DDIT4", "NFIL3", "OAS3", "LONRF3", "TPM2", "TCF7L2", "IRS2", "IL10",
  "DYNLL1", "TPST1", "GPR155", "HEXD", "PARP14", "ENSG00000259316", "SERPINB2", "OR2B11", "LAD1", "THBD",
  "CES1", "PLD4", "AGXT2", "ERICH3", "TTLL8", "CCDC170", "ITGAX", "TP63", "TP53I13", "RSAD2",
  "CREBZF", "HLA-DQB1", "LYVE1", "ZFYVE9", "BTN2A2", "SIGLEC1", "IGKV2-30", "SLAMF7", "LIFR", "TUBB2A",
  "H1-2", "CLEC2D", "PAQR7", "GPR160", "RUVBL1", "ACSM1", "ZNF358", "HOXA10", "FAM156B", "HTRA1",
  "NHLRC3", "ADAM19", "KLRC2", "GBP4", "TBX19", "SAMSN1", "RNASE2", "SHCBP1", "TMEM132C", "MX1",
  "CHST6", "INHBA", "CIITA", "HLA-DOA", "MAP1LC3B", "NAALADL2", "NPIPB4", "SLC8A1", "H2BC12", "HPD",
  "CYSLTR1", "IFI44", "NXPH4", "FOXD1", "EHD1", "RASAL3", "PCDHA4", "DMRTC1", "IGKV2D-24", "ENSG00000267228",
  "IGLV10-54", "FCGBP", "CHKB-CPT1B", "LTB", "CDCA8", "DDB2", "EREG", "IGKV1-17", "ABCG1", "INSM1",
  "SDF2L1", "FOS", "CLVS2", "SLC4A10", "AMPH", "AKR1C1", "CMPK2", "PRSS16", "ECHDC2", "CAMK2D",
  "DUSP15", "HBG1", "SETBP1", "LAMC1", "MYO10", "HOXA5", "PRLR", "ENSG00000285708", "HSPA5", "DGKK",
  "ENPP1", "HMSD", "OTOF", "RASGEF1C", "SLAMF8", "DUSP1", "TNFRSF17", "UBE2J1", "F12", "CLIC4",
  "JAK2", "HSD11B1", "CFAP251", "CYCS", "ACSM3", "TMEM217", "IFI6", "CYP4F3", "TENT5A", "PCDH17",
  "SPSB3", "ADGRB3", "RPRM", "FAM166B", "ENSG00000284337", "SDHAF3", "P2RY6", "NFATC2IP", "SKA1", "ZBTB46",
  "IGLV3-19", "TRAV34", "BPIFB3", "GNG8", "PER1", "CYB561", "HEG1", "SCML1", "STEAP2", "FLT1",
  "MGST1", "C4BPB", "C3AR1", "GABRG1", "ZNF281", "SOD2", "CD38", "C16orf96", "SPRED1", "BTBD19",
  "ME1", "SCYGR4", "HDGF", "FMNL3", "GABRD", "CREB3L1", "TRAV38-2DV8", "MRAP2", "CDH12", "NKX3-1",
  "SNCAIP", "PPARGC1B", "AZIN2", "NCKAP1", "CACFD1", "IFIT3", "THNSL2", "BATF3", "NRARP", "ZNF263",
  "CAMK1D", "GATD3B", "GALNT16", "NETO1", "ELL2", "FAM53B", "MYO7A", "POU2AF1", "IGKV5-2", "ADGRE3",
  "S100A12", "PARP12", "DBN1", "PTGIR", "WHRN", "IL1RAP", "HLA-F", "ROPN1L", "RFX6", "IGHV6-1",
  "IGKV2D-40", "GABARAPL2", "GUCY1A2", "ITGA1", "VASN", "CSMD2", "SELENOP", "TRBV10-3", "CREB5", "CELA1",
  "LDLRAD4", "SDC4", "TSPAN13", "PRR7", "HEY2", "CLDN9", "KIZ", "HOXD11", "C2orf72", "ENSG00000285330",
  "SPDYE2B", "CD70", "RGL3", "KRT80", "PDGFC", "ENSG00000268173", "TESPA1", "INA", "PSPH", "IGLV1-41",
  "PCED1B", "PLAT", "USP9Y", "SEPSECS", "ZWINT", "GLIS3", "DUSP19", "RETREG1", "IFRD1", "TAP1",
  "ANKRD44", "ANOS1", "BRME1", "DCHS2", "SAMD14", "ARHGEF9", "FAM9B", "DDX60", "RNF126", "ZNF627",
  "SPTBN1", "ADRA2C", "ZNF20", "CARMIL1", "PCDH15", "TUBB8B", "BCL11B", "PHGDH", "RUNX3", "TDRD9",
  "NT5DC2", "PDE7B", "HERC5", "POLR3E", "S1PR5", "PILRB", "GLIPR1L1", "SUMF1", "FAM210A", "NDUFA4",
  "MS4A2", "SPR", "FAM200B", "PIGX", "TRAV29DV5", "VANGL1", "KCNQ4", "CC2D1B", "TRERF1", "SPTSSA",
  "EBF1", "IGHG4", "AGPAT1", "KRT72", "MGST3", "CD4", "IFI27", "CDC25C", "PRRT2", "GPRIN3",
  "STX19", "COL6A3", "PRDX5", "CBLN3", "ITGAM", "IL17RD", "OBSL1", "CPEB1", "YBX3", "PARP9",
  "ST20", "PIP4P2", "SLC5A10", "SCGB3A1", "BCKDHA", "SLC16A14", "TTLL3", "MORN3", "DLX1", "ACOD1",
  "STARD3NL", "CBLB", "MT-ATP8", "STIMATE-MUSTN1", "OSR2", "EPX", "LKAAEAR1", "TMEM198", "DYRK3", "MDK",
  "TMEM176A", "FAM110D", "AGFG1", "SLC24A3", "RAB3B", "PLSCR4", "KDM2B", "ENSG00000263620", "RASSF1", "WDR49",
  "CHORDC1", "LRMDA", "COL8A1", "CRIP1", "SAMD9L", "ZNF724", "CATSPERG", "NR2F1", "ZNF781", "VEGFD",
  "HIVEP3", "SNX30", "RPS6KA2", "PNPLA8", "CABIN1", "CD6", "MAP1B", "CHST3", "CYTL1", "SPATA20",
  "CATSPER3", "ENSG00000250803", "PGGHG", "TGFBR3L", "MTFR1", "SLC4A8", "IGLL1", "SRC", "DEFA4", "GPR137C",
  "ENSG00000285938", "IGHV4-61", "FAIM", "PIGL", "EIF3CL", "H3C2", "TJP1", "TBC1D10A", "ID2", "YWHAE",
  "CD244", "PRR5-ARHGAP8", "CLEC10A", "ENSG00000284292", "MYH15", "RIPPLY3", "RCVRN", "UTRN", "FLT3", "SOX7",
  "KCNE5", "BASP1", "CROCC2", "PDX1", "PTGFRN", "FKBP1B", "ORAI3", "KLHL4", "FZD8", "ABCB1",
  "DENND2A", "ATP10A", "ZNF311", "ACMSD", "PCMT1", "NOD1", "FGF11", "NEURL1B", "TAS2R4", "TLR9",
  "NRK", "CA6", "DHPS", "ALDH1A2", "TCF7", "SETDB2", "TMEM91", "PCCA", "NNAT", "STMP1",
  "STMN3", "TBX21", "TMC8", "EFHD1", "YWHAG", "AFF3", "IL12RB1", "ENOSF1", "IL26", "CD5",
  "UVSSA", "ODF3B", "UVRAG", "NPIPB12", "CPLANE1", "RNF207", "SYCP2", "PRICKLE3", "ARF4", "XBP1"
)

# 3. Subset the data for only the requested genes
subset_data <- data[data$Gene %in% gene_list,]

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Define patient order
patient_order <- c(14, 7, 5, 20, 21, 19, 2, 3, 13)

# Create column names based on the patient order
ordered_columns <- paste0("patient_", patient_order, "_log2FC")

# Check that all ordered columns exist in the data
# If any column doesn't exist, you might need to verify the column naming in the CSV.
all(ordered_columns %in% colnames(subset_data))

# columns that exist
ordered_columns <- intersect(ordered_columns, colnames(subset_data))

# Reorder the data columns to match the desired patient order
subset_data <- subset_data[, c("Gene", ordered_columns)]

# Set rownames to Gene, remove Gene column for heatmap matrix
rownames(subset_data) <- subset_data$Gene
subset_data <- subset_data[,-1]  # remove the Gene column

# Convert to a numeric matrix if needed
mat <- as.matrix(subset_data)

patient_ids <- patient_order                              # same order you plotted
group_tag   <- ifelse(patient_ids %in% short_term_survivor_group,
                      "short term", "long term")
new_labels  <- paste0("patient_", patient_ids, " (", group_tag, ")")

colnames(mat) <- new_labels        # replace the matrix column names

# 4. Plot the heatmap
# You can adjust clustering, scaling, color palettes, etc. as needed.
# Run pheatmap and store the result in an object
res <- pheatmap(mat, 
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                scale = "row",
                show_rownames = TRUE, 
                show_colnames = TRUE,
                main = "Heatmap of Selected Genes",
                fontsize_row = 6,
                fontsize_col = 8)

# Extract the hierarchical clustering tree for rows
hc_rows <- res$tree_row

# Choose the number of clusters (e.g., k=5). This number can be adjusted based on 
# the structure you see in the dendrogram.
k <- 20

# Cut the dendrogram into k clusters
clusters <- cutree(hc_rows, k = k)

# Create a data frame with gene and cluster assignments
gene_clusters <- data.frame(
  Gene = rownames(mat),
  Cluster = clusters
)

# Write the cluster assignments to a CSV file so you can inspect them outside R
write.csv(gene_clusters, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_C_Mono_Cells_heatmap_gene_clusters.csv", row.names = FALSE)





########################################################################################################################################################################
# trying to see if there is any correlation between innate and adaptive immunity
########################################################################################################################################################################
library(tidyverse)
# Read survival data
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# Read the log2FC data
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_C_Mono_Cells_scaled_normalized_log2FC_C1_vs_Pre.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

# # Read the log2FC data
# data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_CCR2_Pos_Cells_normalized_log2FC_C1_vs_Pre.csv",
#                  header = TRUE,
#                  stringsAsFactors = FALSE)

################################################################################
# 1. Get the 500 most-significant genes (by FDR) from the DE results
################################################################################
de_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_C_Mono_Cells_differential_expression_results.csv"

de_results <- read.csv(de_file, header = TRUE, stringsAsFactors = FALSE)

sig500 <- de_results %>% 
  arrange(FDR) %>%                   # lowest FDR first
  slice_head(n = 500)                # take top 500 rows

gene_list_500   <- sig500$Gene
gene_list_up    <- sig500 %>% filter(logFC  > 0) %>% pull(Gene)
gene_list_down  <- sig500 %>% filter(logFC  < 0) %>% pull(Gene)

################################################################################
# 2. Helper to collapse a gene set to a per-patient activity
################################################################################
get_activity <- function(expr_df, genes_vec, new_name){
  expr_df %>%
    filter(Gene %in% genes_vec) %>%            # keep only the gene set
    select(-Gene) %>%                          # drop gene column
    pivot_longer(cols = everything(),
                 names_to  = "patient_col",
                 values_to = "log2FC") %>%
    mutate(patient_id = str_extract(patient_col, "\\d+")) %>%
    group_by(patient_id) %>%
    summarize("{new_name}" := mean(log2FC, na.rm = TRUE),
              .groups = "drop")
}

pathway_all  <- get_activity(data, gene_list_500,  "pathway_activity_all")
pathway_up   <- get_activity(data, gene_list_up,   "pathway_activity_up")
pathway_down <- get_activity(data, gene_list_down, "pathway_activity_down")

################################################################################
# 3. Add the three pathway-activity variables to the Cox input
################################################################################

# Extract patient IDs from the cox_input_df 'patient' column
# The patient IDs are the last element after splitting by "-"
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()
cox_input_df <- cox_input_df[cox_input_df$site == "UF", ]
# cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Create a new column 'patient_id' by extracting and cleaning patient IDs
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros


cox_input_df <- cox_input_df %>%
  left_join(pathway_all,  by = "patient_id") %>%
  left_join(pathway_up,   by = "patient_id") %>%
  left_join(pathway_down, by = "patient_id")

# cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9), ]
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(8, 9, 1, 4), ]


# Central Memory CD8 clonal expansion
clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_Memory_Precursor_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C1/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C1.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_All_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
clonal_expansion$patient_id <- rownames(clonal_expansion)
merged_df <- merge(cox_input_df[, c("patient_id", "pathway_activity_all", "pathway_activity_up", "pathway_activity_down", "Arm", "MGMT", "IDH", "Sex", "Age")], clonal_expansion, by = "patient_id")

merged_df <- merged_df[!is.na(merged_df$pathway_activity_all) & !is.na(merged_df$clonal_expansion), ]

# spearman_value <- cor(merged_df$pathway_activity,
#                       merged_df$clonal_expansion,
#                       method = "spearman")
# 
# spearman_value

merged_df[merged_df$IDH == "UNK", "IDH"] = "NEG"
merged_df
merged_df <- merged_df[!is.infinite(merged_df$clonal_expansion), ]
merged_df <- merged_df[merged_df$clonal_expansion != 0, ]

#############################################################################################################
# gemini correlation metrics
#############################################################################################################
# Load the mccr package
library(mccr)
# Initialize list to store results
results <- list()

# 2. Both columns numeric
results$'Pearson (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "pearson")
results$'Spearman (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "spearman")
results$'Kendall (Numeric)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_expansion, method = "kendall")

# 3. Numeric pathway_activity vs. Binary clonal_expansion (> 1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
# Check for variance in the binary column
if (length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- cor(merged_df$pathway_activity_all, merged_df$clonal_binary_gt1, method = "pearson")
} else {
  results$'Point-Biserial (Numeric Pathway vs. Binary Clonal > 1)' <- NA
  warning("Binary clonal expansion column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 4. Binary pathway_activity (> median) vs. Numeric clonal_expansion
median_pathway <- median(merged_df$pathway_activity_all)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity_all > median_pathway)
# Check for variance in the binary column
if (length(unique(merged_df$pathway_binary_median)) > 1) {
  # Point-Biserial is equivalent to Pearson between numeric and 0/1 binary
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson")
} else {
  results$'Point-Biserial (Binary Pathway > Median vs. Numeric Clonal)' <- NA
  warning("Binary pathway activity column has no variance (all values are the same). Cannot calculate Point-Biserial correlation.")
}


# 5. Both columns binary
# Check for variance in both binary columns
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  # Phi coefficient (equivalent to Pearson on 0/1 binary data)
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  
  # Matthews Correlation Coefficient (MCC) using the mccr package
  # Note: mccr expects factors or numeric {0, 1} or {-1, 1}. Our 0/1 format works.
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
  
} else {
  results$'Phi Coefficient (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  results$'MCC (Binary Pathway > Median vs. Binary Clonal > 1)' <- NA
  warning("One or both binary columns lack variance. Cannot calculate Phi or MCC.")
}


# Find the metric with the most negative correlation
best_metric_name <- NULL
highest_inverse_corr <- 0 # Initialize to 0, looking for most negative

for (metric_name in names(results)) {
  value <- results[[metric_name]]
  if (!is.na(value) && value < highest_inverse_corr) {
    highest_inverse_corr <- value
    best_metric_name <- metric_name
  }
}

# Print the results
cat("Correlation Results:\n\n")
for (metric_name in names(results)) {
  cat(sprintf("%s: %.4f\n", metric_name, results[[metric_name]]))
}

cat("\n--------------------------------------------------\n")
cat("Best Inverse Correlation Found:\n")
cat(sprintf("Method: %s\n", best_metric_name))
cat(sprintf("Correlation Value: %.4f\n", highest_inverse_corr))
cat("--------------------------------------------------\n")

# Display the median used for pathway activity binarization
cat(sprintf("\nMedian used for 'pathway_activity' binarization: %.4f\n", median_pathway))



###############################################################################################################
# applying linear model
###############################################################################################################

# 1) Make sure Sex is treated as a factor (categorical) -----------------------
merged_df$Sex <- factor(merged_df$Sex)      # ‘F’ will be the reference level by default

# 2) Fit the model ------------------------------------------------------------
fit <- lm(clonal_expansion ~ pathway_activity_all + Sex + Age, data = merged_df)

# 3) Inspect the results ------------------------------------------------------
summary(fit)            # classic R output
#– or –#
library(broom)
tidy(fit)               # neat tibble with estimates, SEs and p-values
###############################################

# 1. Define the two survivor groups -------------------------------------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# 2. Add a grouping column ----------------------------------------------------------
merged_df$Group <- dplyr::case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                          ~ "Uncategorised"
)

# 3. Plot with ggplot2 --------------------------------------------------------------
library(ggplot2)

ggplot(merged_df, aes(x = pathway_activity_all,
                      y = clonal_expansion,
                      colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Pathway activity vs. clonal expansion",
       subtitle = "Points coloured by survivor group",
       x = "Pathway activity",
       y = "Clonal expansion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor"  = "#4DBBD5",
    "Uncategorised"       = "grey60"
  ))


























########################################################################################################################################################################
# Script: innate_vs_adaptive_correlation.R
# Purpose: Correlate innate (NC Monocytes/DC) pathway activity with adaptive (Central Memory CD8 clonal expansion)
#          using two independent datasets / DE contrasts.  For each pair of (data file, DE file) the top‑500 genes
#          (lowest FDR) are taken, pathway activities are computed per patient, then averaged across the two datasets.
#          Only patients with a pathway activity from **both** datasets are retained.
########################################################################################################################################################################

# ────────────────────────────────────────────────────────────────────────────────
# 0.  Load libraries
# ────────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(survival)
library(data.table)
library(ggfortify)
library(Rcpp)
library(cowplot)
library(stringr)
library(Seurat)
library(mccr)

# ────────────────────────────────────────────────────────────────────────────────
# 1.  Define input files  ▸▸▸  **EDIT THESE PATHS AS NEEDED**
# ────────────────────────────────────────────────────────────────────────────────
# ‑ log2FC matrices (gene × patient) -------------------------------------------
log2fc_file_1 <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/NC_Monocytes_control_exp_scaled_normalized_log2FC_C1_vs_Pre.csv"
log2fc_file_2 <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_scaled_normalized_log2FC_C1_vs_Pre.csv"  # <‑‑ example

# ‑ DE result tables (must contain columns: Gene, FDR) -------------------------
de_file_1     <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/differential_expression_results.csv"
de_file_2     <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/exp_DC_differential_expression_results.csv"        # <‑‑ example

# ‑ Survival metadata -----------------------------------------------------------
survival_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv"

# ‑ Clonal expansion results ----------------------------------------------------
clonal_expansion <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_Memory_Precursor_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt"
# clonal_expansion <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Effector_and_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt"
# clonal_expansion <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/C1_vs_C2.txt"
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C1/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C1.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/Pre_vs_C2/clonal_expansion_simpson_Central_Memory_CD8_T_cells_division/S_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# clonal_expansion <- read.table("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/C1_vs_C2/clonal_expansion_simpson_All_CD8_T_cells_division/C1_vs_C2.txt", header = TRUE, stringsAsFactors = FALSE)
# ────────────────────────────────────────────────────────────────────────────────
# 2.  Helper: given a log2FC matrix + DE table, compute per‑patient pathway score
# ────────────────────────────────────────────────────────────────────────────────
get_pathway_activity <- function(log2fc_path, de_path, top_n = 500) {
  # Read data
  data_df      <- fread(log2fc_path) %>% as_tibble()
  de_df        <- fread(de_path)      %>% as_tibble()
  
  # Top‑N genes by FDR (assumes column names 'Gene' and 'FDR')
  top_genes    <- de_df %>% arrange(FDR) %>% slice_head(n = top_n) %>% pull(Gene)
  
  # Filter expression matrix to those genes
  pathway_mat  <- data_df %>% filter(Gene %in% top_genes)
  
  # Check for missing genes
  missing_genes <- setdiff(top_genes, pathway_mat$Gene)
  if (length(missing_genes) > 0) {
    warning(sprintf("%d / %d top genes were not present in %s",
                    length(missing_genes), top_n, basename(log2fc_path)))
  }
  
  # Reshape to long ⇒ compute mean per patient
  pathway_activity <- pathway_mat %>%
    select(-Gene) %>%
    pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
    mutate(patient_id = str_extract(patient_col, "\\d+")) %>%
    group_by(patient_id) %>%
    summarise(pathway_activity = mean(log2FC, na.rm = TRUE), .groups = "drop")
  
  return(pathway_activity)
}

# ────────────────────────────────────────────────────────────────────────────────
# 3.  Compute pathway activities for BOTH datasets
# ────────────────────────────────────────────────────────────────────────────────
pathway_activity_1 <- get_pathway_activity(log2fc_file_1, de_file_1)
pathway_activity_2 <- get_pathway_activity(log2fc_file_2, de_file_2)

# Merge and keep only patients with scores from BOTH datasets
combined_pathway <- full_join(pathway_activity_1, pathway_activity_2,
                              by = "patient_id", suffix = c("_1", "_2")) %>%
  filter(!is.na(pathway_activity_1) & !is.na(pathway_activity_2)) %>%
  mutate(average_pathway_activity = (pathway_activity_1 + pathway_activity_2) / 2)

# ────────────────────────────────────────────────────────────────────────────────
# 4.  Survival / covariate data
# ────────────────────────────────────────────────────────────────────────────────
survival_data <- read_csv(survival_file)
cox_input_df  <- survival_data %>%
  select(-sample_id, -timepoint, -integrated_sample_type) %>% distinct() %>%
  filter(site == "UF") %>%
  mutate(patient_id = str_remove(str_split(Patient, "-", simplify = TRUE)[,4], "^0+"))

# Merge pathway scores
cox_input_df <- cox_input_df %>% left_join(combined_pathway, by = "patient_id")

# Remove patients lacking averaged activity (already handled) + manual exclusions
cox_input_df <- cox_input_df[!cox_input_df$patient_id %in% c(1,4,8,9), ]

# ────────────────────────────────────────────────────────────────────────────────
# 5.  Clonal expansion data & merge
# ────────────────────────────────────────────────────────────────────────────────
clonal_expansion <- fread(clonal_expansion) %>% as_tibble(rownames = "patient_id")

print(cox_input_df)

merged_df <- cox_input_df %>%
  select(patient_id, average_pathway_activity, Arm, MGMT, IDH, Sex, Age) %>%
  inner_join(clonal_expansion, by = "patient_id") %>%
  rename(pathway_activity = average_pathway_activity)

print(merged_df)
  
merged_df <- merged_df %>%
  filter(!is.na(pathway_activity) & !is.na(clonal_expansion)) %>%
  filter(!is.infinite(clonal_expansion) & clonal_expansion != 0)

# Normalise unknown IDH
merged_df$IDH[merged_df$IDH == "UNK"] <- "NEG"

print(merged_df)

# ────────────────────────────────────────────────────────────────────────────────
# 6.  Correlation metrics (numeric & binary combos)
# ────────────────────────────────────────────────────────────────────────────────
results <- list()
results$"Pearson (Numeric)"  <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "pearson")
results$"Spearman (Numeric)" <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "spearman")
results$"Kendall (Numeric)"  <- cor(merged_df$pathway_activity, merged_df$clonal_expansion, method = "kendall")

# Numeric vs binary (>1)
merged_df$clonal_binary_gt1 <- as.numeric(merged_df$clonal_expansion > 1)
results$"Point‑Biserial (Numeric Pathway vs Binary Clonal >1)" <-
  if (length(unique(merged_df$clonal_binary_gt1)) > 1)
    cor(merged_df$pathway_activity, merged_df$clonal_binary_gt1, method = "pearson") else NA

# Binary pathway (> median) vs numeric clonal
median_pathway <- median(merged_df$pathway_activity)
merged_df$pathway_binary_median <- as.numeric(merged_df$pathway_activity > median_pathway)
results$"Point‑Biserial (Binary Pathway > median vs Numeric Clonal)" <-
  if (length(unique(merged_df$pathway_binary_median)) > 1)
    cor(merged_df$pathway_binary_median, merged_df$clonal_expansion, method = "pearson") else NA

# Both binary
if (length(unique(merged_df$pathway_binary_median)) > 1 && length(unique(merged_df$clonal_binary_gt1)) > 1) {
  results$"Phi (Binary vs Binary)" <- cor(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1, method = "pearson")
  results$"MCC (Binary vs Binary)" <- mccr(merged_df$pathway_binary_median, merged_df$clonal_binary_gt1)
} else {
  results$"Phi (Binary vs Binary)" <- NA
  results$"MCC (Binary vs Binary)" <- NA
}

# Pretty print results
cat("\nCorrelation Results\n====================\n")
walk(names(results), \(nm) cat(sprintf("%s: %.4f\n", nm, results[[nm]])))

best_metric <- names(results)[which.min(results)]
cat("\nBest inverse correlation (most negative):", best_metric,
    sprintf("\nValue: %.4f\n", results[[best_metric]]))
cat(sprintf("\nMedian used for pathway binarisation: %.4f\n", median_pathway))

# ────────────────────────────────────────────────────────────────────────────────
# 7.  Linear model: clonal_expansion ~ pathway_activity + Sex + Age
# ────────────────────────────────────────────────────────────────────────────────
merged_df$Sex <- factor(merged_df$Sex)
fit <- lm(clonal_expansion ~ pathway_activity + Sex + Age, data = merged_df)
print(summary(fit))

# ────────────────────────────────────────────────────────────────────────────────
# 8.  Plot (coloured by survivor group)
# ────────────────────────────────────────────────────────────────────────────────
short_term_survivor_group <- c(7,10,12,14,18)
long_term_survivor_group  <- c(2,3,5,13,19,20,21)
merged_df$Group <- case_when(
  merged_df$patient_id %in% short_term_survivor_group ~ "Short-term survivor",
  merged_df$patient_id %in% long_term_survivor_group  ~ "Long-term survivor",
  TRUE                                                ~ "Uncategorised")

ggplot(merged_df, aes(x = pathway_activity, y = clonal_expansion, colour = Group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Averaged pathway activity vs clonal expansion",
       subtitle = "Patients coloured by survivor group",
       x = "Average pathway activity (two datasets)",
       y = "Clonal expansion (Simpson index)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c(
    "Long-term survivor" = "#E64B35",
    "Short-term survivor" = "#4DBBD5",
    "Uncategorised"      = "grey60"))
