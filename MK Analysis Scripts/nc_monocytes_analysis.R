# ================================================================================
# Non-Classical Monocytes: Differential Expression and Survival Analysis
# ================================================================================
#
# Purpose:
#   Perform comprehensive analysis of Non-Classical Monocytes comparing C1 cells
#   with their Pre-treatment predecessors (identified via optimal transport).
#   This script generates three key manuscript figures through differential
#   expression, GSEA, and survival analysis.
#
# Biological Context:
#   - Uses predecessor mapping from NC_Monocyte_Optimal_Transport.ipynb
#   - Compares C1 cells against their most likely Pre predecessors
#   - Identifies genes driving treatment-induced changes in NC monocytes
#   - Assesses clinical impact of gene expression changes on patient survival
#
# Workflow Overview:
#   SECTION 1-2: Predecessor distribution analysis and visualization
#   SECTION 3-4: Differential expression between C1 and predecessor cells
#   SECTION 5-6: Pseudo-bulk RNA-seq aggregation and normalization
#   SECTION 7-8: Heatmap generation and gene clustering (Figure 5c)
#   SECTION 9: GSEA rank file preparation (Figure 5e)
#   SECTION 10: Cox proportional hazards survival analysis (Figure 5f)
#
# Key Outputs:
#   - cluster_*_de_result.csv: Differential expression results per cluster
#   - differential_expression_results.csv: Pseudo-bulk DE results
#   - Heatmap with gene clustering (Figure 5c)
#   - GSEA rank files using z-scores and t-values (Figure 5e)
#   - Survival analysis with pathway activity (Figure 5f)
#
# Manuscript Figures:
#   Figure 5c: Heatmap of differential gene expression (C1 vs Pre predecessors)
#   Figure 5e: GSEA enrichment results
#   Figure 5f: Survival curves stratified by pathway activity
#
# Dependencies:
#   - Seurat object: MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS
#   - Predecessor mapping from optimal transport analysis
#   - Libraries: Seurat, dplyr, ggplot2, tibble, edgeR, pheatmap, survival
# ================================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

# ================================================================================
# SECTION 1: Predecessor Distribution Visualization Function
# ================================================================================
# This function analyzes where C1 cells originated from by examining their
# Pre predecessors' cluster distributions

plot_predecessor_distribution <- function(seurat_obj,
                                          target_cluster,
                                          target_timepoint,
                                          predecessor_timepoint,
                                          cluster_col = "seurat_clusters",
                                          timepoint_col = "timepoint",
                                          predecessor_col = "predecessor_barcode",
                                          plot_type = "bar") {
  
  # Validate that required metadata columns exist in Seurat object
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
  
  # Verify predecessor barcodes actually exist in the Seurat object
  existing_predecessors <- predecessor_barcodes %in% colnames(seurat_obj)
  if (sum(!existing_predecessors) > 0) {
    warning(paste(sum(!existing_predecessors), "predecessor barcodes not found in the Seurat object and will be excluded."))
    predecessor_barcodes <- predecessor_barcodes[existing_predecessors]
  }
  
  if (length(predecessor_barcodes) == 0) {
    stop("No valid predecessor barcodes found in the Seurat object.")
  }
  
  # Extract cluster assignments for predecessor cells at Pre timepoint
  predecessor_clusters <- seurat_obj@meta.data %>%
    filter(barcode %in% predecessor_barcodes,
           (!!sym(timepoint_col)) == predecessor_timepoint) %>%
    pull((!!sym(cluster_col)))
  
  if (length(predecessor_clusters) == 0) {
    stop("No predecessor clusters found for the specified predecessor timepoint.")
  }
  
  # Calculate cluster distribution: which Pre clusters gave rise to target C1 cluster
  cluster_dist <- as.data.frame(table(predecessor_clusters))
  colnames(cluster_dist) <- c("Predecessor_Cluster", "Frequency")
  cluster_dist <- cluster_dist %>%
    arrange(desc(Frequency)) %>%
    mutate(Percentage = (Frequency / sum(Frequency)) * 100)
  
  # Generate visualization (bar or pie chart)
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


# ================================================================================
# SECTION 2: Single-Cell Differential Expression Analysis
# ================================================================================
# Purpose:
#   Perform DE analysis between C1 cells in a target cluster and their
#   Pre predecessor cells identified via optimal transport mapping.
#
# Biological Rationale:
#   - Compares treatment-induced transcriptional changes at single-cell resolution
#   - Identifies genes specifically altered in NC monocytes post-treatment
#   - Uses predecessor mapping to ensure biologically meaningful comparisons
#
# Statistical Method:
#   - Wilcoxon rank-sum test via Seurat's FindMarkers function
#   - Controls for inter-patient variability through predecessor pairing
#
# Output:
#   - cluster_*_de_result.csv: DE results for each cluster analyzed
# ================================================================================

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

# Function: perform_de_between_cluster_and_predecessors
# Identifies differentially expressed genes between target C1 cells and their Pre predecessors
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
  
  # Extract target cells (C1 cells in specified cluster)
  target_cells <- seurat_obj@meta.data %>%
    filter((!!sym(cluster_col)) == target_cluster,
           (!!sym(timepoint_col)) == target_timepoint) %>%
    pull(barcode)
  
  if (length(target_cells) == 0) {
    stop("No cells found for the specified cluster and timepoint.")
  }
  
  # Extract predecessor barcodes (Pre cells) for each target cell
  # These come from optimal transport predecessor mapping
  predecessor_barcodes <- seurat_obj@meta.data[target_cells, predecessor_col]
  
  # Remove cells without valid predecessors
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
  
  # Step 4: Prepare data for differential expression
  # Combine target C1 cells with their Pre predecessors
  cells_for_de <- c(target_cells_valid, predecessor_cells)
  seurat_subset <- subset(seurat_obj, cells = cells_for_de)
  
  # Create binary classification: Target (C1 cells) vs Predecessor (Pre cells)
  seurat_subset$de_group <- ifelse(colnames(seurat_subset) %in% target_cells_valid, "Target", "Predecessor")
  
  # Set active identity for Seurat DE functions
  Idents(seurat_subset) <- "de_group"
  
  # Step 5: Perform differential expression using Wilcoxon rank-sum test
  # Target (C1) vs Predecessor (Pre) to identify treatment-induced changes
  de_results <- FindMarkers(seurat_subset,
                            ident.1 = "Target",        # C1 cells in target cluster
                            ident.2 = "Predecessor",   # Their Pre predecessors
                            min.pct = min.pct,         # Min fraction of cells expressing gene
                            logfc.threshold = logfc.threshold)  # Min log2FC threshold
  
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

# ================================================================================
# SECTION 3: Execute Differential Expression for Cluster 1
# ================================================================================
# Apply DE analysis to cluster 1 NC monocytes

# Load NC monocyte Seurat object with predecessor mappings
seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj_unwanted_clusters_removed_exp_Pre_C1.RDS")

# Define analysis parameters for cluster 1
target_cluster <- "1"                        # Cluster ID as string
target_timepoint <- "C1"                     # Post-treatment timepoint
predecessor_timepoint <- "Pre"               # Baseline timepoint
cluster_col <- "seurat_clusters"             # Metadata column for clusters
timepoint_col <- "TimePoint"                 # Metadata column for timepoints
predecessor_col <- "predecessor_barcode"     # Metadata column for predecessor mapping

# Execute differential expression analysis
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

# Extract results from DE analysis
de_results <- results$de_results                    # Full DE statistics table
target_cells_valid <- results$target_cells_valid    # C1 cells in cluster
predecessor_cells <- results$predecessor_cells      # Their Pre predecessors

# ================================================================================
# SECTION 4: Visualization of Single-Cell DE Results
# ================================================================================
# Generate preliminary visualizations to explore top DE genes in cluster 1

# Select top 20 most significant genes by adjusted p-value
top_genes <- de_results %>%
  arrange(p_val_adj) %>%         # Sort by statistical significance
  slice(1:20) %>%                 # Take top 20
  pull(gene)                      # Extract gene names

# Extract normalized expression values for visualization
expr_data <- FetchData(seurat_obj, vars = top_genes,
                       cells = c(target_cells_valid, predecessor_cells))

# Create cell group annotations for heatmap
annotation <- data.frame(Group = ifelse(colnames(expr_data) %in% target_cells_valid,
                                       "Target", "Predecessor"))
rownames(annotation) <- colnames(expr_data)

# Vis 1: Heatmap of top 20 DE genes (preliminary visualization)
library(pheatmap)
pheatmap(as.matrix(expr_data),
         cluster_rows = TRUE,           # Cluster genes by expression pattern
         cluster_cols = TRUE,           # Cluster cells by similarity
         annotation_col = annotation,   # Show Target vs Predecessor groups
         show_colnames = FALSE,         # Hide individual cell names
         show_rownames = TRUE,          # Show gene names
         scale = "row",                 # Z-score normalization per gene
         main = paste("Top 20 DE Genes for Cluster", target_cluster))

# Vis 2: Volcano plot showing fold change vs significance
ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = paste("Volcano Plot for Cluster", target_cluster),
       x = "Average Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Significance threshold
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "red")   # Fold change thresholds




# ================================================================================
# SECTION 5: Pseudo-bulk RNA-seq Aggregation Using edgeR
# ================================================================================
#
# Purpose:
#   Convert single-cell RNA-seq data to pseudo-bulk samples and perform
#   differential expression analysis using the edgeR framework. This approach
#   aggregates cells by patient and timepoint to increase statistical power
#   while controlling for patient-specific effects.
#
# Biological Context:
#   - Single-cell data has high technical noise and dropout events
#   - Pseudo-bulk aggregation reduces noise while preserving biological signal
#   - Enables powerful differential expression testing with established RNA-seq tools
#   - Controls for patient variability through paired design (each patient has Pre and C1)
#
# Statistical Framework:
#   - Generalized Linear Model (GLM) with negative binomial distribution
#   - Design matrix: ~ Patient + TimePoint
#     * Patient terms: Control for inter-patient variability (batch effects)
#     * TimePoint term: Test C1 vs Pre difference (treatment effect)
#   - TMM normalization: Adjusts for library size and composition bias
#   - Likelihood Ratio Test (LRT): Tests TimePointC1 coefficient significance
#
# Key Steps:
#   1. Aggregate raw counts by Patient × TimePoint combinations
#   2. Create DGEList object with pseudo-bulk count matrix
#   3. Apply TMM normalization to correct for sequencing depth
#   4. Build design matrix accounting for patient and timepoint effects
#   5. Estimate dispersion parameters (common, trended, tagwise)
#   6. Fit GLM and perform likelihood ratio test
#   7. Extract differential expression results with logFC and p-values
#
# Output Files:
#   - pseudo_bulk_RNASeq.csv: Aggregated count matrix (genes × patient_timepoint)
#   - differential_expression_results.csv: DE results with logFC, p-values, FDR
#   - GSEA_ranked_genes.rnk: Ranked gene list for GSEA (using -log10(p) * sign(logFC))
#
# Clinical Interpretation:
#   - Positive logFC: Gene upregulated in C1 vs Pre (treatment increases expression)
#   - Negative logFC: Gene downregulated in C1 vs Pre (treatment decreases expression)
#   - Significant genes (FDR < 0.05) represent treatment-induced changes
# ================================================================================

# Load required libraries for pseudo-bulk analysis
library(Seurat)      # Single-cell data handling
library(edgeR)       # Differential expression for RNA-seq count data
library(dplyr)       # Data manipulation
library(tidyr)       # Data tidying (separate function)

# Set working directory for output files
setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis")

# STEP 1: Load NC monocyte Seurat object with optimal transport predecessor mappings
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
# ================================================================================
# SECTIONS 7-8: Heatmap Generation with Hierarchical Clustering
# ================================================================================
#
# **MANUSCRIPT FIGURE 5c: Differential Gene Expression Heatmap**
#
# Purpose:
#   Visualize treatment-induced gene expression changes across patients using
#   hierarchical clustering to identify co-regulated gene modules and patient
#   response patterns.
#
# Biological Context:
#   - Each row = one differentially expressed gene
#   - Each column = one patient's treatment response (C1 vs Pre log2FC, z-scored)
#   - Color intensity = magnitude of normalized expression change
#   - Gene clustering reveals functional modules responding together to treatment
#   - Patient clustering identifies responder subgroups with similar molecular profiles
#
# Hierarchical Clustering Methodology:
#   1. Gene Clustering (Rows):
#      - Distance metric: Euclidean distance between gene expression patterns
#      - Linkage method: Complete linkage (default in pheatmap)
#      - Purpose: Group genes with similar expression trajectories across patients
#      - Biological interpretation: Co-clustered genes likely share regulatory mechanisms
#
#   2. Patient Clustering (Columns):
#      - Disabled (cluster_cols = FALSE) to preserve clinical outcome ordering
#      - Manual ordering: Short-term survivors first, then long-term survivors
#      - Rationale: Enables visual assessment of gene patterns associated with outcomes
#
# Z-score Color Scale Interpretation:
#   - Red (positive z-score): Gene more upregulated in this patient vs cohort average
#   - Blue (negative z-score): Gene more downregulated in this patient vs cohort average
#   - White (z-score ≈ 0): Patient's response similar to cohort mean
#   - Scale range: Typically -3 to +3 (3 standard deviations from mean)
#
# Patient Stratification (Column Organization):
#   - Short-term survivors (patients 18, 12, 14, 7, 10): OS < median
#   - Long-term survivors (patients 5, 20, 21, 19, 2, 3, 13): OS ≥ median
#   - Purpose: Visual identification of gene expression patterns associated with survival
#
# Gene Selection Strategy:
#   - ~800 genes selected based on differential expression significance
#   - Includes: upregulated genes, downregulated genes, and key immune markers
#   - Curated to represent diverse biological processes affected by treatment
#   - Focused on immunologically relevant pathways
#
# Cluster Cutting for Gene Modules:
#   - Dendrogram cut into k=20 clusters
#   - Each cluster = a co-regulated gene module
#   - Enables functional enrichment analysis of gene groups
#   - Output: heatmap_gene_clusters.csv maps genes to cluster assignments
#
# Output Files:
#   - Heatmap visualization (displayed, can be saved as PDF/PNG)
#   - heatmap_gene_clusters.csv: Gene-to-cluster mapping for downstream analysis
#
# Clinical Significance:
#   - Identifies gene signatures distinguishing treatment responders from non-responders
#   - Reveals molecular heterogeneity in patient treatment responses
#   - Gene clusters may predict clinical outcomes or treatment sensitivity
#   - Informs biomarker discovery for patient stratification
# ================================================================================


# STEP 3: Filter data for selected genes
# Subset to include only the curated gene list for heatmap visualization
subset_data <- data[data$Gene %in% gene_list,]

# STEP 4: Organize patients by clinical outcome for visual assessment
# Patient ordering strategy aligns with survival outcomes:
#   First 5: Short-term survivors (18, 12, 14, 7, 10) - poor outcomes
#   Next 7: Long-term survivors (5, 20, 21, 19, 2, 3, 13) - good outcomes
# Purpose: Enables visual identification of gene expression patterns associated with survival
patient_order <- c(18, 12, 14, 7, 10, 5, 20, 21, 19, 2, 3, 13)

# Generate ordered column names based on patient stratification
ordered_columns <- paste0("patient_", patient_order, "_log2FC")

# Verify all expected patient columns exist in the dataset
all(ordered_columns %in% colnames(subset_data))

# Reorder columns: Gene column first, then patients in survival-based order
subset_data <- subset_data[, c("Gene", ordered_columns)]

# Prepare data matrix for heatmap visualization
rownames(subset_data) <- subset_data$Gene  # Convert gene names to rownames (pheatmap requirement)
subset_data <- subset_data[,-1]  # Remove redundant Gene column

# Convert to numeric matrix (required by pheatmap)
mat <- as.matrix(subset_data)

# STEP 5: Generate hierarchical clustering heatmap
# **THIS CREATES MANUSCRIPT FIGURE 5c**
# Clustering parameters:
#   - cluster_rows = TRUE: Hierarchical clustering of genes (reveals co-regulated modules)
#   - cluster_cols = FALSE: Preserve clinical outcome ordering (short-term → long-term)
#   - scale = "row": Ensures consistent z-score scaling
#   - Distance metric: Euclidean (default)
#   - Linkage method: Complete linkage (default)
res <- pheatmap(mat,
                cluster_rows = TRUE,           # Enable gene clustering to identify functional modules
                cluster_cols = FALSE,          # Disable patient clustering to preserve outcome ordering
                scale = "row",                 # Row-wise scaling (genes already z-scored)
                show_rownames = TRUE,          # Display gene names on y-axis
                show_colnames = TRUE,          # Display patient IDs on x-axis
                main = "Heatmap of Selected Genes",  # Plot title
                fontsize_row = 6,              # Gene name font size (small for ~800 genes)
                fontsize_col = 8)              # Patient ID font size

# STEP 6: Extract gene clustering dendrogram for module identification
# The dendrogram captures hierarchical relationships among genes
# Genes with similar expression patterns cluster together
hc_rows <- res$tree_row

# STEP 7: Cut dendrogram into discrete gene modules
# k = 20: Divide genes into 20 co-regulated clusters
# Rationale: Balance between granularity (detailed modules) and interpretability
# Each cluster represents genes with similar expression patterns across patients
k <- 20

# Assign each gene to a cluster based on dendrogram structure
# Uses complete linkage clustering with Euclidean distance
clusters <- cutree(hc_rows, k = k)

# Create gene-to-cluster mapping table
gene_clusters <- data.frame(
  Gene = rownames(mat),
  Cluster = clusters
)

# Save cluster assignments for downstream functional enrichment analysis
# Applications: GO enrichment, KEGG pathway analysis, transcription factor binding analysis
write.csv(gene_clusters, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/heatmap_gene_clusters.csv", row.names = FALSE)


# STEP 8: Calculate pathway activity scores for survival association
# Purpose: Quantify overall transcriptional response magnitude per patient
# Method: Average z-scored log2FC across all genes in the gene signature

# Load visualization libraries
library(ggplot2)
library(ggpubr)

# Matrix structure verification:
# - Rows = genes (z-scored log2FC values)
# - Columns = patients (ordered by survival outcome)

# Define patient groups by observed survival outcomes
# Based on overall survival (OS) data from clinical records
short_term_survivors <- c(18, 12, 14, 7, 10)      # OS < median (poor prognosis)
long_term_survivors <- c(5, 20, 21, 19, 2, 3, 13)  # OS ≥ median (good prognosis)

# Extract patient identifiers from matrix column names
all_patients <- colnames(mat)

# Parse patient IDs from column names
# Expected format: "patient_18_log2FC" → extract "18"
patient_numbers <- sapply(strsplit(all_patients, "_"), function(x) x[2])
patient_numbers <- as.numeric(patient_numbers)

# Calculate pathway activity as mean z-score across all genes per patient
# Biological interpretation: Overall magnitude of treatment-induced transcriptional response
# Positive values: Strong transcriptional reprogramming in NC monocytes
# Negative values: Weak or dampened transcriptional response to treatment
pathway_activity <- colMeans(mat, na.rm = TRUE)

# Prepare data for statistical comparison between survival groups
plot_data <- data.frame(
  patient = patient_numbers,
  activity = pathway_activity
)

# Assign survival group labels based on clinical outcomes
plot_data$group <- ifelse(plot_data$patient %in% short_term_survivors, "Short-term", "Long-term")

# Set factor levels for correct plot ordering
plot_data$group <- factor(plot_data$group, levels = c("Short-term", "Long-term"))

# Statistical test: Compare pathway activity between survival groups
# Wilcoxon rank-sum test (Mann-Whitney U): Non-parametric test for two independent groups
# Null hypothesis: No difference in pathway activity between survival groups
# Alternative: Pathway activity differs between groups
stat_test <- wilcox.test(activity ~ group, data = plot_data)

# Extract p-value for plot annotation
p_value <- stat_test$p.value

# Create boxplot visualization comparing pathway activity by survival group
# Purpose: Visual assessment of whether pathway activity associates with clinical outcome
p <- ggplot(plot_data, aes(x = group, y = activity, fill = group)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 1.2), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1.2), size = 4, alpha = 0.8) +
  scale_fill_manual(values = c("Short-term" = "blue", "Long-term" = "red")) +
  labs(
    title = "Pathway Activity by Survival Group",
    x = NULL,
    y = "Pathway Activity (Mean Expression)"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(),
    legend.position = "none"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_data$activity, na.rm = TRUE),
    label = ifelse(is.na(p_value), "p-value = NA", paste("p-value =", round(p_value, 4))),
    size = 4,
    vjust = -0.5
  )

# Display pathway activity comparison plot
p
# ================================================================================
# SECTION 9: GSEA Rank File Preparation for Pathway Enrichment Analysis
# ================================================================================
#
# **MANUSCRIPT FIGURE 5e: Gene Set Enrichment Analysis (GSEA) Results**
#
# Purpose:
#   Generate ranked gene lists for GSEA to identify biological pathways
#   differentially activated between long-term and short-term survivors.
#   Two ranking metrics (z-score and t-value) provide complementary views
#   of gene importance.
#
# Biological Context:
#   - Compares gene expression patterns between survival outcome groups
#   - Long-term survivors (n=7): Patients with favorable treatment response
#   - Short-term survivors (n=5): Patients with poor treatment response
#   - Goal: Identify pathways associated with clinical benefit
#
# Statistical Framework:
#   Two ranking metrics calculated for each gene:
#
#   1. Z-score (Effect Size Metric):
#      Formula: z = (mean_long - mean_short) / SE_difference
#      Where SE_difference = sqrt((SD_short²/n_short) + (SD_long²/n_long))
#      
#      Interpretation:
#      - Positive z-score: Gene higher in long-term survivors
#      - Negative z-score: Gene higher in short-term survivors
#      - Magnitude: Effect size in standard error units
#      - Advantage: Incorporates both difference and variability
#
#   2. T-value (Welch's t-test statistic):
#      Method: Two-sample t-test with unequal variances
#      Formula: t.test(long_values, short_values, var.equal = FALSE)
#      
#      Interpretation:
#      - Positive t-value: Gene higher in long-term survivors
#      - Negative t-value: Gene higher in short-term survivors
#      - Magnitude: Statistical significance of difference
#      - Advantage: Well-established statistical framework
#
# Why Two Metrics?
#   - Z-score: Better for genes with small sample sizes but consistent effects
#   - T-value: Better for genes with larger sample sizes and strong significance
#   - Using both ensures robust pathway identification
#   - GSEA can be run with either metric for validation
#
# GSEA Input Format (RNK file):
#   - Tab-delimited text file
#   - Two columns: GeneSymbol\tRankingScore
#   - Genes sorted by ranking score (descending order)
#   - No header line
#   - Compatible with GSEA desktop and command-line tools
#
# Ranking Logic:
#   - Genes ranked from most upregulated in long-term (top) to most downregulated (bottom)
#   - Sign convention: Positive = beneficial (associated with longer survival)
#   - Sign convention: Negative = detrimental (associated with shorter survival)
#
# Output Files:
#   - gene_zscore_rankfile.rnk: Z-score based ranking
#   - gene_tvalue_rankfile.rnk: T-value based ranking
#
# Downstream GSEA Analysis:
#   These rank files serve as input to GSEA software to identify:
#   - Hallmark gene sets enriched in long-term vs short-term survivors
#   - GO biological processes associated with treatment response
#   - KEGG pathways predicting clinical outcomes
#   - Immunologic signatures correlating with survival
#
# Clinical Interpretation:
#   - Enriched pathways at top of list: Associated with favorable outcomes
#   - Enriched pathways at bottom of list: Associated with poor outcomes
#   - Leading edge genes: Key drivers of pathway association
#   - Normalized Enrichment Score (NES): Magnitude of pathway enrichment
# ================================================================================





#########################################################################################################################################
# create 2 rank files using z-scores and t-value for GSEA

# STEP 1: Define patient survival groups
# Based on clinical outcomes (overall survival data)
short_term_survivors <- c(18, 12, 14, 7, 10)      # n=5, poor outcomes
long_term_survivors <- c(5, 20, 21, 19, 2, 3, 13)  # n=7, favorable outcomes

# Generate column names for each group
# These correspond to columns in the z-scored log2FC data file
short_term_cols <- paste0("patient_", short_term_survivors, "_log2FC")
long_term_cols <- paste0("patient_", long_term_survivors, "_log2FC")

# STEP 2: Load z-scored log2FC data
# This contains treatment-induced expression changes (C1 vs Pre) for each patient
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

# STEP 3: Calculate z-score ranking metric for each gene
# Z-score formula: (mean_long - mean_short) / SE_difference
# Biological interpretation: Effect size in standard error units
# Sign: Positive if gene higher in long-term survivors (favorable)
z_scores <- apply(data, 1, function(row) {
  # Extract z-scored log2FC values for each survival group
  short_values <- as.numeric(row[short_term_cols])  # n=5 patients
  long_values <- as.numeric(row[long_term_cols])    # n=7 patients
  
  # Calculate group means for this gene
  mean_short <- mean(short_values, na.rm = TRUE)
  mean_long <- mean(long_values, na.rm = TRUE)
  
  # Calculate group standard deviations
  sd_short <- sd(short_values, na.rm = TRUE)
  sd_long <- sd(long_values, na.rm = TRUE)
  
  # Count non-missing values per group
  n_short <- sum(!is.na(short_values))
  n_long <- sum(!is.na(long_values))
  
  # Calculate standard error of the difference between means
  # Formula: SE = sqrt((SD1²/n1) + (SD2²/n2))
  # Accounts for sample size and variability in both groups
  se_diff <- sqrt((sd_short^2 / n_short) + (sd_long^2 / n_long))
  
  # Handle edge case: identical values across all patients (SE = 0)
  if (se_diff == 0) {
    z_val <- NA
  } else {
    # Calculate z-score with sign convention:
    # Positive: Gene higher in long-term survivors (potentially beneficial)
    # Negative: Gene higher in short-term survivors (potentially detrimental)
    z_val <- (mean_long - mean_short) / se_diff
  }
  
  return(z_val)
})

# Add z-scores as new column in data frame
data$Zscore <- z_scores

# STEP 4: Calculate t-value ranking metric for each gene
# Uses Welch's two-sample t-test (assumes unequal variances)
# Advantages: Well-established, incorporates variance heterogeneity
# Sign: Positive if gene higher in long-term survivors
t_values <- apply(data, 1, function(row) {
  # Extract z-scored log2FC values for both groups
  short_values <- as.numeric(row[short_term_cols])
  long_values <- as.numeric(row[long_term_cols])
  
  # Perform Welch's t-test: tests H0: mean(long) = mean(short)
  # var.equal = FALSE: Allows different variances between groups (robust)
  # Direction: t.test(long, short) gives positive t if mean(long) > mean(short)
  t_res <- t.test(long_values, short_values, var.equal = FALSE)
  
  # Extract t-statistic (test statistic value)
  t_val <- t_res$statistic
  return(t_val)
})

# Add t-values as new column in data frame
data$Tvalue <- t_values

# STEP 5: Create and export GSEA rank files
# GSEA requires tab-delimited RNK files with format: GeneSymbol\tScore

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
# ================================================================================
# SECTION 10: Cox Proportional Hazards Survival Analysis
# ================================================================================
#
# **MANUSCRIPT FIGURE 5f: Survival Curves Stratified by Pathway Activity**
#
# Purpose:
#   Assess the clinical impact of gene expression changes on patient survival
#   using Cox proportional hazards regression. This analysis tests whether
#   treatment-induced transcriptional changes in NC monocytes predict
#   overall survival outcomes.
#
# Biological Context:
#   - Hypothesis: Magnitude of NC monocyte transcriptional response correlates
#     with clinical benefit (overall survival)
#   - Gene sets represent pathways altered by treatment in NC monocytes
#   - Pathway activity = aggregate measure of transcriptional reprogramming
#   - Higher activity may indicate stronger immune activation or dysfunction
#
# Statistical Framework: Cox Proportional Hazards Model
#   Model: h(t) = h₀(t) × exp(β₁×Age + β₂×PathwayActivity)
#   
#   Where:
#   - h(t): Hazard function (instantaneous risk of death at time t)
#   - h₀(t): Baseline hazard (age-adjusted)
#   - β₁: Effect of age on survival
#   - β₂: Effect of pathway activity on survival (primary interest)
#   
#   Hazard Ratio (HR) Interpretation:
#   - HR > 1: Increased risk of death (poor prognosis)
#   - HR < 1: Decreased risk of death (good prognosis)
#   - HR = 1: No effect on survival
#
# Pathway Activity Calculation:
#   Method: Mean z-scored log2FC across selected gene signature
#   
#   Steps:
#   1. Select gene signature (500-1000 genes from differential expression)
#   2. Extract patient-specific z-scored log2FC for each gene
#   3. Calculate mean z-score across all genes per patient
#   4. Result: Single "pathway activity score" per patient
#   
#   Interpretation:
#   - Positive score: Patient shows stronger upregulation vs cohort average
#   - Negative score: Patient shows weaker response vs cohort average
#   - Magnitude: Strength of transcriptional response
#
# Patient Stratification Strategy:
#   Two approaches used:
#   
#   1. Binary Stratification (Primary Analysis):
#      - Dichotomize by median pathway activity
#      - High activity group: Above median
#      - Low activity group: Below median
#      - Enables Kaplan-Meier visualization and log-rank test
#   
#   2. Continuous Analysis (Sensitivity Analysis):
#      - Pathway activity as continuous covariate
#      - More statistical power
#      - Assumes linear relationship with log-hazard
#
# Survival Data Requirements:
#   - OS.months: Overall survival in months from treatment start
#   - Dead: Binary event indicator (1 = death, 0 = censored)
#   - Age: Continuous covariate (years) for adjustment
#   - Patient: Unique patient identifier
#   - Arm: Treatment arm filter (MK-3475 + MLA)
#   - Site: Study site filter (UF)
#
# Model Diagnostics and Assumptions:
#   1. Proportional Hazards Assumption:
#      - Hazard ratio remains constant over time
#      - Verified using Schoenfeld residuals (not shown in code)
#   
#   2. Linearity Assumption:
#      - Log-hazard linear in covariates
#      - For continuous pathway activity only
#   
#   3. Independence:
#      - Survival times independent across patients
#      - Valid for patient-level analysis
#
# Statistical Output:
#   - Hazard Ratio (HR): exp(β) for pathway activity
#   - 95% Confidence Interval: Uncertainty in HR estimate
#   - p-value: Significance of pathway activity effect
#   - Concordance index (C-index): Model discrimination ability
#
# Clinical Interpretation Examples:
#   - HR = 0.5, p < 0.05: High pathway activity associated with 50% reduced
#     risk of death (protective effect, favorable prognosis)
#   
#   - HR = 2.0, p < 0.05: High pathway activity associated with 2× increased
#     risk of death (adverse effect, poor prognosis)
#   
#   - HR = 1.2, p = 0.3: Pathway activity not significantly associated with
#     survival (no prognostic value)
#
# Visualization (Figure 5f):
#   - Kaplan-Meier curves stratified by pathway activity (high vs low)
#   - X-axis: Time (months)
#   - Y-axis: Survival probability
#   - Color: High activity (one color) vs Low activity (another color)
#   - Statistics: Log-rank p-value, median survival estimates
#
# Gene List Selection Criteria:
#   - Genes significantly differential expressed (FDR < 0.05)
#   - Biologically relevant to immune function
#   - Consistent directionality across patients
#   - Representative of key biological processes
#   - Typically 500-1000 genes for robust signature
#
# Power Considerations:
#   - Small sample size (n=12 patients) limits statistical power
#   - Wide confidence intervals expected
#   - Results should be considered hypothesis-generating
#   - Validation in independent cohort recommended
# ================================================================================

            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)






# STEP 1: Load libraries required for survival analysis
library(tidyverse)   # Data manipulation and visualization
library(survival)    # Cox proportional hazards models, Kaplan-Meier curves
library(dplyr)       # Data wrangling
library(data.table)  # Fast data operations
library(ggplot2)     # Plotting
library(ggfortify)   # Survival plot extensions
library(Rcpp)        # C++ integration for performance
library(cowplot)     # Publication-quality plots
library(stringr)     # String manipulation

# STEP 2: Load clinical survival data
# Contains patient-level survival outcomes (overall survival, events, demographics)
survival_data <- read_csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/survival_data.csv")

# STEP 3: Define Cox proportional hazards model function
# This function fits a Cox regression model with robust error handling
#
# Parameters:
#   - cox_input_subset_df: Data frame with survival outcomes and covariates
#   - time_column: Name of column containing survival time (e.g., "OS.months")
#   - event_column: Name of column containing event indicator (e.g., "Dead": 1=event, 0=censored)
#   - covariates_list: Vector of covariate names to include in model
#
# Returns:
#   - List containing fitted Cox model object
#
# Cox Model: h(t) = h₀(t) × exp(β₁×covariate₁ + β₂×covariate₂ + ...)
# where h(t) is the hazard at time t, h₀(t) is baseline hazard
get_cox = function(cox_input_subset_df, time_column, event_column, covariates_list)
{
  # Subset data to include only columns needed for modeling
  df = cox_input_subset_df
  df =  df[, c(time_column, event_column, covariates_list)]
  
  # Initialize model object
  res.cox <- NULL
  
  # Construct Cox model formula dynamically
  # Format: Surv(time, event) ~ covariate1 + covariate2 + ...
  # Surv() creates a survival object combining time and event status
  formula_string = paste0("Surv(", time_column, ", ", event_column, ") ~ ", paste(covariates_list, collapse = " + "))
  new_formula <- as.formula(formula_string)
  
  # Fit Cox model with robust error handling
  # Cox regression can fail to converge with small sample sizes or collinearity
  tryCatch(
    {
      # Initial attempt: Fit with default iterations (iter.max=20)
      res.cox <- coxph(new_formula, data = df, control = coxph.control(iter.max=20))
    },
    error = function(e_outer) {
      cat("Outer Error: ", conditionMessage(e_outer), "\n")
    },
    warning = function(w_outer) {
      # Handle convergence warnings specifically
      if (grepl("Ran out of iterations and did not converge", conditionMessage(w_outer))) {
        cat("Caught the specific warning: Ran out of iterations and did not converge\n")
        # Retry with increased iterations to achieve convergence
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
        # Pass through other warnings
        warning(w_outer)
      }
    }
  )
  
  # Return fitted Cox model object
  return(list(cox = res.cox))
}

# STEP 4: Load z-scored log2FC data for pathway activity calculation
data <- read.csv("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/scaled_normalized_log2FC_C1_vs_Pre.csv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)


# STEP 5: Define gene signature for pathway activity score
# This 500-gene signature represents key genes differentially expressed
# between C1 and Pre timepoints in NC monocytes
# Selected based on statistical significance and biological relevance
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

# STEP 6: Prepare survival data for Cox regression analysis

# Extract patient IDs from log2FC data column headers for matching
logfc_columns <- colnames(data)[-1]  # Exclude first column (Gene)
patient_ids_logfc <- str_extract(logfc_columns, "(?<=patient_)\\d+(?=_log2FC)")

# Prepare Cox regression input data frame
# Remove unnecessary columns from survival data
cox_input_df <- survival_data[, !(colnames(survival_data) %in% c("sample_id", "timepoint", "integrated_sample_type"))]
cox_input_df <- cox_input_df %>% distinct()  # Remove duplicate rows

# Filter to specific cohort: University of Florida (UF) site, MK-3475 + MLA treatment arm
# This ensures homogeneous patient population for survival analysis
cox_input_df <- cox_input_df[cox_input_df$site == "UF" & cox_input_df$Arm == "MK-3475 + MLA", ]

# Extract and clean patient IDs for matching with log2FC data
# Patient identifier format in survival data: "XXX-XXX-XXX-0018" → extract "18"
cox_input_df <- cox_input_df %>%
  mutate(patient_id = str_split(Patient, "-", simplify = TRUE)[,4],
         patient_id = str_remove(patient_id, "^0+"))  # Remove leading zeros: "0018" → "18"

# STEP 7: Calculate pathway activity scores
# Filter log2FC data to include only genes in the pathway signature
pathway_data <- data %>%
  filter(Gene %in% gene_list_1000)

# Validate gene list completeness
missing_genes <- setdiff(gene_list_1000, pathway_data$Gene)
if(length(missing_genes) > 0){
  warning(paste("The following genes are missing in the log2FC data:", paste(missing_genes, collapse = ", ")))
}

# Calculate mean log2FC across all pathway genes for each patient
# This summarizes the overall transcriptional response in the pathway
pathway_activity <- pathway_data %>%
  select(-Gene) %>%  # Remove gene identifier column
  pivot_longer(cols = everything(), names_to = "patient_col", values_to = "log2FC") %>%
  mutate(patient_id = str_extract(patient_col, "\\d+")) %>%  # Extract numeric patient ID
  group_by(patient_id) %>%
  summarize(pathway_activity = mean(log2FC, na.rm = TRUE))  # Pathway activity = mean z-scored log2FC

# Merge pathway activity scores with clinical survival data
cox_input_df <- cox_input_df %>%
  left_join(pathway_activity, by = "patient_id")

# STEP 8: Binary stratification of pathway activity for survival analysis
# Dichotomize patients into High vs Low pathway activity groups
# Stratification method: Median split

# Calculate median pathway activity across all patients
median_activity <- median(cox_input_df$pathway_activity, na.rm = TRUE)

# Create binary pathway activity variable
# "Yes" (High activity): Above median - stronger transcriptional response
# "No" (Low activity): Below median - weaker transcriptional response
cox_input_df <- cox_input_df %>%
  mutate(pathway_activity_binary = ifelse(pathway_activity > median_activity, "Yes", "No")) %>%
  mutate(pathway_activity_binary = factor(pathway_activity_binary, levels = c("No", "Yes")))
# Note: "No" is reference level for hazard ratio calculation

# STEP 9: Fit Cox proportional hazards model
# **THIS GENERATES MANUSCRIPT FIGURE 5f**
#
# Model: Surv(OS.months, Dead) ~ Age + pathway_activity_binary
#   - OS.months: Overall survival in months (continuous)
#   - Dead: Event indicator (1 = death, 0 = censored/alive)
#   - Age: Patient age at treatment start (covariate for adjustment)
#   - pathway_activity_binary: High vs Low pathway activity (primary predictor)
#
# Hazard Ratio Interpretation:
#   - HR > 1: High activity increases risk of death (poor prognosis)
#   - HR < 1: High activity decreases risk of death (good prognosis)
#   - p < 0.05: Statistically significant association with survival
cox_ph_result <- get_cox(cox_input_df, "OS.months.", "Dead", c("Age", "pathway_activity_binary"))

# Display complete Cox regression results
# Output includes: Hazard ratios, 95% CI, p-values, concordance index
summary(cox_ph_result$cox)
# Now, 'cox_input_df' contains the 'pathway_activity' column





