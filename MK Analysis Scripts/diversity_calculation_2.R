################################################################################
# CLONAL EXPANSION ANALYSIS - STEP 2: DIVERSITY CALCULATION
################################################################################
#
# PURPOSE:
#   This is the SECOND step (2/4) in the clonal expansion analysis pipeline.
#   It calculates clonal diversity indices (Shannon and Simpson) for each sample
#   to quantify the degree of clonal expansion in T cell populations.
#
# ANALYSIS PIPELINE:
#   Step 1: clonal_expansion_analysis_1.R - Prepare clonotype tables
#   Step 2 (THIS SCRIPT): Calculate diversity indices
#   Step 3: clonal_expansion_calculation_3.R - Calculate expansion metrics
#   Step 4: clonal_expansion_plots_4.R - Generate visualizations (Figure 6d)
#
# BIOLOGICAL INTERPRETATION:
#   LOWER diversity = MORE clonal expansion = HIGHER antigen specificity
#   
#   When T cells recognize specific antigens (e.g., tumor antigens), those
#   particular clones undergo massive proliferation, resulting in:
#   - Few dominant clones (low diversity)
#   - Oligoclonal expansion (focused immune response)
#   - High antigen specificity
#   
#   Conversely, high diversity indicates:
#   - Many equally-sized clones
#   - Polyclonal response (less focused)
#   - Lower antigen specificity
#
# DIVERSITY INDICES:
#   1. Shannon Index (H): Measures both richness (number of clones) and evenness
#      - Formula: H = -Σ(pi × ln(pi))
#      - Range: 0 to ln(N), where N is number of clones
#      - Higher H = more diverse (many equal-sized clones)
#      - Lower H = less diverse (few dominant clones)
#
#   2. Simpson Index (SDI): Emphasizes dominance of abundant clones
#      - Formula: SDI = 1 - Σ(pi²)  [Gini-Simpson variant]
#      - Range: 0 to 1
#      - Higher SDI = more diverse
#      - Lower SDI = less diverse (more expansion)
#
# INPUT FILES:
#   - clonotype_df_proportion_*.csv (from Step 1)
#     Proportional clone frequencies per sample
#
# OUTPUT FILES (per cell type):
#   - shannon_clonal_diversity_*.txt: Shannon index per sample
#   - simpson_clonal_diversity_*.txt: Simpson index per sample
#
################################################################################

################################################################################
# FUNCTION: simpson_diversity_index
################################################################################
# Calculate Simpson's Diversity Index (Gini-Simpson)
#
# PARAMETERS:
#   clonal_proportions - Numeric vector of clone proportions (must sum to 1)
#
# RETURNS:
#   SDI value (0 to 1, higher = more diverse)
#
# FORMULA:
#   SDI = 1 - Σ(pi²)
#   where pi is the proportion of clone i
#
################################################################################
simpson_diversity_index = function(clonal_proportions) {
  # Calculate Simpson's Diversity Index (Gini-Simpson variant)
  # SDI = 1 - Σ(ni / N)²
  # Higher values indicate more diversity
  SDI = 1 - sum(clonal_proportions^2)
  return(SDI)  
}

################################################################################
# FUNCTION: shannon_diversity_index
################################################################################
# Calculate Shannon Diversity Index (Shannon Entropy)
#
# PARAMETERS:
#   clonal_proportions - Numeric vector of clone proportions (must sum to 1)
#
# RETURNS:
#   H value (0 to ln(N), higher = more diverse)
#
# FORMULA:
#   H = -Σ(pi × ln(pi))
#   where pi is the proportion of clone i
#
################################################################################
shannon_diversity_index = function(clonal_proportions) {
  # Remove zeros to avoid log(0) errors
  clonal_proportions = clonal_proportions[clonal_proportions > 0]
  
  # Calculate Shannon entropy
  # H = -Σ(pi × ln(pi))
  # Higher values indicate more diversity
  h = (clonal_proportions * log(clonal_proportions))
  h = -sum(h)
  return(h)  
}

################################################################################
# FUNCTION: get_clonotype_diversity
################################################################################
# Calculate diversity index for all samples in a clonotype table
#
# PARAMETERS:
#   trackClonotypes_df_relative - Dataframe with clone proportions (first col = CDR3, rest = samples)
#   celltype - Name of the T cell subpopulation
#   diversity_index - Type of index: "shannon" or "simpson"
#   outfolder - Output directory for results
#
# RETURNS:
#   Dataframe with diversity index for each sample
#
################################################################################
get_clonotype_diversity = function(trackClonotypes_df_relative, celltype, diversity_index = "shannon", outfolder=".") {
  # Dictionary mapping index names to functions
  diversity_dict = list(shannon = shannon_diversity_index, simpson = simpson_diversity_index)
  
  # Calculate diversity for each sample (columns 2 onwards)
  # Column 1 is CDR3 sequence, so we skip it
  clonotype_diversity = sapply(trackClonotypes_df_relative[, 2:ncol(trackClonotypes_df_relative)], 
                               diversity_dict[[diversity_index]])
  
  # Convert to dataframe
  clonotype_diversity = as.data.frame(clonotype_diversity)
  
  # Save results to file
  write.table(clonotype_diversity, 
             paste0(outfolder, "/", diversity_index, "_clonal_diversity_", celltype, "_T_cells.txt"), 
             sep="\t", row.names=TRUE, quote=FALSE)
  
  return(clonotype_diversity)
}

################################################################################
# FUNCTION: run
################################################################################
# Main function to calculate both Shannon and Simpson indices for a cell type
#
# PARAMETERS:
#   celltype - Name of the T cell subpopulation
#
# RETURNS:
#   List containing shannon and simpson diversity dataframes
#
################################################################################
run = function(celltype) {
  # Set working directory to output location
  setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature")
  
  # Load proportional clonotype table from Step 1
  file_path = paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/clonotype_df_proportion_", celltype, "_T_cells.csv")
  trackClonotypes_df_relative <- read.csv(file_path, check.names = FALSE)
  
  # Remove CDR3.aa column (keep only sample columns with proportions)
  trackClonotypes_df_relative <- trackClonotypes_df_relative[, !(colnames(trackClonotypes_df_relative) %in% c("CDR3.aa"))]
  
  # Calculate Shannon diversity index
  shannon = get_clonotype_diversity(trackClonotypes_df_relative = trackClonotypes_df_relative, 
                                    celltype, 
                                    diversity_index = "shannon")
  
  # Calculate Simpson diversity index
  simpson = get_clonotype_diversity(trackClonotypes_df_relative = trackClonotypes_df_relative, 
                                    celltype, 
                                    diversity_index = "simpson")
  
  return(list(shannon = shannon, simpson = simpson))
}

################################################################################
# SECTION: PROCESS ALL T CELL SUBPOPULATIONS
################################################################################
# Define the same T cell subpopulation mappings as in Step 1
# This ensures consistency across the analysis pipeline

mapping <- list(
  "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
  "Activated_CD4" = c(0),
  "Effector_CD8" = c(1),
  "Effector_Memory_Precursor_CD8" = c(2),
  "Exhausted_T" = c(3),
  "Gamma_Delta_T" = c(4),
  "Active_CD4" = c(5),
  "Naive_CD4" = c(6, 9, 18),
  "Memory_CD4" = c(7),
  "Stem_Like_CD8" = c(8),
  "Effector_Memory_CD8" = c(10),
  "Central_Memory_CD8" = c(12),
  "Effector_and_Central_Memory_CD8" = c(1, 12),
  "Effector_Memory_Precursor_and_Central_Memory_CD8" = c(2, 10, 12),
  "Memory_Precursor_and_Central_Memory_CD8" = c(2, 12),
  "Effector_Memory_and_Central_Memory_CD8" = c(10, 12),
  "GZMK_Effector_Memory_CD8" = c(13),
  "Proliferating_Effector" = c(14, 16, 17),
  "All_CD8" = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
)

# Convert to dataframe
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# Get unique cell types
celltypes <- unique(celltype_to_cluster$celltype)

# Calculate diversity indices for all cell types
for (celltype in celltypes) {
  out = run(celltype)
  print(out[[1]])  # Shannon diversity results
  print(out[[2]])  # Simpson diversity results
}

################################################################################
# END OF SCRIPT
################################################################################
# 
# OUTPUT SUMMARY:
# For each T cell subpopulation, two diversity index files are created:
#   1. shannon_clonal_diversity_*.txt - Shannon diversity per sample
#   2. simpson_clonal_diversity_*.txt - Simpson diversity per sample
#
# INTERPRETATION GUIDE:
#   - LOW diversity = HIGH expansion = STRONG antigen specificity
#   - HIGH diversity = LOW expansion = WEAK/BROAD antigen response
#
# These diversity metrics serve as input for:
#   - Step 3: clonal_expansion_calculation_3.R (expansion metrics)
#   - Step 4: clonal_expansion_plots_4.R (Figure 6d visualization)
#   - Correlation analyses with clinical outcomes
#
################################################################################
