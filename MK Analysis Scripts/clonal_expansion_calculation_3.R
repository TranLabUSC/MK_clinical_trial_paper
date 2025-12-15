################################################################################
# CLONAL EXPANSION ANALYSIS - STEP 3: EXPANSION METRIC CALCULATION
################################################################################
#
# PURPOSE:
#   This is the THIRD step (3/4) in the clonal expansion analysis pipeline.
#   It calculates clonal expansion metrics for each patient by comparing diversity
#   indices between two timepoints (typically C1 vs C2).
#
# ANALYSIS PIPELINE:
#   Step 1: clonal_expansion_analysis_1.R - Prepare clonotype tables
#   Step 2: diversity_calculation_2.R - Calculate diversity indices
#   Step 3 (THIS SCRIPT): Calculate expansion metrics (diversity ratios)
#   Step 4: clonal_expansion_plots_4.R - Generate visualizations (Figure 6d)
#
# BIOLOGICAL RATIONALE:
#   Clonal expansion is measured by changes in diversity over time:
#   
#   - DECREASING diversity (ratio < 1) = EXPANSION occurred
#     → Specific clones proliferated, reducing overall diversity
#     → Indicates antigen-specific immune response
#   
#   - INCREASING diversity (ratio > 1) = NO expansion / diversification
#     → Clones became more evenly distributed
#     → Indicates polyclonal response
#   
#   - Ratio = Diversity(later timepoint) / Diversity(earlier timepoint)
#   - Lower ratio = More expansion = Stronger antigen-specific response
#
# CALCULATION METHODS:
#   1. Division (ratio): diversity_t2 / diversity_t1
#      - Values < 1: expansion (diversity decreased)
#      - Values > 1: diversification (diversity increased)
#      - Preferred method for proportional changes
#
#   2. Subtraction (difference): diversity_t2 - diversity_t1
#      - Negative values: expansion
#      - Positive values: diversification
#      - Useful for absolute change magnitude
#
# INPUT FILES:
#   - shannon_clonal_diversity_*.txt (from Step 2)
#   - simpson_clonal_diversity_*.txt (from Step 2)
#   - Sample metadata (from get_metadata.R)
#
# OUTPUT FILES (per cell type × diversity type):
#   - clonal_expansion_{shannon|simpson}_{celltype}_T_cells/
#     └── {timepoint1}_vs_{timepoint2}.txt (expansion metrics per patient)
#
################################################################################

################################################################################
# FUNCTION: get_clonal_expansion_timepoint_comparison
################################################################################
# Calculate clonal expansion by comparing diversity between two timepoints
#
# PARAMETERS:
#   compared_pair - Vector of two timepoint names to compare [earlier, later]
#                   Example: c("C1", "C2") or c("Pre", "C1")
#   metadata - Dataframe with Patient, timepoint, and sample_id columns
#   clonal_diversity_df - Dataframe with diversity values (rows=samples, col=diversity)
#   outdir - Output directory name
#   method - Calculation method: "division" (ratio) or "subtraction" (difference)
#
# RETURNS:
#   NULL (results saved to file)
#
# OUTPUT:
#   Text file with clonal expansion value per patient
#
################################################################################
get_clonal_expansion_timepoint_comparison = function(compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                      metadata,
                                                      clonal_diversity_df, 
                                                      outdir = "clonal_expansion_shannon_all_T_cells",
                                                      method = "division") {
  # Create output directory with method suffix
  outdir = paste0(outdir, "_", method) 
  dir.create(outdir)
  
  # Get list of all patients
  patients = unique(metadata$Patient)
  
  ################################################################################
  # NESTED FUNCTION: get_clonal_expansion_per_patient
  ################################################################################
  # Calculate expansion for a single patient
  get_clonal_expansion_per_patient = function(patient) {
    # Get sample ID for first timepoint (earlier)
    sample1 = metadata %>% 
      filter(Patient == patient, timepoint == compared_pair[[1]]) %>% 
      select(sample_id)
    sample1 = sample1$sample_id

    # Get sample ID for second timepoint (later)
    sample2 = metadata %>% 
      filter(Patient == patient, timepoint == compared_pair[[2]]) %>% 
      select(sample_id)
    sample2 = sample2$sample_id

    # Skip if either timepoint is missing for this patient
    if (length(sample1) == 0 | length(sample2) == 0) {
      return(NULL)
    }
    
    # Skip if diversity data is missing for either sample
    if (!(sample1 %in% rownames(clonal_diversity_df)) | !(sample2 %in% rownames(clonal_diversity_df))) {
      return(NULL)
    }
    
    # Extract diversity values for both timepoints
    clonal_diversity_sample_1 = clonal_diversity_df[sample1, "clonotype_diversity"]
    clonal_diversity_sample_2 = clonal_diversity_df[sample2, "clonotype_diversity"]
    
    # Calculate expansion metric based on method
    if (method == "division") {
      # Ratio method: later / earlier
      # Values < 1 indicate expansion (diversity decreased)
      clonal_expansion <- clonal_diversity_sample_2 / clonal_diversity_sample_1
    } else {
      # Difference method: later - earlier
      # Negative values indicate expansion (diversity decreased)
      clonal_expansion <- clonal_diversity_sample_2 - clonal_diversity_sample_1
    }
    
    return(clonal_expansion)
  }
  
  # Initialize results dataframe
  clonal_expansion_df = data.frame(matrix(ncol = 1, nrow = length(patients)))
  rownames(clonal_expansion_df) = patients
  colnames(clonal_expansion_df) = "clonal_expansion"
  
  # Calculate expansion for each patient
  for (i in 1:length(patients)) {
    patient = patients[i]
    out = get_clonal_expansion_per_patient(patient = patient)
    if (!is.null(out)) {
      clonal_expansion_df[i, "clonal_expansion"] = out
    }
  }
  
  # Save results to file
  outfile = paste0(outdir, "/", compared_pair[[1]], "_vs_", compared_pair[[2]], ".txt")
  outfile = gsub(" ", "_", outfile)
  write.table(clonal_expansion_df, outfile, sep="\t", quote = FALSE, row.names = TRUE)
}

################################################################################
# FUNCTION: run_example
################################################################################
# Execute expansion calculation for a specific cell type and diversity type
#
# PARAMETERS:
#   celltype - T cell subpopulation name
#   diversity_type - "shannon" or "simpson"
#
# WORKFLOW:
#   1. Load diversity index file from Step 2
#   2. Load sample metadata
#   3. Calculate expansion (C1 vs C2 comparison)
#   4. Save results per patient
#
################################################################################
run_example = function(celltype, diversity_type) {
  # Load diversity index file from Step 2
  clonal_diversity_df_file = paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/", diversity_type, "_clonal_diversity_", celltype, "_T_cells.txt")
  
  # Set output directory
  outdir = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/"
  outdir = paste0(outdir, "C1_vs_C2/") 
  dir.create(outdir)
  setwd(outdir)
  
  # Load metadata (patient-sample-timepoint mappings)
  metadata = get_metadata_sorted()
  
  # Load diversity index data
  clonal_diversity_df <- read.delim(clonal_diversity_df_file, check.names = FALSE)
  
  # Calculate expansion by comparing C1 vs C2 timepoints
  # C1 = Cycle 1 (earlier), C2 = Cycle 2 (later)
  out = get_clonal_expansion_timepoint_comparison(
    compared_pair = c("C1", "C2"),
    metadata = metadata,
    clonal_diversity_df = clonal_diversity_df,
    outdir = paste0("clonal_expansion_", diversity_type, "_", celltype, "_T_cells")
  )
}

################################################################################
# SECTION: PROCESS ALL T CELL SUBPOPULATIONS
################################################################################
# Define the same T cell subpopulation mappings as in previous steps
# Ensures consistency across the entire analysis pipeline

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

celltypes <- unique(celltype_to_cluster$celltype)

# Load metadata utility function
source("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/get_metadata.R")

# Calculate expansion for all combinations of cell type × diversity index
for (celltype in celltypes) {
  for (diversity_type in c("shannon", "simpson")) {
    run_example(celltype, diversity_type) 
  }
}

################################################################################
# END OF SCRIPT
################################################################################
# 
# OUTPUT SUMMARY:
# For each combination of cell type and diversity index:
#   Directory: clonal_expansion_{shannon|simpson}_{celltype}_T_cells_division/
#   File: C1_vs_C2.txt
#   Content: Expansion metric per patient (rows=patients, col=clonal_expansion)
#
# INTERPRETATION:
#   - Ratio < 1: Clonal expansion occurred (diversity decreased from C1 to C2)
#   - Ratio > 1: Clonal diversification occurred (diversity increased)
#   - The magnitude indicates strength of expansion/diversification
#
# These expansion metrics serve as input for:
#   - Step 4: clonal_expansion_plots_4.R (Figure 6d visualization)
#   - Statistical analyses correlating expansion with clinical outcomes
#   - Comparison between responders vs non-responders
#
################################################################################
