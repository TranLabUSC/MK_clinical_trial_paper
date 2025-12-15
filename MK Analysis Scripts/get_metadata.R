################################################################################
# get_metadata.R - Utility Functions for Metadata Processing
################################################################################
#
# PURPOSE:
# Supporting script for clonal expansion analysis providing utility functions
# for standardizing timepoint nomenclature and processing sample metadata.
#
# MAIN FUNCTIONS:
# 1. get_metadata_sorted() - Extracts and standardizes timepoint labels from
#    sample IDs, ensuring consistent ordering across analyses
# 2. to_bool() - Converts "Y"/"N" strings to TRUE/FALSE boolean values
#
# USAGE:
# This script is typically sourced by other analysis scripts that need to:
# - Parse timepoint information from sample identifiers
# - Standardize timepoint ordering (S/Pre → C1 → C2 → ... → C36)
# - Convert categorical metadata to boolean format
#
# TIMEPOINT NOMENCLATURE MAPPING:
# - "S" or "Pre" = Baseline/Pre-treatment sample
# - "C1" through "C36" = Treatment cycle numbers (C1 = Cycle 1, etc.)
# - Standardized order ensures consistent temporal analysis across all scripts
#
################################################################################

library(dplyr)

# Define standardized timepoint order for consistent temporal analysis
# "S" represents baseline/pre-treatment samples (also called "Pre" in some contexts)
# C1-C36 represent treatment cycle numbers in chronological order
sample_order <- c(
  "S",      # Baseline/Pre-treatment
  "C1",     # Cycle 1
  "C2",     # Cycle 2
  "C4",     # Cycle 4
  "C6",     # Cycle 6
  "C9",     # Cycle 9
  "C18",    # Cycle 18
  "C36"     # Cycle 36
)


################################################################################
# FUNCTION: get_metadata_sorted
################################################################################
#
# PURPOSE:
# Loads and processes sample metadata, extracting standardized timepoint labels
# from sample identifiers and sorting samples chronologically within each patient.
#
# PARAMETERS:
# - meta_file: Path to CSV file containing sample metadata (default: MK2_aggregate.csv)
#
# WORKFLOW:
# 1. Read metadata CSV file with sample information
# 2. Extract timepoint labels by removing patient ID prefix from sample origin
#    - Example: "8C1" (Patient 8, Cycle 1) → "C1"
#    - Example: "3S" (Patient 3, Baseline) → "S"
# 3. Convert timepoint to ordered factor using sample_order
#    - Ensures consistent temporal ordering: S < C1 < C2 < C4 < ... < C36
# 4. Sort metadata by Patient ID and timepoint chronologically
#
# RETURNS:
# Data frame with processed metadata containing:
# - All original columns from input CSV
# - New 'timepoint' column with standardized, ordered factor values
# - Rows sorted by Patient (ascending) then timepoint (chronological)
#
# BIOLOGICAL RATIONALE:
# Standardized timepoint ordering is critical for:
# - Temporal trajectory analysis (tracking changes over treatment course)
# - Paired statistical comparisons (e.g., C1 vs Pre, C2 vs C1)
# - Clonal expansion calculations (requires consistent time ordering)
# - Visualization (plotting samples in correct chronological sequence)
#
################################################################################
get_metadata_sorted <- function(meta_file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/MK2_aggregate.csv") {
  # Load metadata CSV
  combined_meta <- read.csv(meta_file)
  
  # Extract timepoint by removing patient ID prefix from origin column
  # mapply applies the substitution function to each (origin, Patient) pair
  combined_meta$timepoint <- mapply(
    function(orig, pat) sub(paste0("^", pat), "", orig),  # Remove "PatientID" prefix
    combined_meta$origin,
    combined_meta$Patient
  )
  
  # Convert timepoint to ordered factor for chronological sorting
  combined_meta$timepoint <- factor(combined_meta$timepoint, levels = sample_order)
  
  # Sort metadata by Patient ID (primary) and timepoint (secondary)
  metadata <- combined_meta %>% arrange(Patient, timepoint)
  
  return(metadata)
}

################################################################################
# FUNCTION: to_bool
################################################################################
#
# PURPOSE:
# Converts string-based yes/no indicators to boolean TRUE/FALSE values.
#
# PARAMETERS:
# - x: Character string, typically "Y" or "N"
#
# RETURNS:
# - TRUE if x == "Y"
# - FALSE otherwise (including "N", NA, or any other value)
#
# USAGE:
# Commonly used to convert categorical metadata fields (e.g., "IDH_mutant": "Y"/"N")
# into boolean format suitable for logical operations and filtering.
#
# EXAMPLE:
# - to_bool("Y") → TRUE (patient has IDH mutation)
# - to_bool("N") → FALSE (patient is IDH wildtype)
#
################################################################################
to_bool <- function(x) {
  if (x == "Y") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}