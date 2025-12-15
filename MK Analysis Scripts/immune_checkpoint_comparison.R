#########################################################################################
# IMMUNE CHECKPOINT EXPRESSION COMPARISON ACROSS COHORTS AND TIMEPOINTS
#########################################################################################
#
# PURPOSE:
# Compare expression of 21 immune checkpoint genes across patient cohorts
# (Control, Short-term survivors, Long-term survivors) and treatment timepoints
# (Pre-treatment baseline, Cycle 1, Cycle 2) in T cells from glioblastoma patients
# receiving anti-PD-1 (MK-3475) + MLA checkpoint blockade therapy.
#
# GENERATES MANUSCRIPT FIGURE:
# - Figure 8c: 21-page PDF with 2D kernel density contours showing spatial distribution
#              of immune checkpoint expression in UMAP space, colored by cohort
#
# BIOLOGICAL RATIONALE:
# Immune checkpoints are regulatory molecules that modulate T cell activation and function.
# In cancer, tumor cells exploit these pathways to evade immune destruction. Multiple
# checkpoint proteins work synergistically:
# - Inhibitory checkpoints (PD-1, CTLA-4, LAG3, TIM-3, TIGIT, BTLA, VISTA): suppress
#   T cell responses and contribute to exhaustion phenotypes
# - Co-stimulatory checkpoints (CD28, ICOS, 4-1BB, OX40, GITR, B7-H3): enhance T cell
#   activation and survival
# - Ligands (PD-L1, PD-L2, B7-H3, B7-H4, HHLA2): expressed on tumor/immune cells
# - Metabolic checkpoints (IDO1, IDO2, CD47): regulate immunosuppressive metabolism
#
# Comparing checkpoint expression patterns between survivor cohorts reveals:
# 1. Baseline checkpoint profiles predictive of therapeutic response
# 2. Evolution of checkpoint expression during anti-PD-1 therapy
# 3. Compensatory upregulation of alternative checkpoints post-treatment
# 4. Spatial localization of checkpoint-expressing cells in tumor microenvironment
#
# CLINICAL SIGNIFICANCE:
# Understanding checkpoint expression dynamics guides:
# - Patient stratification for checkpoint blockade therapy
# - Identification of resistance mechanisms
# - Selection of combination therapy targets
# - Biomarker development for treatment monitoring
#
# WORKFLOW:
# SECTION 1 (Lines 45-254): FeaturePlot Generation & Data Extraction
#   1. Load survival data and assign cohorts (Control/ShortTerm/LongTerm)
#   2. Load T cell Seurat object (all CD4+ and CD8+ subpopulations)
#   3. Subset to matching patients and timepoints (Pre, C1, C2)
#   4. Downsample: equal cells across timepoints within each cohort (prevents bias)
#   5. Generate FeaturePlots for 21 checkpoint genes split by cohort×timepoint
#   6. Extract underlying data (UMAP coordinates + expression) for density analysis
#   7. Export: 21-row PDF (genes) × 9 columns (3 cohorts × 3 timepoints)
#
# SECTION 2 (Lines 256-358): 2D Kernel Density Visualization
#   1. Load FeaturePlot data from Section 1
#   2. For each gene: filter to top 10% expressers (reduce noise from low-expressers)
#   3. Generate 2D density contours weighted by expression level
#   4. Three panels per gene (Pre, C1, C2) with cohort-specific colored contours
#   5. Gray background shows all cells; colored contours show high-expressers
#   6. Export: 21-page PDF (one gene per page, 3 timepoint panels per page)
#
# STATISTICAL APPROACH:
# - Downsampling: Equal representation across timepoints (min cells across Pre/C1/C2)
# - Density estimation: Kernel density with expression-weighted contribution
# - Expression filtering: Top 10% (90th percentile) to focus on biologically relevant cells
# - Adaptive binning: Contour count scales with cell number (3-10 bins)
#
# OUTPUT FILES:
# 1. FeaturePlots PDF: 21 genes × 9 cohort-timepoint combinations
# 2. 2D Density PDF: 21 pages (Figure 8c), each with 3 timepoint panels
# 3. Intermediate RDS: UMAP coordinates + expression for all cells/genes
#
# DEPENDENCIES:
# - Seurat: Feature extraction, subsetting, plotting
# - dplyr: Data manipulation, filtering
# - ggplot2: Base plotting framework
# - cowplot: Multi-panel figure assembly
#
# INPUT DATA:
# - T cell Seurat object: MK_T_Cells_seurat_obj.RDS (all CD4+/CD8+ clusters)
# - Survival metadata: patient_survival_data.csv (cohort assignments, OS)
#
# TECHNICAL NOTES:
# - Cluster 13 excluded (likely doublets or low-quality cells)
# - UMAP coordinates: visualization space for density contours
# - Expression values: normalized counts from Seurat's RNA assay
# - Downsampling seed: 123 (reproducibility)
# - FeaturePlot percentiles: q10-q90 (10th-90th percentile cutoffs)
#
#########################################################################################

#########################################################################################
# SECTION 1: FEATUREPLOT GENERATION & DATA EXTRACTION
#########################################################################################

# Load necessary libraries for Seurat analysis, data manipulation, and visualization
library(Seurat)    # Single-cell analysis framework
library(dplyr)     # Data manipulation (filter, mutate, join)
library(ggplot2)   # Grammar of graphics plotting
library(cowplot)   # Multi-panel figure assembly

# -----------------------
# 1. Define Survivor Groups
# -----------------------
control_group <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# -----------------------
# 2. Load & Prepare Survival Data
# -----------------------
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df$patient_id <- as.integer(survival_df$patient_id)

survival_df <- survival_df %>%
  filter(site == "UF") %>%
  mutate(SurvivalGroup = case_when(
    patient_id %in% short_term_survivor_group ~ "ShortTerm",
    patient_id %in% long_term_survivor_group  ~ "LongTerm",
    patient_id %in% control_group             ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(SurvivalGroup))

required_columns <- c("OS.months.", "patient_id")
missing_columns <- setdiff(required_columns, colnames(survival_df))
if (length(missing_columns) > 0) {
  stop(paste("ERROR: Missing columns in survival_df:", paste(missing_columns, collapse = ", ")))
}
survival_df <- survival_df %>% arrange(OS.months.)

# -----------------------
# 3. LOAD T CELL SEURAT OBJECT AND SUBSET TO ANALYSIS COHORT
# -----------------------
# Load pre-processed T cell Seurat object containing all CD4+ and CD8+ T cell
# subpopulations with UMAP embeddings and cluster annotations (resolution 1.0).
# This object includes cells from all timepoints across the clinical trial.

seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)

# QUALITY CONTROL: Remove cluster 13 (likely doublets or low-quality cells)
# This cluster was identified as problematic during initial QC and is excluded
# from all downstream analyses to prevent technical artifacts
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters != 13)

# Validate required metadata column exists
if (!"Patient" %in% colnames(seurat_object_t_cells@meta.data)) {
  stop("ERROR: 'Patient' column not found in Seurat object metadata.")
}

# Find overlapping patients between scRNA-seq data and survival metadata
# This ensures we only analyze patients with complete clinical annotations
common_patients <- intersect(
  unique(seurat_object_t_cells@meta.data$Patient),
  unique(survival_df$patient_id)
)
if (length(common_patients) == 0) {
  stop("ERROR: No overlapping patients found between Seurat object and survival data.")
} else {
  message(paste("Number of overlapping patients:", length(common_patients)))
}

# Subset Seurat object to only patients in our defined cohorts (Control/ShortTerm/LongTerm)
# This removes patients without survival annotations or outside UF site
patients_to_keep <- survival_df$patient_id
seurat_subset <- subset(
  seurat_object_t_cells,
  cells = WhichCells(seurat_object_t_cells, expression = Patient %in% patients_to_keep)
)

message(paste("Number of patients after subsetting:", length(unique(seurat_subset@meta.data$Patient))))
        
# -----------------------
# 4. ANNOTATE CELLS WITH SURVIVAL COHORT INFORMATION
# -----------------------
# Add SurvivalGroup metadata to each cell based on patient assignment.
# This enables cohort-level comparisons in downstream visualizations.

# Preserve cell barcodes as rownames (required for Seurat metadata structure)
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  mutate(CellBarcode = rownames(.))

# Create mapping table: Patient ID → Survival Cohort
mapping_df <- survival_df %>%
  select(patient_id, SurvivalGroup) %>%
  rename(Patient = patient_id)  # Match column name in Seurat metadata

# Left join: Assign each cell its patient's survival cohort
seurat_subset@meta.data <- seurat_subset@meta.data %>%
  left_join(mapping_df, by = "Patient")

# Restore rownames (critical for Seurat object integrity)
rownames(seurat_subset@meta.data) <- seurat_subset@meta.data$CellBarcode
seurat_subset@meta.data$CellBarcode <- NULL

# Convert to ordered factor for consistent visualization order (Control → ShortTerm → LongTerm)
seurat_subset@meta.data$SurvivalGroup <- factor(
  seurat_subset@meta.data$SurvivalGroup,
  levels = c("Control", "ShortTerm", "LongTerm")
)

# -----------------------
# 5. FILTER TO KEY TIMEPOINTS AND PERFORM BALANCED DOWNSAMPLING
# -----------------------
# TIMEPOINT SELECTION: Focus on Pre (baseline), C1 (first treatment), C2 (second treatment)
# These capture: pre-treatment state, early response, and sustained response patterns
#
# DOWNSAMPLING RATIONALE:
# Different timepoints often have unequal cell counts due to:
# 1. Sample availability (not all patients sampled at all timepoints)
# 2. Sequencing depth variation
# 3. Cell viability differences across timepoints
#
# Without downsampling, timepoints with more cells would dominate density plots
# and statistical comparisons, creating bias. Equal representation ensures fair
# comparison of checkpoint expression evolution across treatment timeline.

valid_timepoints <- c("Pre","C1","C2")  # Baseline and first two treatment cycles
seurat_subset <- subset(seurat_subset, subset = TimePoint %in% valid_timepoints)

# DOWNSAMPLING STRATEGY: Equal cells per timepoint within each cohort
# For each cohort (Control/ShortTerm/LongTerm):
#   1. Count cells at Pre, C1, C2
#   2. Find minimum count (n_min)
#   3. Randomly sample n_min cells from each timepoint
# This ensures equal representation: Control_Pre = Control_C1 = Control_C2, etc.

set.seed(123)  # Reproducibility for random sampling
final_cells <- c()

for(grp in levels(seurat_subset@meta.data$SurvivalGroup)) {
  # Get all cells for this cohort
  grp_cells <- WhichCells(seurat_subset, expression = SurvivalGroup == grp)
  
  # Separate by timepoint
  pre_cells <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "Pre"))
  c1_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C1"))
  c2_cells  <- intersect(grp_cells, WhichCells(seurat_subset, expression = TimePoint == "C2"))
  
  # Count cells per timepoint
  n_pre <- length(pre_cells)
  n_c1  <- length(c1_cells)
  n_c2  <- length(c2_cells)
  
  # Find minimum count (limiting timepoint)
  n_min <- min(n_pre, n_c1, n_c2)
  if(n_min == 0) {
    warning(paste("Skipping group", grp, " - zero cells in at least one timepoint."))
    next
  }
  
  # Randomly sample n_min cells from each timepoint
  final_cells <- c(
    final_cells,
    sample(pre_cells, n_min),  # Pre-treatment baseline
    sample(c1_cells, n_min),   # Cycle 1 (early response)
    sample(c2_cells, n_min)    # Cycle 2 (sustained response)
  )
}

# Create downsampled Seurat object with balanced timepoint representation
seurat_timepoint_downsampled <- subset(seurat_subset, cells = final_cells)

# Verify balanced counts (should be equal within each cohort)
cat("Final cell counts after timepoint-downsampling:\n")
print(table(
  seurat_timepoint_downsampled@meta.data$SurvivalGroup,
  seurat_timepoint_downsampled@meta.data$TimePoint
))

# -----------------------
# 6. CREATE COMBINED COHORT × TIMEPOINT METADATA FOR SPLIT PLOTTING
# -----------------------
# Combine SurvivalGroup and TimePoint into a single factor for FeaturePlot splitting.
# This creates 9 distinct categories (3 cohorts × 3 timepoints) enabling side-by-side
# comparison of checkpoint expression patterns across all conditions simultaneously.
#
# VISUALIZATION STRATEGY:
# Instead of separate plots for each cohort or timepoint, split.by="Group_TimePoint"
# generates a 9-panel grid showing how checkpoint expression evolves within each cohort
# across treatment timeline: Control Pre→C1→C2, ShortTerm Pre→C1→C2, LongTerm Pre→C1→C2

seurat_timepoint_downsampled$Group_TimePoint <- paste(
  seurat_timepoint_downsampled$SurvivalGroup,
  seurat_timepoint_downsampled$TimePoint,
  sep="_"
)

# Define explicit factor order for consistent left-to-right layout in plots
# Order: Control (baseline) → ShortTerm (non-responders) → LongTerm (responders)
# Within each: Pre (baseline) → C1 (early) → C2 (sustained)
desired_order <- c(
  "Control_Pre","Control_C1","Control_C2",       # Standard chemotherapy timeline
  "ShortTerm_Pre","ShortTerm_C1","ShortTerm_C2", # Poor responder timeline
  "LongTerm_Pre","LongTerm_C1","LongTerm_C2"     # Good responder timeline
)
seurat_timepoint_downsampled$Group_TimePoint <- factor(
  seurat_timepoint_downsampled$Group_TimePoint,
  levels=desired_order
)

# -----------------------
# 7. DEFINE IMMUNE CHECKPOINT GENE PANEL (21 GENES)
# -----------------------
# Comprehensive panel covering major checkpoint pathways validated in cancer immunotherapy.
# Includes inhibitory checkpoints, co-stimulatory molecules, ligands, and metabolic regulators.
#
# GENE CATEGORIES:
#
# INHIBITORY CHECKPOINTS (suppress T cell function):
#   PDCD1 (PD-1): Primary target of checkpoint blockade therapy
#   CTLA4: Competes with CD28 for B7 binding, early checkpoint
#   LAG3: Binds MHC-II, synergizes with PD-1 in exhaustion
#   HAVCR2 (TIM-3): T cell immunoglobulin mucin-3, exhaustion marker
#   TIGIT: T cell immunoreceptor with Ig and ITIM domains
#   BTLA: B and T lymphocyte attenuator, HVEM ligand
#   VSIR (VISTA): V-domain Ig suppressor of T cell activation
#
# CHECKPOINT LIGANDS (expressed on tumor/immune cells):
#   CD274 (PD-L1): Primary PD-1 ligand, predictive biomarker
#   PDCD1LG2 (PD-L2): Alternative PD-1 ligand
#   CD276 (B7-H3): Immunoglobulin superfamily, dual function
#   VTCN1 (B7-H4): Coinhibitory ligand
#   NCR3LG1 (B7-H6): NK cell activating ligand
#   HHLA2: HERV-H LTR-associating 2, alternative checkpoint
#
# CO-STIMULATORY MOLECULES (enhance T cell activation):
#   CD28: Primary co-stimulator, binds B7-1/B7-2
#   ICOS: Inducible T cell co-stimulator, memory formation
#   TNFRSF4 (OX40): TNF receptor, T cell survival
#   TNFRSF9 (4-1BB): TNF receptor, effector function
#   TNFRSF18 (GITR): Glucocorticoid-induced TNFR
#
# METABOLIC/IMMUNOSUPPRESSIVE:
#   IDO1: Indoleamine 2,3-dioxygenase 1, tryptophan catabolism
#   IDO2: Indoleamine 2,3-dioxygenase 2
#   CD47: "Don't eat me" signal, macrophage checkpoint
#
# THERAPEUTIC RELEVANCE:
# - PD-1/PD-L1: Current trial target (MK-3475 anti-PD-1)
# - Multiple genes: Potential combination therapy targets
# - Expression patterns: Predict response vs resistance mechanisms

immune_checkpoint_genes <- c(
  "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",  # Core checkpoints
  "BTLA","VSIR","CD276","VTCN1","IDO1","IDO2","CD47",          # Extended panel
  "CD28","ICOS","TNFRSF4","TNFRSF9","TNFRSF18",                # Co-stimulatory
  "NCR3LG1","HHLA2"                                             # Additional ligands
)

# -----------------------
# 8. GENERATE FEATUREPLOTS AND EXTRACT UNDERLYING DATA
# -----------------------
# Create comprehensive FeaturePlot visualizations for all 21 checkpoint genes
# across all 9 cohort×timepoint combinations, then extract the raw data for
# subsequent 2D density analysis (Section 2).
#
# VISUALIZATION STRUCTURE:
# - 21-row PDF (one row per gene)
# - 9 columns per row (3 cohorts × 3 timepoints)
# - Color scale: lightgrey (low expression) → red (high expression)
# - Expression cutoffs: q10-q90 (10th-90th percentile) to reduce outlier effects
#
# DATA EXTRACTION STRATEGY:
# FeaturePlot with combine=FALSE returns a list of ggplot objects, one per split panel.
# Each ggplot contains a $data dataframe with:
#   - UMAP_1, UMAP_2: Coordinates for visualization
#   - <gene>: Expression values (normalized counts)
#   - Cell barcodes: Rownames linking back to metadata
# We extract and annotate this data with split_label (Group_TimePoint) for density plots.

# Create top row labels for the 9-column grid
top_labels <- levels(seurat_timepoint_downsampled$Group_TimePoint)
top_label_plots <- lapply(top_labels, function(txt) {
  ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=txt), size=4)
})
top_row <- plot_grid(plotlist=top_label_plots, nrow=1)

# Initialize storage for plots and extracted data
all_rows_list  <- list()  # Stores one row per gene for final PDF assembly
all_feature_df <- list()  # Stores extracted data frames for all genes

# Loop through each checkpoint gene
for (gene in immune_checkpoint_genes) {
  cat("Default Assay is:", DefaultAssay(seurat_timepoint_downsampled), "\n")
  
  # GENERATE FEATUREPLOT with combine=FALSE
  # This returns a list of 9 ggplot objects (one per Group_TimePoint split)
  # Each plot shows expression of the current gene in UMAP space
  p_list <- FeaturePlot(
    object    = seurat_timepoint_downsampled,
    features  = gene,                   # Current checkpoint gene
    split.by  = "Group_TimePoint",      # Split into 9 panels
    cols      = c("lightgrey","red"),   # Low → High expression gradient
    combine   = FALSE,                  # Return list of plots (not combined)
    pt.size   = 0.6,                    # Point size for single cells
    order     = TRUE,                   # Plot high-expressers on top (visibility)
    min.cutoff= "q10",                  # 10th percentile minimum (reduce noise)
    max.cutoff= "q90"                   # 90th percentile maximum (reduce outliers)
  )
  
  # EXTRACT UNDERLYING DATA FROM EACH SPLIT PANEL
  # Each element of p_list is a ggplot object with a $data dataframe
  # containing UMAP coordinates and expression values for that panel's cells
  splitted_data_dfs <- lapply(seq_along(p_list), function(i) {
    # Extract data frame from ggplot object
    df_i <- p_list[[i]]$data  # Columns: UMAP_1, UMAP_2, <gene>, etc.
    df_i$panel_index <- i     # Track which panel (1-9) this data came from
    df_i$gene        <- gene  # Add gene name for combined dataset
    
    # DETERMINE WHICH GROUP_TIMEPOINT THIS PANEL REPRESENTS
    # Rownames of df_i are cell barcodes - look them up in Seurat metadata
    these_cells <- rownames(df_i)
    
    # All cells in a split panel should have the same Group_TimePoint value
    # (that's what split.by does - separates cells into homogeneous groups)
    gtp <- unique(seurat_timepoint_downsampled$Group_TimePoint[these_cells])
    
    if (length(gtp) == 1) {
      # Expected case: all cells belong to one Group_TimePoint
      df_i$split_label <- as.character(gtp)  # e.g., "Control_Pre", "LongTerm_C2"
    } else {
      # Unexpected case: multiple or zero Group_TimePoint values found
      # This shouldn't happen with proper split.by, but handle gracefully
      df_i$split_label <- NA
    }
    
    return(df_i)
  })
  
  # Combine all 9 panels for this gene into a single dataframe
  splitted_data_combined <- dplyr::bind_rows(splitted_data_dfs)
  all_feature_df[[gene]] <- splitted_data_combined  # Store for later density analysis
  
  # PREPARE PLOTS FOR PDF EXPORT
  # Apply cosmetic modifications: flip coordinates, remove titles/legends, minimal theme
  p_list <- lapply(p_list, function(x) {
    x + coord_flip() + ggtitle(NULL) + NoLegend() + theme_void()
  })
  
  # Arrange 9 panels in a single horizontal row
  gene_row_plots <- plot_grid(plotlist=p_list, nrow=1)
  
  # Create gene name label (vertical text on left side)
  gene_label_plot <- ggplot() + theme_void() +
    geom_text(aes(x=0.5, y=0.5, label=gene), angle=90, size=5)
  
  # Combine gene label + 9 panels into final row
  row_for_gene <- plot_grid(
    gene_label_plot,  # Left: gene name (4% width)
    gene_row_plots,   # Right: 9 panels (96% width)
    nrow=1, rel_widths=c(0.04,1)
  )
  
  # Store this row for final PDF assembly
  all_rows_list[[gene]] <- row_for_gene
}

# ASSEMBLE FINAL MULTI-PAGE FIGURE
# Combine all gene rows (21) into a single vertical matrix
main_matrix <- plot_grid(plotlist=all_rows_list, ncol=1)

# Add top row with cohort×timepoint labels
final_plot <- plot_grid(top_row, main_matrix, ncol=1, rel_heights=c(0.05,1))

# EXPORT COMPREHENSIVE FEATUREPLOT PDF
# Dimensions: 20 inches wide × 42 inches tall (2×21 genes)
# Layout: 21 rows (genes) × 9 columns (cohort×timepoint combinations)
output_pdf_original <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/immune_checkpoint_timepoint_splitted_featureplots.pdf"
num_genes <- length(immune_checkpoint_genes)  # 21
pdf_width  <- 2*10   # 20 inches (accommodates 9 columns)
pdf_height <- 2*num_genes  # 42 inches (21 rows × 2 inches each)

pdf(output_pdf_original, width=pdf_width, height=pdf_height)
print(final_plot)
dev.off()
cat("Saved FeaturePlots to:", output_pdf_original, "\n")

# EXPORT EXTRACTED DATA FOR DENSITY ANALYSIS
# Combine all extracted data frames (21 genes × 9 panels) into single RDS file
# This data will be used in Section 2 for 2D density contour visualization
featureplot_data <- dplyr::bind_rows(all_feature_df)
saveRDS(featureplot_data,
        file="/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds"
)
cat("Saved underlying FeaturePlot data to: /project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds\n")
        
        


#########################################################################################
# SECTION 2: 2D KERNEL DENSITY VISUALIZATION (FIGURE 8c)
#########################################################################################
#
# PURPOSE:
# Transform the FeaturePlot data from Section 1 into publication-quality 2D kernel
# density contour plots showing spatial distribution of immune checkpoint expression
# in UMAP space. This visualization method reveals regional enrichment patterns that
# simple scatter plots obscure, particularly for moderately-expressed genes.
#
# VISUALIZATION STRATEGY:
# - Filter to top 10% expressers per gene (reduces noise from low/absent expression)
# - Generate expression-weighted 2D density contours using kernel density estimation
# - Three panels per gene (Pre, C1, C2) with cohort-specific colored contours
# - Gray background displays all cells for spatial context
# - Colored contours (yellow=Control, blue=ShortTerm, red=LongTerm) show high-expressers
#
# BIOLOGICAL INTERPRETATION:
# Density contours reveal:
# 1. Spatial localization: Which T cell subpopulations express each checkpoint
# 2. Cohort differences: Differential expression patterns between responders/non-responders
# 3. Temporal dynamics: How expression patterns evolve across treatment timeline
# 4. Co-localization: Whether different checkpoints mark the same or distinct cell states
#
# ADVANTAGES OVER SCATTER PLOTS:
# - Reduces visual clutter from thousands of overlapping points
# - Highlights density patterns (where expression is concentrated)
# - Weighted contours reflect both cell count AND expression level
# - More informative for moderately-expressed genes (avoids "salt and pepper" appearance)
#
# OUTPUT:
# 21-page PDF (one gene per page), each page with 3 panels (Pre, C1, C2)
# This becomes Manuscript Figure 8c
#
################################################################################

################################################################################
# 1. LOAD LIBRARIES AND INPUT DATA
################################################################################
library(dplyr)     # Data manipulation
library(ggplot2)   # Plotting framework
library(cowplot)   # Multi-panel assembly

# Load extracted data from Section 1 (UMAP coordinates + expression values)
# columns: UMAP_1, UMAP_2, gene, split_label (e.g., "Control_Pre"), rna_<gene>
df_all <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/featureplot_data.rds")

# Define output path for 2D density PDF (Figure 8c)
pdf_output <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/NC_Monocytes_Analysis/ImmuneCheckpoint_2D_Density_from_FeaturePlotData.pdf"


################################################################################
# 2. DEFINE 2D DENSITY PLOTTING FUNCTION
################################################################################
#
# FUNCTION: plot_2d_density()
# Creates a 3-panel figure (Pre, C1, C2) with expression-weighted density contours
# for a single checkpoint gene, colored by survival cohort.
#
# PARAMETERS:
# - df_all: Complete FeaturePlot data from Section 1 (all genes, all cells)
# - target_gene: Name of checkpoint gene to visualize (e.g., "PDCD1")
#
# PROCESSING STEPS:
# 1. Subset to target gene's data
# 2. Identify expression column (handles special case for IDO1)
# 3. Filter to top 10% expressers (90th percentile cutoff)
# 4. Parse timepoint and cohort from split_label
# 5. Generate 3 timepoint panels with cohort-colored density overlays
#
# VISUALIZATION COMPONENTS PER PANEL:
# - Gray background: All cells (spatial context)
# - Yellow contours: Control cohort high-expressers
# - Blue contours: ShortTerm cohort high-expressers
# - Red contours: LongTerm cohort high-expressers
# - Contour density weighted by expression level (not just cell count)
# - Adaptive binning: 3-10 contour lines based on cell count
#
# RETURNS:
# A cowplot grid with 3 panels (Pre, C1, C2) arranged horizontally
#
plot_2d_density <- function(df_all, target_gene) {
  
  # STEP 1: Subset to current gene's data
  df_gene <- df_all %>% filter(gene == target_gene)
  
  # STEP 2: Identify expression column name
  # Most genes: "rna_<GENE>" (e.g., "rna_PDCD1")
  # Exception: IDO1 stored as "IDO1" without "rna_" prefix
  if (target_gene == "IDO1") {
    expr_col <- "IDO1"
  } else {
    expr_col <- paste0("rna_", target_gene)
  }
  
  # Validate column exists
  if (! expr_col %in% colnames(df_gene)) {
    stop("No column named ", expr_col, " for gene=", target_gene)
  }
  
  # Copy expression values to standardized "feature" column for downstream use
  df_gene$feature <- df_gene[[expr_col]]
  
  # STEP 3: Create two datasets for layered visualization
  # df_all_cells: All cells (gray background showing spatial context)
  df_all_cells <- df_gene
  
  # df_gene: Only top 10% expressers (colored density contours)
  # RATIONALE FOR 90TH PERCENTILE CUTOFF:
  # - Reduces noise from cells with low/absent expression
  # - Focuses on biologically meaningful expression levels
  # - Prevents density contours dominated by zero-inflation
  # - More interpretable visualization of true checkpoint-expressing populations
  expr_threshold <- quantile(df_gene$feature, probs = 0.9, na.rm = TRUE)
  df_gene <- df_gene %>% filter(feature > expr_threshold)
  
  # STEP 4: Parse metadata from split_label
  # split_label format: "SurvivalGroup_TimePoint" (e.g., "Control_Pre", "LongTerm_C2")
  # Extract timepoint (after underscore) and cohort (before underscore)
  df_all_cells$TimePoint <- sub("^.+_", "", df_all_cells$split_label)  # "Pre", "C1", or "C2"
  df_all_cells$Group     <- sub("_.*$", "", df_all_cells$split_label)  # "Control", "ShortTerm", or "LongTerm"
  df_gene$TimePoint      <- sub("^.+_", "", df_gene$split_label)
  df_gene$Group          <- sub("_.*$", "", df_gene$split_label)
  
  # STEP 5: Generate density plots for each timepoint
  timepoints   <- c("Pre", "C1", "C2")  # Treatment timeline
  group_colors <- c("Control" = "yellow", "ShortTerm" = "blue", "LongTerm" = "red")
  
  p_list <- list()
  for (tp in timepoints) {
    # Subset to current timepoint
    df_tp_all <- df_all_cells %>% filter(TimePoint == tp)  # All cells (background)
    df_tp_top <- df_gene %>% filter(TimePoint == tp)       # Top 10% expressers (contours)
    
    # CREATE BASE PLOT: Gray background showing all cells in UMAP space
    p <- ggplot(df_tp_all, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(color = "gray70", size = 0.4, alpha = 0.5) +  # Low alpha for subtle background
      coord_fixed() +    # Preserve UMAP aspect ratio
      coord_flip() +     # Rotate for consistent orientation
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank()) +  # Remove grid lines (cleaner visualization)
      ggtitle(paste(target_gene, "-", tp))
    
    # ADD DENSITY CONTOURS: One layer per survival cohort
    for (grp in c("Control", "ShortTerm", "LongTerm")) {
      df_grp <- filter(df_tp_top, Group == grp)
      
      # Need at least 2 cells to compute density
      if (nrow(df_grp) >= 2) {
        # ADAPTIVE BINNING: Scale contour count with cell number
        # - Fewer cells → fewer bins (3 minimum, avoids over-smoothing)
        # - More cells → more bins (10 maximum, avoids clutter)
        # - Heuristic: 1 bin per 20 cells, bounded [3, 10]
        n_cells <- nrow(df_grp)
        n_bins <- max(3, min(10, floor(n_cells / 20)))
        
        # EXPRESSION-WEIGHTED DENSITY:
        # aes(weight = feature) gives cells with higher expression more influence
        # Result: contours concentrated where expression is BOTH dense AND high
        p <- p + geom_density_2d(
          data = df_grp,
          aes(weight = feature),           # Weight by expression level
          bins = n_bins,                   # Adaptive contour count
          color = group_colors[grp],       # Cohort-specific color
          size  = 0.7                      # Line thickness
        )
      }
    }
    
    # Store panel for this timepoint
    p_list[[tp]] <- p
  }
  
  # Combine 3 panels into horizontal row (Pre | C1 | C2)
  return(plot_grid(plotlist = p_list, nrow = 1))
}

################################################################################
# 3. GENERATE 21-PAGE DENSITY PDF (FIGURE 8c)
################################################################################
#
# Loop through all 21 checkpoint genes and create one page per gene.
# Each page contains 3 panels (Pre, C1, C2) showing temporal evolution
# of checkpoint expression patterns across survival cohorts.
#
# PDF STRUCTURE:
# - Page 1: PDCD1 (PD-1) - Primary therapeutic target
# - Page 2: CD274 (PD-L1) - PD-1 ligand, biomarker
# - Pages 3-21: Remaining 19 checkpoint genes
# - Dimensions: 15 inches wide × 5 inches tall (accommodates 3 panels)
#
# OUTPUT: Figure 8c (Manuscript)
#

# Open PDF device for multi-page output
pdf(pdf_output, width=15, height=5)

# Generate one page per checkpoint gene (21 pages total)
for (g in unique(df_all$gene)) {
  cat("Making 2D density for gene:", g, "\n")
  
  # Call plot_2d_density() to generate 3-panel figure for this gene
  p_out <- plot_2d_density(df_all, target_gene=g)
  
  # Print to PDF (adds new page)
  print(p_out)
}

# Close PDF device
dev.off()

cat("Saved 2D Density Plots to:", pdf_output, "\n")

#########################################################################################
# END OF SCRIPT
#########################################################################################
# SUMMARY OF OUTPUTS:
# 1. FeaturePlot PDF (21 rows × 9 columns) - Comprehensive expression overview
# 2. FeaturePlot data RDS - Extracted coordinates and expression for density analysis
# 3. 2D Density PDF (Figure 8c) - 21 pages showing cohort-specific spatial patterns
#
# KEY FINDINGS VISUALIZED:
# - Baseline checkpoint expression differences between cohorts
# - Treatment-induced changes in checkpoint expression (Pre → C1 → C2)
# - Spatial localization of checkpoint-expressing T cell subsets
# - Cohort-specific patterns predictive of therapeutic response
#########################################################################################


