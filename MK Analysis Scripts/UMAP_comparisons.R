# ================================================================================
# UMAP Chemotaxis Gene Expression Comparison in Non-Classical Monocytes
# ================================================================================
#
# Purpose:
#   Compare chemotaxis gene expression in Non-Classical Monocytes between
#   timepoints Pre and C1 for experimental cohort patients to assess immune
#   cell recruitment and migration potential changes after treatment initiation.
#
# Biological Context:
#   - Non-Classical Monocytes: CD14+CD16++ monocyte subset with patrolling behavior
#   - Chemotaxis: Directed cell migration in response to chemical gradients
#   - Chemotaxis genes: Include chemokines (CCL, CXCL), receptors (CCR, CXCR),
#     and regulators of cell migration/recruitment (400+ genes total)
#   - Analysis compares Pre-treatment (baseline) vs C1 (Cycle 1 after treatment)
#
# Workflow:
#   1. Load Non-Classical Monocyte Seurat object and filter clusters
#   2. Separate into control vs experimental cohorts
#   3. Downsample to equal cell numbers (Pre vs C1) for balanced comparison
#   4. Generate UMAP visualizations for both cohorts
#   5. Filter experimental cohort cells within specific UMAP boundaries
#   6. Create comprehensive chemotaxis gene expression FeaturePlots
#   7. Generate IL (interleukin) gene expression FeaturePlots
#
# Key Parameters:
#   - Control patients: 1, 4, 8, 9
#   - Experimental patients: 2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21
#   - Retained clusters: 0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16
#   - UMAP limits: X [-7, 7], Y [-6, 8]
#   - Expression cutoffs: q10 (min) to q90 (max) quantiles
#
# Outputs:
#   - NC_Mono_Control_Pre_C1.pdf: Control cohort UMAP (Pre vs C1)
#   - NC_Mono_Experiment_Pre_C1.pdf: Experimental cohort UMAP (Pre vs C1)
#   - chemotaxis_genes.pdf: 400+ chemotaxis genes split by timepoint (Figure 5d)
#   - IL_genes_test.pdf: Interleukin genes split by timepoint
#
# Manuscript Figure:
#   Figure 5d: Chemotaxis gene expression comparison (Pre vs C1, experimental cohort)
#
# Dependencies:
#   - Seurat object: MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS
#   - Libraries: Seurat, ggplot2, dplyr, patchwork
# ================================================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# Load Non-Classical Monocyte Seurat object
MK_NC_Monocyte_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")

# ================================================================================
# SECTION 1: Data Preparation and Cohort Separation
# ================================================================================

# Remove patient-specific/artifact clusters, retain only high-quality NC monocyte clusters
MK_NC_Monocyte_Cells_seurat_obj <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16))

# Separate control cohort (patients 1, 4, 8, 9) - received no experimental treatment
MK_NC_Monocyte_Cells_seurat_obj_control <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(1, 4, 8, 9))

# Factor timepoints in chronological order for control cohort
MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint, levels = c("Pre", "C1", "C2"))

# Separate experimental cohort (patients receiving MK-3475 + MLA treatment)
MK_NC_Monocyte_Cells_seurat_obj_exp <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21))

# Factor timepoints in chronological order for experimental cohort
MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint, levels = c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36"))

# ================================================================================
# SECTION 2: Control Cohort UMAP Visualization (Pre vs C1)
# ================================================================================

# Extract Pre and C1 cells from control cohort
MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("C1"), ])

# Downsample to equal cell numbers for balanced comparison
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells))
control_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells, min_num_cells)
control_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells, min_num_cells)
control_Pre_C1_cells <- c(control_Pre_sampled, control_C1_sampled)

# Subset to downsampled Pre and C1 cells
MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_control, cells = control_Pre_C1_cells)

# Generate UMAP plot colored by timepoint (blue=Pre, red=C1)
dim_plot <- DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1, group.by = "TimePoint", cols = c("blue", "red")) +
  scale_x_continuous(limits = c(-7, 7)) +  # Standardize x-axis limits
  scale_y_continuous(limits = c(-6, 8)) +  # Standardize y-axis limits
  ggtitle("Non Classical Monocytes in Control")

# Save control cohort UMAP plot
ggsave(
  filename = "NC_Mono_Control_Pre_C1.pdf",
  plot = dim_plot,
  device = "pdf",
  path = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material",
  width = 8, height = 6,
  units = "in",
  dpi = 300
)

# ================================================================================
# SECTION 3: Experimental Cohort UMAP Visualization (Pre vs C1)
# ================================================================================

# Extract Pre and C1 cells from experimental cohort
MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("C1"), ])

# Downsample to equal cell numbers for balanced comparison
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells))
exp_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells, min_num_cells)
exp_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells, min_num_cells)
exp_Pre_C1_cells <- c(exp_Pre_sampled, exp_C1_sampled)

# Subset to downsampled Pre and C1 cells
MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_exp, cells = exp_Pre_C1_cells)

# Generate UMAP plot colored by timepoint (blue=Pre, red=C1)
dim_plot <- DimPlot(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, group.by = "TimePoint", cols = c("blue", "red")) +
  scale_x_continuous(limits = c(-7, 7)) +  # Standardize x-axis limits
  scale_y_continuous(limits = c(-6, 8)) +  # Standardize y-axis limits
  ggtitle("Non Classical Monocytes in Experiment")

# Save experimental cohort UMAP plot
ggsave(
  filename = "NC_Mono_Experiment_Pre_C1.pdf",
  plot = dim_plot,
  device = "pdf",
  path = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Publication_Material",
  width = 8, height = 6,
  units = "in",
  dpi = 300
)

# ================================================================================
# SECTION 4: Filter Experimental Cohort Cells Within UMAP Boundaries
# ================================================================================

# Extract UMAP coordinates
umap_embeddings <- Embeddings(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, reduction = "umap")

# Filter cells within standardized UMAP boundaries to remove outliers
cells_to_keep <- rownames(umap_embeddings)[
  umap_embeddings[, "UMAP_1"] >= -7 & umap_embeddings[, "UMAP_1"] <= 7 &
    umap_embeddings[, "UMAP_2"] >= -6 & umap_embeddings[, "UMAP_2"] <= 8
]

# Subset to only cells within UMAP boundaries
MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, cells = cells_to_keep)

# ================================================================================
# SECTION 5: Define Gene Lists for Expression Analysis
# ================================================================================

# Extract all gene names from RNA assay
gene_names <- MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1@assays$RNA@counts@Dimnames[[1]]

# Filter interleukin genes (genes starting with "IL")
IL_genes <- gene_names[grep("^IL", gene_names)]

library(patchwork)

# Define comprehensive chemotaxis gene list (400+ genes)
# Includes: chemokines (CCL, CXCL), receptors (CCR, CXCR), adhesion molecules,
# signaling proteins, and migration regulators
chemotaxis_genes <- c(
  "ABCC1", "ACKR2", "ACKR3", "ACKR4", "ADAM10", "ADAM17", "ADAM8",
  "ADGRE2", "AGTR1", "AIF1", "AKIRIN1", "ALOX5", "ANO6", "ANXA1",
  "ARHGEF16", "ARHGEF5", "ARRB2", "AZU1", "BCAR1", "BIN2", "BSG",
  "BST1", "C1QBP", "C3AR1", "C5", "C5AR1", "C5AR2", "CALCA",
  "CALR", "CAMK1D", "CCL1", "CCL11", "CCL13", "CCL14", "CCL15",
  "CCL16", "CCL17", "CCL18", "CCL19", "CCL2", "CCL20", "CCL21",
  "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28",
  "CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5", "CCL7", "CCL8",
  "CCN3", "CCR1", "CCR10", "CCR2", "CCR3", "CCR4", "CCR5",
  "CCR6", "CCR7", "CCR8", "CCR9", "CCRL2", "CD300H", "CD74",
  "CH25H", "CHGA", "CKLF", "CMKLR1", "CNR2", "CORO1A", "CORO1B",
  "CREB3", "CRK", "CRKL", "CSF1", "CSF1R", "CSF3R", "CTSG",
  "CX3CL1", "CX3CR1", "CXADR", "CXCL1", "CXCL10", "CXCL11",
  "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17", "CXCL2",
  "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCR1", "CXCR2",
  "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CYP19A1", "CYP7B1",
  "DAPK2", "DDT", "DEFA1", "DEFA1B", "DEFB103A", "DEFB103B",
  "DEFB104A", "DEFB104B", "DEFB109B", "DEFB110", "DEFB114",
  "DEFB124", "DEFB130A", "DEFB130B", "DEFB133", "DEFB4A",
  "DNM1L", "DOCK4", "DPEP1", "DPP4", "DUSP1", "EDN1", "EDN2",
  "EDN3", "EDNRB", "EGR3", "ELMO2", "EPHA2", "EPHB1", "F2RL1",
  "F7", "FCER1G", "FFAR2", "FGF1", "FGF16", "FGF18", "FGF2",
  "FGF4", "FGFR1", "FLT1", "FOLR2", "FPR2", "GAB1", "GAS6",
  "GBF1", "GPR15LG", "GPR18", "GPR183", "GPSM3", "GREM1",
  "HBEGF", "HGF", "HMGB1", "HMGB2", "HOXB9", "HRG", "HSD3B7",
  "HSPB1", "IL10", "IL12A", "IL16", "IL17RA", "IL17RC", "IL1B",
  "IL23A", "IL34", "IL6", "IL6R", "ITGA1", "ITGA9", "ITGB2",
  "JAM3", "JAML", "KDR", "KIT", "KLRC4-KLRK1", "KLRK1", "LBP",
  "LEF1", "LGALS3", "LGALS9", "LGMN", "LOX", "LPAR1", "LYN",
  "LYST", "MAPK1", "MAPK3", "MCU", "MDK", "MET", "MICOS10-NBL1",
  "MIF", "MIR149", "MIR15A", "MIR16-1", "MIR223", "MIR34A",
  "MIR424", "MMP2", "MMP28", "MOSPD2", "MPP1", "MSMP", "MSTN",
  "MTUS1", "NBL1", "NCKAP1L", "NEDD9", "NINJ1", "NOD2",
  "NOTCH1", "NR4A1", "NRP1", "NUP85", "OXSR1", "P2RX4", "PADI2",
  "PARVA", "PDE4B", "PDGFB", "PDGFD", "PDGFRA", "PDGFRB",
  "PERP", "PF4", "PF4V1", "PGF", "PIK3CD", "PIK3CG", "PIKFYVE",
  "PIP5K1A", "PIP5K1C", "PLA2G1B", "PLA2G7", "PLEC", "PLEKHG5",
  "PLXNB3", "PPBP", "PPIA", "PPIB", "PREX1", "PRKCD", "PRKCQ",
  "PRKD1", "PRKD2", "PRSS56", "PTK2", "PTK2B", "PTN",
  "PTPRJ", "PTPRO", "RAB13", "RAC1", "RAC2", "RAC3", "RARRES2",
  "RHOG", "RIN3", "RIPOR2", "RPL13A", "RPS19", "S100A12",
  "S100A14", "S100A7", "S100A8", "S100A9", "S1PR1", "SAA1",
  "SBDS", "SCG2", "SEMA5A", "SERPINE1", "SFTPD", "SLAMF1",
  "SLAMF8", "SLC12A2", "SLC8B1", "SLIT2", "SMOC2", "SRP54",
  "STAP1", "STK39", "SWAP70", "SYK", "TAFA4", "TGFB2",
  "THBS1", "THBS4", "TIAM1", "TIRAP", "TMEM102", "TMSB4X",
  "TNFAIP6", "TNFRSF11A", "TNFSF11", "TNFSF14", "TPBG",
  "TREM1", "TRPM2", "TRPM4", "TRPV4", "VAV1", "VAV3",
  "VCAM1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "WNK1",
  "WNT5A", "XCL1", "XCL2", "XCR1", "ZNF580"
)

# ================================================================================
# SECTION 6: Generate Chemotaxis Gene Expression FeaturePlots (Figure 5d)
# ================================================================================

# Create FeaturePlot showing all chemotaxis genes split by timepoint
p <- FeaturePlot(
  MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1,
  features = chemotaxis_genes,          # 400+ chemotaxis-related genes
  split.by = "TimePoint",                # Separate panels for Pre vs C1
  label = FALSE,
  cols = c("lightblue", "red"),          # Blue (low) to red (high) expression
  raster = FALSE,                         # Vector graphics for quality
  pt.size = 0.2,
  min.cutoff = 'q10',                     # Bottom 10% set to min color
  max.cutoff = 'q90'                      # Top 10% set to max color
)

# Save chemotaxis genes FeaturePlot (Figure 5d)
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/NC_Monocyte_res_1_Marker_Plots"
pdf_filename <- file.path(output_dir, paste0("chemotaxis_genes"))
ggsave(pdf_filename, plot = p, device = "pdf", width = 10, height = 500, limitsize = FALSE)

# ================================================================================
# SECTION 7: Generate Interleukin Gene Expression FeaturePlots
# ================================================================================

# Create FeaturePlot showing all IL genes split by timepoint
# Order=TRUE plots high-expressing cells on top for better visualization
p <- FeaturePlot(
  MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1,
  features = IL_genes,                    # All interleukin genes found in dataset
  split.by = "TimePoint",                 # Separate panels for Pre vs C1
  label = FALSE,
  cols = c("lightblue", "red"),           # Blue (low) to red (high) expression
  raster = FALSE,                          # Vector graphics for quality
  pt.size = 0.4,                           # Larger points to highlight expressing cells
  order = TRUE,                            # Plot high-expressing cells on top
  min.cutoff = 'q10',                      # Bottom 10% set to min color
  max.cutoff = 'q90'                       # Top 10% set to max color
)

# Save IL genes FeaturePlot
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/NC_Monocyte_res_1_Marker_Plots"
pdf_filename <- file.path(output_dir, paste0("IL_genes_test.pdf"))
ggsave(pdf_filename, plot = p, device = "pdf", width = 10, height = 330, limitsize = FALSE)

# ================================================================================
# End of Script
# ================================================================================

