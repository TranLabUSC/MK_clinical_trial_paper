library(Seurat)
library(ggplot2)
library(dplyr)


MK_NC_Monocyte_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_seurat_obj.RDS")
MK_NC_Monocyte_Cells_seurat_obj_control <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(1, 4, 8, 9, 11))
MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint, levels = c("Pre", "C1", "C2"))
MK_NC_Monocyte_Cells_seurat_obj_exp <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21))
MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint, levels = c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36"))


MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("Pre", "C1"), ])
MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("Pre", "C1"), ])

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control, cells = MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1_cells, group.by = "TimePoint", label = TRUE)
DimPlot(MK_NC_Monocyte_Cells_seurat_obj_exp, cells = MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1_cells, group.by = "TimePoint", label = TRUE)

MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("C1"), ])
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells))
control_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells, min_num_cells)
control_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells, min_num_cells)
control_Pre_C1_cells <- c(control_Pre_sampled, control_C1_sampled)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control, cells = control_Pre_C1_cells, group.by = "TimePoint", label = TRUE)


MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("C1"), ])
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells))
exp_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells, min_num_cells)
exp_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells, min_num_cells)
exp_Pre_C1_cells <- c(exp_Pre_sampled, exp_C1_sampled)

MK_NC_Monocyte_Cells_seurat_obj_control_Pre <- subset(MK_NC_Monocyte_Cells_seurat_obj_control, cells = MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control_Pre, label = TRUE, split.by = "Patient")

MK_NC_Monocyte_Cells_seurat_obj_control_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_control, cells = MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control_C1, label = TRUE, split.by = "Patient")


MK_NC_Monocyte_Cells_seurat_obj_exp_Pre <- subset(MK_NC_Monocyte_Cells_seurat_obj_exp, cells = MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre, label = TRUE, split.by = "Patient")

MK_NC_Monocyte_Cells_seurat_obj_exp_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_exp, cells = MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_exp_C1, label = TRUE, split.by = "Patient")


############################################################################################################################################################################
MK_NC_Monocyte_Cells_seurat_obj <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_NC_Monocyte_Cells_res_1_seurat_obj.RDS")
MK_NC_Monocyte_Cells_seurat_obj_control <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(1, 4, 8, 9, 11))
MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint, levels = c("Pre", "C1", "C2"))
MK_NC_Monocyte_Cells_seurat_obj_exp <- subset(MK_NC_Monocyte_Cells_seurat_obj, subset = Patient %in% c(2, 3, 5, 7, 10, 12, 13, 14, 18, 19, 20, 21))
MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint <- factor(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint, levels = c("Pre", "C1", "C2", "C4", "C6", "C9", "C18", "C36"))


MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_control@meta.data[MK_NC_Monocyte_Cells_seurat_obj_control@meta.data$TimePoint %in% c("C1"), ])
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells))
control_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_cells, min_num_cells)
control_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_control_C1_cells, min_num_cells)
control_Pre_C1_cells <- c(control_Pre_sampled, control_C1_sampled)

MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_control, cells = control_Pre_C1_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_control_Pre_C1, group.by = "TimePoint", label = TRUE, split.by = "Patient")


MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("Pre"), ])
MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data[MK_NC_Monocyte_Cells_seurat_obj_exp@meta.data$TimePoint %in% c("C1"), ])
min_num_cells <- min(length(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells), length(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells))
exp_Pre_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_cells, min_num_cells)
exp_C1_sampled <- sample(MK_NC_Monocyte_Cells_seurat_obj_exp_C1_cells, min_num_cells)
exp_Pre_C1_cells <- c(exp_Pre_sampled, exp_C1_sampled)

MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1 <- subset(MK_NC_Monocyte_Cells_seurat_obj_exp, cells = exp_Pre_C1_cells)

DimPlot(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, group.by = "TimePoint", label = TRUE, split.by = "Patient")

grep("IL*", MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1@assays$RNA@counts@Dimnames[[1]])

FeaturePlot(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, features = c("IL6", "IL12A"), split.by = "TimePoint", label = TRUE, cols = c("lightblue", "red"), raster = FALSE)


# Extract gene names using rownames()
gene_names <- rownames(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1)

# Filter and display IL gene names
IL_genes <- gene_names[grep("^IL", gene_names)]
print(IL_genes)


# Step 1: Extract gene names
gene_names <- MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1@assays$RNA@counts@Dimnames[[1]]

# Step 2: Filter gene names starting with "IL"
IL_genes <- gene_names[grep("^IL", gene_names)]

# Step 3: View the results
print(IL_genes)


p <- FeaturePlot(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, features = IL_genes, split.by = "TimePoint", label = TRUE, cols = c("lightblue", "red"), raster = FALSE)
library(patchwork)


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
# Step 1: Create the initial FeaturePlot
p <- FeaturePlot(
  MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1,
  features = chemotaxis_genes,
  split.by = "TimePoint",
  label = TRUE,
  cols = c("lightblue", "red"),
  raster = FALSE,
  pt.size = 0.2,
  min.cutoff = 'q10',
  max.cutoff = 'q90'
)

# # Step 2: Define the function to adjust transparency
# adjust_transparency <- function(plot) {
#   plot$layers[[1]]$aes_params$alpha <- 0.5  # Set transparency to 0.5
#   plot + guides(alpha = FALSE)              # Remove alpha legend
# }
# 
# # Step 3: Apply the function to the plot(s)
# if (is.list(p)) {
#   p <- lapply(p, adjust_transparency)
# } else {
#   p <- adjust_transparency(p)
# }
# 
# # Combine the list of plots into a single plot layout
# if (is.list(p)) {
#   combined_plot <- wrap_plots(p)
# } else {
#   combined_plot <- p
# }
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/NC_Monocyte_res_1_Marker_Plots"
pdf_filename <- file.path(output_dir, paste0("chemotaxis_genes.pdf"))
ggsave(pdf_filename, plot = p, device = "pdf", width = 10, height = 1200, limitsize = FALSE)


# Use scaled data for expression values
p <- FeaturePlot(
  MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1,
  features = c("IL22RA1"),
  split.by = "TimePoint",
  label = TRUE,
  cols = c("lightblue", "red"),
  raster = FALSE,
  pt.size = 0.5,
  min.cutoff = -2,      # Adjust based on scaled data range
  max.cutoff = 2,
  slot = "scale.data"   # Use scaled data
)




# Calculate exact cutoff values (replace 'gene_of_interest' with an actual gene)
expression_values <- FetchData(MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1, vars = IL_genes)
min_value <- 0          # Set to desired minimum expression value
max_value <- 5          # Set to desired maximum expression value

# Adjusted FeaturePlot with numerical cutoffs
p <- FeaturePlot(
  MK_NC_Monocyte_Cells_seurat_obj_exp_Pre_C1,
  features = c("IL22RA1"),
  split.by = "TimePoint",
  label = TRUE,
  cols = c("lightblue", "red"),
  raster = FALSE,
  pt.size = 0.2,
  min.cutoff = min_value,
  max.cutoff = max_value,
)



######################################################################################
pre_subset <- subset(seurat_obj, subset = TimePoint == "Pre")
c1_subset <- subset(seurat_obj, subset = TimePoint == "C1")