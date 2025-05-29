library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(broom)
library(ggfortify)
library("survminer")
library("Rcpp")
library(cowplot)
library(tidyr)
library(rlang)


# subsetting for CD8 T cells
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS")
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = seurat_clusters %in% c(1, 2, 3, 8, 10, 12, 14, 16, 17))
# only considering UF data as WUSTL TCR data is not good quality
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = Site == "UF")
seurat_metadata <- seurat_object_t_cells@meta.data
clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
clone_cell_barcode_df <- read.csv(clone_cell_barcode_file)
clone_cell_barcode_df <- clone_cell_barcode_df[clone_cell_barcode_df$chain == "TRB",]
clone_cell_barcode_df <- clone_cell_barcode_df[, c("barcode", "cdr3")]


# Run the table command on the specified column
value_freq_table <- table(clone_cell_barcode_df$cdr3)
# Convert the table output to a dataframe
value_freq_df <- as.data.frame(value_freq_table)
# Rename the columns for clarity
colnames(value_freq_df) <- c("cdr3", "Frequency")

clone_cell_barcode_df <- merge(clone_cell_barcode_df, value_freq_df, by = "cdr3")
# Filter the dataframe to keep the row with higher frequency for duplicated barcodes
clone_cell_barcode_df <- clone_cell_barcode_df %>%
  group_by(barcode) %>%
  arrange(desc(Frequency)) %>%
  slice(1) %>%
  ungroup()

seurat_object_tcr_cells <- subset(seurat_object_t_cells, cells = clone_cell_barcode_df$barcode)

seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, clone_cell_barcode_df, by.x = "row.names", by.y = "barcode")
rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells@meta.data

seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
  group_by(origin, cdr3) %>%
  mutate(clone_frequency = n()) %>%
  ungroup()

# Calculate the total number of cells in each origin
total_cells_in_origin <- seurat_object_tcr_cells_metadata %>%
  group_by(origin) %>%
  summarise(total_cells = n())

# Join the total cell count back to the metadata
seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
  left_join(total_cells_in_origin, by = "origin")

# Calculate clone proportion
seurat_object_tcr_cells_metadata <- seurat_object_tcr_cells_metadata %>%
  mutate(clone_proportion = clone_frequency / total_cells) %>%
  select(-total_cells)  # Optionally remove the total_cells column

seurat_object_tcr_cells@meta.data <- merge(seurat_object_tcr_cells@meta.data, seurat_object_tcr_cells_metadata[, c("Row.names", "clone_frequency", "clone_proportion")], by.x = "row.names", by.y = "Row.names")
rownames(seurat_object_tcr_cells@meta.data) <- seurat_object_tcr_cells@meta.data$Row.names
seurat_object_tcr_cells_metadata <- as.data.frame(seurat_object_tcr_cells_metadata)
rownames(seurat_object_tcr_cells_metadata) <- seurat_object_tcr_cells_metadata$Row.names

##########################################################################################################################################################
# subset seurat object for cell barcodes where barcodes match the TCR barcodes and tag cell barcodes as "shared" or "not shared"
# map TCR barcodes
blood_clone_cell_barcode_file <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/vdj_t/filtered_contig_annotations.csv"
blood_clone_cell_barcode_df <- read.csv(blood_clone_cell_barcode_file)
blood_clone_cell_barcode_df <- blood_clone_cell_barcode_df[blood_clone_cell_barcode_df$chain == "TRB",]
cell_barcodes <- unique(blood_clone_cell_barcode_df$barcode)

seurat_object_tcr_barcodes <- subset(seurat_object_t_cells, cells = cell_barcodes)

# add an attribute to the seurat object metadata, "clone_status"
# Initialize the clone_status column in metadata
seurat_object_tcr_barcodes@meta.data$clone_proportion_status <- "non_dominant"

# Get the list of donors
donors <- unique(seurat_object_tcr_barcodes@meta.data$donor)

# 0.01 is the threshold for our frequent clones
high_proportion_clone_cells <- rownames(seurat_object_tcr_cells_metadata[seurat_object_tcr_cells_metadata$clone_proportion > 0.01,])

# Iterate over each donor
for (donor in donors) {
  # Get the cell barcodes for the current donor
  donor_cells <- rownames(seurat_object_tcr_barcodes@meta.data[seurat_object_tcr_barcodes@meta.data$donor == donor, ])
  patient_id <- donor
  # Find the shared cell barcodes for the current donor
  hp_cells <- intersect(donor_cells, high_proportion_clone_cells)
  
  # Update the clone_status for shared cells
  seurat_object_tcr_barcodes@meta.data[hp_cells, "clone_proportion_status"] <- "dominant"
}




####################################################################################################################
# neoantigen pathway activity
# ────────────────────────────────────────────────────────────────
# 0.  Libraries  -------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
})

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





# #####################################################################################################################
# # two cohort version - control and experiment
# # ────────────────────────────────────────────────────────────────
# # 1.  Inputs  ----------------------------------------------------
# 
# ## 1b) Survival / clinical table
# survival_data_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
# 
# ## 1c) Pathways & GMT
# gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
# pathways <- c(
#   "VanderLeun_2020",
#   "Caushi_2021",
#   "Lowery_2022",
#   "Hanada_2022​",
#   "Oliveira_2021​",
#   "Combined_Neoantigen"
# )
# 
# ## 1d) Output directory
# out_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Neoantigen_Analysis/"
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# 
# ## 1e) Choose “ratio” or “difference”
# method_choice <- "difference"   # or "difference"
# 
# # ────────────────────────────────────────────────────────────────
# # 2.  Split object by clone status  ------------------------------
# dominant_obj      <- subset(seurat_object_tcr_barcodes,
#                             subset = clone_proportion_status == "dominant")
# non_dominant_obj  <- subset(seurat_object_tcr_barcodes,
#                             subset = clone_proportion_status == "non_dominant")
# 
# # ────────────────────────────────────────────────────────────────
# # 3.  Helper: build long‑format df for one subset & pathway  -----
# make_long_df <- function(seurat_subset, clone_status_label, pathway) {
#   
#   df <- create_survival_data(
#     gmt_file     = gmt_file,
#     pathway_name = pathway,
#     seurat_obj   = seurat_subset,
#     survival_data= survival_data_path
#   )
#   
#   ## keep UF site & non‑IDH‑POS, mirroring your earlier filters
#   df <- df %>%
#     dplyr::filter(site == "UF", IDH != "POS") %>%
#     dplyr::select(patient_id, Arm,
#                   Mean_Expr_C1, Mean_Expr_C2) %>%
#     dplyr::filter(!is.na(Mean_Expr_C1), !is.na(Mean_Expr_C2)) %>%
#     dplyr::mutate(
#       signal_change = if (method_choice == "ratio")
#         Mean_Expr_C2 / Mean_Expr_C1
#       else
#         Mean_Expr_C2 - Mean_Expr_C1,
#       clone_status  = clone_status_label
#     )
#   
#   return(df)
# }
# 
# # ────────────────────────────────────────────────────────────────
# # 4.  Loop through pathways  -------------------------------------
# for (pw in pathways) {
#   message("Processing pathway: ", pw)
#   
#   df_dom  <- make_long_df(dominant_obj,     "dominant",      pw)
#   df_nond <- make_long_df(non_dominant_obj, "non_dominant",  pw)
#   
#   plot_df <- bind_rows(df_dom, df_nond)
#   
#   ## ── 4a. ensure desired ordering ─────────────────────────────
#   plot_df$Arm          <- factor(plot_df$Arm,
#                                  levels = c("MK-3475 Alone",
#                                             "MK-3475 + MLA"))
#   plot_df$clone_status <- factor(plot_df$clone_status,
#                                  levels = c("non_dominant", "dominant"))
#   
#   ## ── 4b. colours ─────────────────────────────────────────────
#   box_cols <- c(non_dominant = "blue",
#                 dominant      = "red")
#   
#   ## ── 4c. Wilcoxon p‑values per cohort ────────────────────────
#   pvals <- plot_df %>%
#     group_by(Arm) %>%
#     summarise(p = tryCatch(
#       wilcox.test(signal_change ~ clone_status)$p.value,
#       error = function(e) NA_real_),
#       y = max(signal_change, na.rm = TRUE) * 1.05,
#       .groups = "drop")
#   
#   ## ── 4d. BOX‑&‑JITTER plot with p‑values ---------------------
#   p_box <- ggplot(plot_df,
#                   aes(x    = clone_status,
#                       y    = signal_change,
#                       fill = clone_status)) +
#     geom_boxplot(width = 0.5,
#                  position = position_dodge(width = 1.2),
#                  outlier.shape = NA) +
#     geom_jitter(position = position_jitterdodge(jitter.width = 1,
#                                                 dodge.width  = 1.2),
#                 size = 4, alpha = 0.8) +
#     facet_wrap(~ Arm, nrow = 1) +
#     scale_fill_manual(values = box_cols) +
#     labs(title = paste0(pw, "  |  C1 → C2"),
#          x     = NULL,
#          y     = ifelse(method_choice == "ratio",
#                         "Fold‑change (C2 / C1)",
#                         "Mean difference (C2 – C1)")) +
#     theme_minimal(base_size = 14) +
#     theme(axis.line        = element_line(color = "black"),
#           axis.title.x     = element_blank(),
#           axis.title.y     = element_text(size = 20),
#           plot.title       = element_text(size = 20),
#           axis.text.y      = element_text(size = 15),
#           axis.text.x      = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.ticks.y     = element_line(),
#           legend.position  = "none") +
#     geom_text(data = pvals,
#               aes(x = 1.5, y = y,
#                   label = paste0("p = ",
#                                  formatC(p, digits = 3, format = "g"))),
#               inherit.aes = FALSE,
#               size = 4)
#   
#   ## ── 4e. LONG data for dot‑plots -----------------------------
#   dot_df <- plot_df %>%
#     pivot_longer(cols = starts_with("Mean_Expr_"),
#                  names_to  = "TimePoint",
#                  values_to = "Expr") %>%
#     mutate(TimePoint   = recode(TimePoint,
#                                 "Mean_Expr_C1" = "C1",
#                                 "Mean_Expr_C2" = "C2"),
#            TimePoint   = factor(TimePoint,
#                                 levels = c("C1", "C2")))
#   
#   ## ── 4f. DOT‑PLOTS (facet grid clone_status × Arm) -----------
#   p_dot <- ggplot(dot_df,
#                   aes(x = TimePoint,
#                       y = Expr,
#                       group = patient_id)) +
#     geom_line(aes(color = clone_status), size = 1, alpha = 0.6) +
#     geom_point(aes(color = clone_status), size = 3) +
#     facet_grid(clone_status ~ Arm) +
#     scale_color_manual(values = box_cols, guide = "none") +
#     labs(title = paste0(pw, "  |  Per‑patient pathway activity"),
#          x = NULL,
#          y = "Mean pathway expression") +
#     theme_minimal(base_size = 14) +
#     theme(axis.line        = element_line(color = "black"),
#           strip.text       = element_text(size = 12, face = "bold"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
#   
#   ## ── 4g. Combine panels & save -------------------------------
#   combined <- cowplot::plot_grid(p_box, p_dot,
#                                  ncol = 1,
#                                  rel_heights = c(1, 1.2))
#   
#   ggsave(filename = file.path(out_dir,
#                               paste0("CloneDominance_", pw,
#                                      "_C1vsC2_", method_choice, ".pdf")),
#          plot     = combined,
#          device   = "pdf",
#          width    = 8, height = 10)
# }



















#####################################################################################################################
# three cohort version - control, short-term and long-term
# ────────────────────────────────────────────────────────────────
# 1.  Inputs  ----------------------------------------------------

# ── NEW: define the three cohorts ───────────────────────────────
control_group            <- c(1, 4, 8, 9, 11)
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

## 1b) Survival / clinical table
survival_data_path <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"

## 1c) Pathways & GMT
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
pathways <- c(
  "VanderLeun_2020",
  "Caushi_2021",
  "Lowery_2022",
  "Hanada_2022​",
  "Oliveira_2021​",
  "Combined_Neoantigen"
)

## 1d) Output directory
out_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/Finale/Neoantigen_Analysis/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 1e) Choose “ratio” or “difference”
method_choice <- "ratio"   # or "difference"

# ────────────────────────────────────────────────────────────────
# 2.  Split object by clone status  ------------------------------
dominant_obj      <- subset(seurat_object_tcr_barcodes,
                            subset = clone_proportion_status == "dominant")
non_dominant_obj  <- subset(seurat_object_tcr_barcodes,
                            subset = clone_proportion_status == "non_dominant")

# ────────────────────────────────────────────────────────────────
# 3.  Helper: build long‑format df for one subset & pathway  -----
make_long_df <- function(seurat_subset, clone_status_label, pathway) {
  
  df <- create_survival_data(
    gmt_file     = gmt_file,
    pathway_name = pathway,
    seurat_obj   = seurat_subset,
    survival_data= survival_data_path
  )
  
  ## keep UF + non‑IDH‑POS and add three‑way cohort label
  df <- df %>%
    dplyr::filter(site == "UF", IDH != "POS") %>%
    dplyr::mutate(
      SurvivalGroup = dplyr::case_when(
        patient_id %in% short_term_survivor_group ~ "ShortTerm",
        patient_id %in% long_term_survivor_group  ~ "LongTerm",
        patient_id %in% control_group             ~ "Control",
        TRUE                                      ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(SurvivalGroup)) %>%
    dplyr::select(patient_id, SurvivalGroup,
                  Mean_Expr_C1, Mean_Expr_C2) %>%
    dplyr::filter(!is.na(Mean_Expr_C1), !is.na(Mean_Expr_C2)) %>%
    dplyr::mutate(
      signal_change = if (method_choice == "ratio")
        Mean_Expr_C2 / Mean_Expr_C1
      else
        Mean_Expr_C2 - Mean_Expr_C1,
      clone_status  = clone_status_label
    )
  
  return(df)
}

# ────────────────────────────────────────────────────────────────
# 4.  Loop through pathways  -------------------------------------
for (pw in pathways) {
  message("Processing pathway: ", pw)
  
  df_dom  <- make_long_df(dominant_obj,     "dominant",      pw)
  df_nond <- make_long_df(non_dominant_obj, "non_dominant",  pw)
  
  plot_df <- dplyr::bind_rows(df_dom, df_nond)
  
  ## ordering ---------------------------------------------------------------------------------
  plot_df$SurvivalGroup <- factor(plot_df$SurvivalGroup,
                                  levels = c("Control", "ShortTerm", "LongTerm"))
  plot_df$clone_status  <- factor(plot_df$clone_status,
                                  levels = c("non_dominant", "dominant"))
  
  ## colours ----------------------------------------------------------------------------------
  box_cols <- c(non_dominant = "blue",
                dominant      = "red")
  
  ## Wilcoxon p‑values per cohort -------------------------------------------------------------
  pvals <- plot_df %>%
    group_by(SurvivalGroup) %>%
    summarise(p = tryCatch(
      wilcox.test(signal_change ~ clone_status)$p.value,
      error = function(e) NA_real_),
      y = max(signal_change, na.rm = TRUE) * 1.05,
      .groups = "drop")
  
  ## 4‑d  BOX‑AND‑JITTER ----------------------------------------------------------------------
  p_box <- ggplot(plot_df,
                  aes(x    = clone_status,
                      y    = signal_change,
                      fill = clone_status)) +
    geom_boxplot(width = 0.5,
                 position = position_dodge(width = 1.2),
                 outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 1,
                                                dodge.width  = 1.2),
                size = 4, alpha = 0.8) +
    facet_wrap(~ SurvivalGroup, nrow = 1) +
    scale_fill_manual(values = box_cols) +
    labs(title = paste0(pw, "  |  C1 → C2"),
         x     = NULL,
         y     = ifelse(method_choice == "ratio",
                        "Fold‑change (C2 / C1)",
                        "Mean difference (C2 – C1)")) +
    theme_minimal(base_size = 14) +
    theme(axis.line        = element_line(color = "black"),
          axis.title.x     = element_blank(),
          axis.title.y     = element_text(size = 20),
          plot.title       = element_text(size = 20),
          axis.text.y      = element_text(size = 15),
          axis.text.x      = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y     = element_line(),
          legend.position  = "none") +
    geom_text(data = pvals,
              aes(x = 1.5, y = y,
                  label = paste0("p = ",
                                 formatC(p, digits = 3, format = "g"))),
              inherit.aes = FALSE,
              size = 4)
  
  ## 4‑e  DOT‑PLOT grid -----------------------------------------------------------------------
  dot_df <- plot_df %>%
    tidyr::pivot_longer(cols = starts_with("Mean_Expr_"),
                        names_to  = "TimePoint",
                        values_to = "Expr") %>%
    dplyr::mutate(TimePoint = dplyr::recode(TimePoint,
                                            "Mean_Expr_C1" = "C1",
                                            "Mean_Expr_C2" = "C2"),
                  TimePoint = factor(TimePoint,
                                     levels = c("C1", "C2")))
  
  p_dot <- ggplot(dot_df,
                  aes(x     = TimePoint,
                      y     = Expr,
                      group = patient_id)) +
    geom_line(aes(color = clone_status), size = 1, alpha = 0.6) +
    geom_point(aes(color = clone_status), size = 3) +
    facet_grid(clone_status ~ SurvivalGroup) +
    scale_color_manual(values = box_cols, guide = "none") +
    labs(title = paste0(pw, "  |  Per‑patient pathway activity"),
         x = NULL,
         y = "Mean pathway expression") +
    theme_minimal(base_size = 14) +
    theme(axis.line        = element_line(color = "black"),
          strip.text       = element_text(size = 12, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ## 4‑g  assemble & save ---------------------------------------------------------------------
  combined <- cowplot::plot_grid(p_box, p_dot, ncol = 1,
                                 rel_heights = c(1, 1.2))
  
  ggsave(file.path(out_dir,
                   paste0("CloneDominance_", pw,
                          "_C1vsC2_", method_choice, "_extended.pdf")),
         combined, device = "pdf",
         width = 9, height = 11)
}