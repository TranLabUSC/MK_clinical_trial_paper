# -----------------------
# 0. Libraries
# -----------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# -----------------------
# 1. Define groups
# -----------------------
short_term_survivor_group <- c(7, 10, 12, 14, 18)
long_term_survivor_group  <- c(2, 3, 5, 13, 19, 20, 21)

# Read survival data
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/Analysis/UF_WUSTL_combined/Aggregated/outs/NETZEN_analysis/new_analysis/patient_survival_data.csv"
survival_df <- read.csv(survival_data)
survival_df <- survival_df[survival_df$site == "UF", ]

# Create a column "SurvivalGroup" indicating short-term vs. long-term
survival_df$SurvivalGroup <- ifelse(
  survival_df$patient_id %in% short_term_survivor_group, "ShortTerm",
  ifelse(survival_df$patient_id %in% long_term_survivor_group,  "LongTerm", NA)
)

# OPTIONAL: Filter out patients not in either group
survival_df <- survival_df[!is.na(survival_df$SurvivalGroup), ]

# Make sure ShortTerm appears first in any factor
survival_df$SurvivalGroup <- factor(
  survival_df$SurvivalGroup,
  levels = c("ShortTerm", "LongTerm")
)

# -----------------------
# 2. Load Seurat metadata
# -----------------------
seurat_object_t_cells <- readRDS(
  "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/nextflow_pipeline/Aggregated/outs/NETZEN_analysis/clustering_res_1/MK_T_Cells_seurat_obj.RDS"
)
seurat_metadata_t_cells <- seurat_object_t_cells@meta.data

# -----------------------
# 3. Define mapping
# -----------------------
mapping <- list(
  "Exhausted_T"            = c(3)
)

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

# -----------------------
# 4. Define timepoints
# -----------------------
timepoints       <- c("Pre", "C1", "C2", "C4")
timepoint_pairs  <- list(c("Pre","C1"), c("C1","C2"), c("C2","C4"))

# -----------------------
# 5. Helper function:
#    Compute cluster proportion at one timepoint
# -----------------------
get_cluster_proportion_single_timepoint <- function(metadata,
                                                    patient_col,
                                                    timepoint_col,
                                                    cluster_col,
                                                    cluster_ids,
                                                    timepoint_to_use) {
  # Filter by the timepoint
  df_filtered <- metadata[metadata[[timepoint_col]] == timepoint_to_use, ]
  
  # For each patient, compute proportion of cells in specified cluster_ids
  proportion_df <- df_filtered %>%
    group_by(.data[[patient_col]]) %>%
    summarize(
      total_cells   = n(),
      cluster_cells = sum(.data[[cluster_col]] %in% cluster_ids)
    ) %>%
    ungroup() %>%
    mutate(
      Cluster_Proportion = cluster_cells / total_cells
    )
  
  # Rename the grouping column to "patient_id"
  colnames(proportion_df)[colnames(proportion_df) == patient_col] <- "patient_id"
  
  return(proportion_df)
}


# -----------------------
# 6. Loop over cell types
# -----------------------
output_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Tracking_Results/final_nomenclature"

for (celltype in unique(celltype_to_cluster$celltype)) {
  
  # The cluster IDs for the current cell type
  tcell_cluster_list <- as.vector(celltype_to_cluster[
    celltype_to_cluster$celltype == celltype, 
    "cluster"
  ])
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6A. Compare cluster proportion at each individual timepoint
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plot_list_proportion <- list()
  
  for (tp in timepoints) {
    # 1) Get cluster proportion at timepoint
    timepoint_df <- get_cluster_proportion_single_timepoint(
      metadata      = seurat_metadata_t_cells,
      patient_col   = "Patient",
      timepoint_col = "TimePoint",
      cluster_col   = "seurat_clusters",
      cluster_ids   = tcell_cluster_list,
      timepoint_to_use = tp
    )
    
    # 2) Merge with survival data (short vs long)
    merged_df <- merge(timepoint_df,
                       survival_df[, c("patient_id", "SurvivalGroup")],
                       by = "patient_id", all.x = TRUE)
    # Filter out rows with NA
    merged_df <- merged_df[!is.na(merged_df$SurvivalGroup), ]
    
    # 3) Wilcoxon test
    short_count <- sum(merged_df$SurvivalGroup == "ShortTerm")
    long_count  <- sum(merged_df$SurvivalGroup == "LongTerm")
    
    if (short_count >= 2 && long_count >= 2) {
      test_result <- wilcox.test(Cluster_Proportion ~ SurvivalGroup, data = merged_df)
      p_value     <- test_result$p.value
    } else {
      p_value <- NA_real_
    }
    
    # 4) Plot
    # Make sure short term is factor level 1
    merged_df$SurvivalGroup <- factor(merged_df$SurvivalGroup,
                                      levels = c("ShortTerm","LongTerm"))
    custom_colors <- c("ShortTerm" = "blue", "LongTerm" = "red")
    
    p <- ggplot(merged_df, aes(x = SurvivalGroup, y = Cluster_Proportion, fill = SurvivalGroup)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +
      geom_jitter(position = position_jitter(width = 0.1), size = 2) +
      scale_fill_manual(values = custom_colors) +
      labs(
        title = paste(celltype, "at", tp),
        x     = "Survival Group",
        y     = "Cluster Proportion"
      ) +
      theme_minimal() +
      theme(
        axis.line        = element_line(color = "black"),
        axis.title.x     = element_text(size = 16),
        axis.title.y     = element_text(size = 16),
        plot.title       = element_text(size = 18, face = "bold"),
        axis.text.y      = element_text(size = 12),
        axis.text.x      = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y     = element_line()
      ) +
      annotate(
        "text",
        x = 1.5,
        y = max(merged_df$Cluster_Proportion, na.rm = TRUE),
        label = paste("p =", signif(p_value, 4)),
        size  = 4,
        vjust = -0.5
      )
    
    plot_list_proportion[[tp]] <- p
  }
  
  # Combine the 4 timepoint plots (2x2)
  final_plot_proportion <- plot_grid(plotlist = plot_list_proportion, ncol = 2, nrow = 2, align = "v")
  
  # Save the PDF for the cluster proportion at each timepoint
  file_name_proportion <- file.path(
    output_dir,
    paste0(celltype, "_T_Cells_Cluster_Proportion_Comparison.pdf")
  )
  ggsave(filename = file_name_proportion,
         plot     = final_plot_proportion,
         device   = "pdf",
         width    = 12, 
         height   = 10,
         limitsize= FALSE)
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 6B. Compare cluster proportion CHANGE (difference & ratio)
  #     between Pre_vs_C1, C1_vs_C2, C2_vs_C4
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # We will create two separate sets of plots:
  #   1) difference
  #   2) ratio
  # and then save them to separate PDFs.
  
  # -----------------------
  # 6B-i. "difference" plots
  # -----------------------
  plot_list_difference <- list()
  
  for (tp_pair in timepoint_pairs) {
    tp1 <- tp_pair[1]
    tp2 <- tp_pair[2]
    
    # Get cluster proportion for each timepoint
    df_tp1 <- get_cluster_proportion_single_timepoint(
      metadata      = seurat_metadata_t_cells,
      patient_col   = "Patient",
      timepoint_col = "TimePoint",
      cluster_col   = "seurat_clusters",
      cluster_ids   = tcell_cluster_list,
      timepoint_to_use = tp1
    )
    df_tp2 <- get_cluster_proportion_single_timepoint(
      metadata      = seurat_metadata_t_cells,
      patient_col   = "Patient",
      timepoint_col = "TimePoint",
      cluster_col   = "seurat_clusters",
      cluster_ids   = tcell_cluster_list,
      timepoint_to_use = tp2
    )
    
    # Merge them together
    merged_df <- merge(df_tp1, df_tp2, by = "patient_id", suffixes = c(paste0("_", tp1), paste0("_", tp2)))
    merged_df <- merge(merged_df, survival_df[, c("patient_id", "SurvivalGroup")],
                       by = "patient_id", all.x = TRUE)
    
    merged_df <- merged_df[!is.na(merged_df$SurvivalGroup), ]
    
    # Compute difference
    merged_df$Difference <- merged_df[[paste0("Cluster_Proportion_", tp2)]] - 
      merged_df[[paste0("Cluster_Proportion_", tp1)]]
    
    # Wilcoxon test
    short_count <- sum(merged_df$SurvivalGroup == "ShortTerm")
    long_count  <- sum(merged_df$SurvivalGroup == "LongTerm")
    
    if (short_count >= 2 && long_count >= 2) {
      test_result <- wilcox.test(Difference ~ SurvivalGroup, data = merged_df)
      p_value     <- test_result$p.value
    } else {
      p_value <- NA_real_
    }
    
    # Plot
    merged_df$SurvivalGroup <- factor(merged_df$SurvivalGroup, levels = c("ShortTerm","LongTerm"))
    custom_colors <- c("ShortTerm" = "blue", "LongTerm" = "red")
    
    p <- ggplot(merged_df, aes(x = SurvivalGroup, y = Difference, fill = SurvivalGroup)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +
      geom_jitter(position = position_jitter(width = 0.1), size = 2) +
      scale_fill_manual(values = custom_colors) +
      labs(
        title = paste(celltype, ": Difference (", tp1, "->", tp2, ")"),
        x     = "Survival Group",
        y     = "Proportion Difference"
      ) +
      theme_minimal() +
      theme(
        axis.line        = element_line(color = "black"),
        axis.title.x     = element_text(size = 16),
        axis.title.y     = element_text(size = 16),
        plot.title       = element_text(size = 18, face = "bold"),
        axis.text.y      = element_text(size = 12),
        axis.text.x      = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y     = element_line()
      ) +
      annotate(
        "text",
        x = 1.5,
        y = ifelse(all(is.na(merged_df$Difference)), 0, max(merged_df$Difference, na.rm = TRUE)),
        label = paste("p =", signif(p_value, 4)),
        size  = 4,
        vjust = -0.5
      )
    
    plot_list_difference[[paste(tp1,tp2,sep="_vs_")]] <- p
  }
  
  # Combine difference plots in a single row or 3x1 grid
  final_plot_difference <- plot_grid(plotlist = plot_list_difference, ncol = 3)
  
  file_name_difference <- file.path(
    output_dir,
    paste0(celltype, "_T_Cells_Cluster_Proportion_Difference_Comparison.pdf")
  )
  ggsave(filename = file_name_difference,
         plot     = final_plot_difference,
         device   = "pdf",
         width    = 18, 
         height   = 6,
         limitsize= FALSE)
  
  
  # -----------------------
  # 6B-ii. "ratio" plots
  # -----------------------
  plot_list_ratio <- list()
  
  for (tp_pair in timepoint_pairs) {
    tp1 <- tp_pair[1]
    tp2 <- tp_pair[2]
    
    # Get cluster proportion for each timepoint
    df_tp1 <- get_cluster_proportion_single_timepoint(
      metadata      = seurat_metadata_t_cells,
      patient_col   = "Patient",
      timepoint_col = "TimePoint",
      cluster_col   = "seurat_clusters",
      cluster_ids   = tcell_cluster_list,
      timepoint_to_use = tp1
    )
    df_tp2 <- get_cluster_proportion_single_timepoint(
      metadata      = seurat_metadata_t_cells,
      patient_col   = "Patient",
      timepoint_col = "TimePoint",
      cluster_col   = "seurat_clusters",
      cluster_ids   = tcell_cluster_list,
      timepoint_to_use = tp2
    )
    
    # Merge them together
    merged_df <- merge(df_tp1, df_tp2, by = "patient_id", suffixes = c(paste0("_", tp1), paste0("_", tp2)))
    merged_df <- merge(merged_df, survival_df[, c("patient_id", "SurvivalGroup")],
                       by = "patient_id", all.x = TRUE)
    
    merged_df <- merged_df[!is.na(merged_df$SurvivalGroup), ]
    
    # Compute ratio
    # Avoid division by zero if proportion is 0
    eps <- 1e-9
    merged_df$Ratio <- (merged_df[[paste0("Cluster_Proportion_", tp2)]] + eps) /
      (merged_df[[paste0("Cluster_Proportion_", tp1)]] + eps)
    
    # Wilcoxon test
    short_count <- sum(merged_df$SurvivalGroup == "ShortTerm")
    long_count  <- sum(merged_df$SurvivalGroup == "LongTerm")
    
    if (short_count >= 2 && long_count >= 2) {
      test_result <- wilcox.test(Ratio ~ SurvivalGroup, data = merged_df)
      p_value     <- test_result$p.value
    } else {
      p_value <- NA_real_
    }
    
    # Plot
    merged_df$SurvivalGroup <- factor(merged_df$SurvivalGroup, levels = c("ShortTerm","LongTerm"))
    custom_colors <- c("ShortTerm" = "blue", "LongTerm" = "red")
    
    p <- ggplot(merged_df, aes(x = SurvivalGroup, y = Ratio, fill = SurvivalGroup)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +
      geom_jitter(position = position_jitter(width = 0.1), size = 2) +
      scale_fill_manual(values = custom_colors) +
      labs(
        title = paste(celltype, ": Ratio (", tp1, "->", tp2, ")"),
        x     = "Survival Group",
        y     = "Proportion Ratio"
      ) +
      theme_minimal() +
      theme(
        axis.line        = element_line(color = "black"),
        axis.title.x     = element_text(size = 16),
        axis.title.y     = element_text(size = 16),
        plot.title       = element_text(size = 18, face = "bold"),
        axis.text.y      = element_text(size = 12),
        axis.text.x      = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y     = element_line()
      ) +
      annotate(
        "text",
        x = 1.5,
        y = ifelse(all(is.na(merged_df$Ratio)), 1, max(merged_df$Ratio, na.rm = TRUE)),
        label = paste("p =", signif(p_value, 4)),
        size  = 4,
        vjust = -0.5
      )
    
    plot_list_ratio[[paste(tp1,tp2,sep="_vs_")]] <- p
  }
  
  # Combine ratio plots in a single row or 3x1 grid
  final_plot_ratio <- plot_grid(plotlist = plot_list_ratio, ncol = 3)
  
  file_name_ratio <- file.path(
    output_dir,
    paste0(celltype, "_T_Cells_Cluster_Proportion_Ratio_Comparison.pdf")
  )
  ggsave(filename = file_name_ratio,
         plot     = final_plot_ratio,
         device   = "pdf",
         width    = 18, 
         height   = 6,
         limitsize= FALSE)
  
}
