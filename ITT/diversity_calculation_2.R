

# for each vector:
 simpson_diversity_index = function(clonal_proportions)
 {
   
   # Calculate
   #Simpson's Diversity Index (SDI) 
   #SDI = 1 - ∑(ni / N)^2
   # a vector, where each value is the proportion of each clone. Total proportion is 1
   # Gini-Simpson
   SDI = 1 - sum(clonal_proportions^2)
    return(SDI)  
 }
 
 shannon_diversity_index = function(clonal_proportions)
 {
   
   # Calculate
  
   #H = - sum(p*ln(pi))
   clonal_proportions = clonal_proportions[clonal_proportions > 0]
   h = (clonal_proportions * log(clonal_proportions))
   h = -sum(h)
   return(h)  
 }
 

get_clonotype_diversity = function(trackClonotypes_df_relative, celltype, diversity_index = "shannon", outfolder=".")
{
  diversity_dict = list(shannon = shannon_diversity_index,   simpson = simpson_diversity_index )
  
clonotype_diversity = sapply(trackClonotypes_df_relative[,2:ncol(trackClonotypes_df_relative)],diversity_dict[[diversity_index]]
  )
clonotype_diversity = as.data.frame(clonotype_diversity)
write.table(clonotype_diversity, paste0(outfolder, "/", diversity_index, "_clonal_diversity_", celltype, "_T_cells.txt"), sep="\t", row.names=TRUE, quote=FALSE)
return(clonotype_diversity)
}

run = function(celltype)
{
  setwd("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature")
  file_path = paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/clonotype_df_proportion_", celltype, "_T_cells.csv")
  trackClonotypes_df_relative <- read.csv(file_path, check.names = FALSE)
  trackClonotypes_df_relative <- trackClonotypes_df_relative[, !(colnames(trackClonotypes_df_relative) %in% c("CDR3.aa"))]
  shannon = diversity = get_clonotype_diversity(trackClonotypes_df_relative = trackClonotypes_df_relative, celltype, diversity_index = "shannon")
  simpson =diversity = get_clonotype_diversity(trackClonotypes_df_relative = trackClonotypes_df_relative, celltype, diversity_index = "simpson")
  return(list(shannon = shannon, simpson = simpson))
}

# mapping <- list(
#   "all" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18),
#   "Activated_CD4" = c(0),
#   "Effector_CD8" = c(1),
#   "Effector_Memory_Precursor_CD8" = c(2),
#   "Exhausted_T" = c(3),
#   "Gamma_Delta_T" = c(4),
#   "Active_CD4" = c(5),
#   "Naive_CD4" = c(6, 9, 18),
#   "Memory_CD4" = c(7),
#   "Stem_Like_CD8" = c(8),
#   "Effector_Memory_CD8" = c(10),
#   "Central_Memory_CD8" = c(12),
#   "GZMK_Effector_Memory_CD8" = c(13),
#   "Proliferating_Effector" = c(14, 16, 17),
#   "All_CD8" = c(1, 2, 3, 8, 10, 12, 14, 16, 17)
# )

# Define the mapping
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

# Convert the list to a data frame
celltype_to_cluster <- do.call(rbind, lapply(names(mapping), function(celltype) {
  data.frame(celltype = gsub(" ", "_", celltype), cluster = mapping[[celltype]])
}))

celltypes <- unique(celltype_to_cluster$celltype)


for (celltype in celltypes) {
  out = run(celltype)
  print(out[[1]])
  print(out[[2]])
}


