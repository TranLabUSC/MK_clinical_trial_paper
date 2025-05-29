
get_clonal_expansion_timepoint_comparison = function( compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                      metadata,
                                                      clonal_diversity_df, outdir= "clonal_expansion_shannon_all_T_cells",
                                                      method = "division")
{
  outdir = paste0(outdir, "_", method) 
  dir.create(outdir)
  patients = unique(metadata$Patient)
  get_clonal_expansion_per_patient = function(patient)
  {
    sample1 = metadata %>% filter(Patient == patient,
                                  timepoint == compared_pair[[1]]) %>% select(sample_id)
    sample1 = sample1$sample_id

    sample2 = metadata %>% filter(Patient == patient,
                                  timepoint == compared_pair[[2]]) %>% select(sample_id)
    sample2 = sample2$sample_id

    if (length(sample1) == 0 | length(sample2) == 0)
    {
      return(NULL)
    }
    if (!(sample1 %in% rownames(clonal_diversity_df)) | !(sample2 %in% rownames(clonal_diversity_df)))
    {
      return(NULL)
    }
    clonal_diversity_sample_1  = clonal_diversity_df[sample1, "clonotype_diversity"]
    clonal_diversity_sample_2  = clonal_diversity_df[sample2, "clonotype_diversity"]
    
    if (method =="division"){
      clonal_expansion <- clonal_diversity_sample_2/clonal_diversity_sample_1
    } else {
      clonal_expansion <- clonal_diversity_sample_2 - clonal_diversity_sample_1
    }
    
    return(clonal_expansion)
  }
  
  clonal_expansion_df = data.frame(matrix(ncol = 1, nrow = length(patients)))
  rownames(clonal_expansion_df) = patients
  colnames(clonal_expansion_df) = "clonal_expansion"
  
  for (i in 1:length(patients))
  {
    patient = patients[i]
    out = get_clonal_expansion_per_patient(patient = patient)
    if (!is.null(out))
    {
      clonal_expansion_df[i, "clonal_expansion"] =  out
    }
    
  }
  
  outfile = paste0(outdir,"/",compared_pair[[1]], "_vs_", compared_pair[[2]], ".txt" )
  outfile =gsub(" ", "_", outfile)
  write.table(clonal_expansion_df, outfile, sep="\t", quote= FALSE, row.names = TRUE)
}

run_example = function(celltype, diversity_type)
{
  clonal_diversity_df_file = paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/", diversity_type, "_clonal_diversity_", celltype, "_T_cells.txt")
  outdir = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/Clonal_Expansion_Results/final_nomenclature/"
  outdir = paste0(outdir, "C1_vs_C2/") 
  dir.create(outdir)
  setwd(
    outdir
  )
  metadata = get_metadata_sorted()
  clonal_diversity_df <-
    read.delim(clonal_diversity_df_file, check.names = FALSE)
  out = get_clonal_expansion_timepoint_comparison( compared_pair = c("C1", "C2"),
                                                   metadata= metadata,
                                                   clonal_diversity_df = clonal_diversity_df,
                                                   outdir= paste0("clonal_expansion_", diversity_type, "_", celltype, "_T_cells"))
  
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


source("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/ITT/get_metadata.R")
for (celltype in celltypes) {
  for (diversity_type in c("shannon", "simpson")){
    run_example(celltype, diversity_type) 
  }
}
