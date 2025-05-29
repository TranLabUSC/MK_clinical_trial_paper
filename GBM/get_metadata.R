library(dplyr)
sample_order =  c (
  "S",
  "C1",
  "C2",
  "C4",
  "C6",
  "C9",
  "C18",
  "C36"
)


get_metadata_sorted = function(meta_file = "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/MK/new_sequencing_data/MK2_aggregate.csv")
{
  combined_meta <- read.csv(meta_file)
  combined_meta$timepoint <- mapply(function(orig, pat) sub(paste0("^", pat), "", orig), combined_meta$origin, combined_meta$Patient)
  combined_meta$timepoint = factor(combined_meta$timepoint, levels = sample_order)
  metadata = combined_meta %>% arrange(Patient, timepoint)
  
  return(metadata)
}

to_bool = function(x)
{
  if (x == "Y")
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}