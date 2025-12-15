library(optparse)
library(dplyr)
library(purrr)

make_rnk = function(infile, outfile) {
  df <- read.csv(infile, sep = "\t")
  df$rnk = df$logFC

  rnk = data.frame("#gene" = df$gene, rnk = df$rnk)
  write.table(rnk, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
}

process_files_in_directory = function(root_dir) {
  # Find all 'dea.txt' files in the directory and its subdirectories
  dea_files <- list.files(path = root_dir, pattern = "dea.txt", full.names = TRUE, recursive = TRUE)

  # Iterate over each file, process it, and save the output
  for (infile in dea_files) {
    outfile <- gsub("dea.txt", "dea.rnk", infile)
    make_rnk(infile, outfile)
  }
}

# Define options for command line arguments
option_list = list(
  make_option(c("-d", "--dir"), action="store", default=getwd(), type='character',
              help="Directory to search for dea.txt files")
)

# Parse arguments
opt = parse_args(OptionParser(option_list=option_list))

# Process files
process_files_in_directory(opt$dir)
