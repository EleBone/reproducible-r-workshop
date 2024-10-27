# DGE analysis workflow

# Optionally have an argument that removes the results directory
# for a clean slate run, using
# unlink("results", recursive = TRUE)

# Output results/save directory
save_dir <- "results/"
dir.create(save_dir)

source("src/load_data.R")
source("src/data_prep.R")
source("src/dge_analysis.R")
source("src/plots.R")
