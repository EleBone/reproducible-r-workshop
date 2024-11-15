# DGE analysis workflow

# Optionally have an argument that removes the results directory
# for a clean slate run, using
unlink("results", recursive = TRUE)

source("src/helpers.R")
source("src/load_data.R")
source("src/data_prep.R")
source("src/dge_analysis.R")
source("src/plots.R")
rmarkdown::render("src/report.Rmd", output_dir = save_dir)
