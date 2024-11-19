# DGE analysis workflow
here::i_am('main.R')

source("src/helpers.R")

# Optionally have an argument that removes the results directory
# for a clean slate run, using
unlink(save_dir, recursive = TRUE)

source("src/load_data.R")
source("src/data_prep.R")
source("src/dge_analysis.R")
source("src/plots.R")
rmarkdown::render("src/report.Rmd", output_dir = save_dir)
