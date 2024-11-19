# Constants used in several files
library(here)

# Output results/save directory
save_dir <- "results/"
dir.create(here(save_dir))

# Gene filtering 
min.cpm <- 1
min.cpm.fraction <- 1/4