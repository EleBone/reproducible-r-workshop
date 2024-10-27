# Data preparation ----

library(airway)

## Load and inspect experiment data ----

# The 'airway' package containes an example RNAseq experiment.
# Himes et al. 'RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.'
# PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.

data(airway, package="airway")
str(airway)
str(airway@assays@data@listData[["counts"]])

# Inspect the samples and experiment groupings (dex treatment); 4 treated vs 4 untreated
table(airway$Sample,airway$dex)
treatment.groups <- factor(airway$dex)

# Could edit here to rename the treatment groups "Control" and "DEX" to match the paper, then see below that later code needs to be updated?
  # treatment.groups <- vector()
  # treatment.groups[airway$dex=="untrt"] <- "Control"
  # treatment.groups[airway$dex=="trt"] <- "DEX"
  # treatment.groups <- factor(treatment.groups)

# Obtain the raw counts as a data matrix: genes(rows) x sample(columns)
raw.counts.matrix <- assay(airway, "counts")
str(raw.counts.matrix)
# Rows are genes, columns are samples, values are integer expression counts

# Save outputs
saveRDS(treatment.groups, file = file.path(save_dir, "treatmentGroups.rds"))
saveRDS(raw.counts.matrix, file = file.path(save_dir, "rawCountsMatrix.rds"))
