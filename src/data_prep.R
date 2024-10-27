# Data preparation

# load packages
library(biomaRt)
library(edgeR)

# Analysis variables ----

# Gene filtering
min.cpm <- 1
min.cpm.fraction <- 1/4

# load data
treatment.groups <- readRDS(file.path(save_dir, "treatmentGroups.rds"))
raw.counts.matrix <- readRDS(file.path(save_dir, "rawCountsMatrix.rds"))

# Continue analysis with loaded objects


## Filter out low-expression counts ----

# Counts per million (cpm)
cpm <- cpm(raw.counts.matrix)
dim(cpm)

# Require genes have at least 'min.cpm' in at least 'min.cpm.fraction' of the samples
keepGenes <- rowSums(cpm >= min.cpm) >= ncol(cpm) * min.cpm.fraction
sum(keepGenes)
filteredCountsMat <- data.matrix(raw.counts.matrix[keepGenes,])
dim(filteredCountsMat)
print(paste0("Keeping ", sum(keepGenes)," of ", nrow(cpm)," genes."))



# Convert Ensembl ID numbers to Gene Symbols ----
# This section covers a typical issue, where genes are referred to by an accession number (in this case Ensembl ID).
# It is easier to interpret them if they are named with their gene symbol. The biomaRt package can be used.

# get an object that references the relevant biomart database (human gene names)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# # Could specify which version of the ensembl gene names for better reproducibility?
# # Some gene names changed over time (with new research, 
# # or fixing problematic gene names eg. "MARCH7" text name autocorrecting to a 
# # date value if you ever view it in Excel... https://doi.org/10.1371/journal.pcbi.1008984)
# ensembl_110 <- useEnsembl(biomart = 'genes', 
#                           dataset = 'hsapiens_gene_ensembl',
#                           version = 110)

ensemble.ids=row.names(filteredCountsMat)
listAttributes(mart = ensembl)
gene.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
                    filters="ensembl_gene_id", values = ensemble.ids,
                    mart = ensembl, uniqueRows = TRUE)

# an easy mistake would be to mismatch the rows, or not check for 1:1 relationship. (manual checking steps)
dim(gene.table)
table(is.na(gene.table$hgnc_symbol)) 
table(duplicated(gene.table$hgnc_symbol)) 

# inspect some of the problem rows
which(duplicated(gene.table$hgnc_symbol))
gene.table[2809,]
gene.table$hgnc_symbol[which(duplicated(gene.table$hgnc_symbol))]

# There is a blank value for gene symbol in a few hundred cases.
# They can fall back to Ensemble Id.
gene.table$hgnc_symbol[gene.table$hgnc_symbol==""] <- gene.table$ensembl_gene_id[gene.table$hgnc_symbol==""]
table(is.na(gene.table$hgnc_symbol)) #ok
table(duplicated(gene.table$hgnc_symbol)) #ok

gene.table$hgnc_symbol[grep(x=gene.table$hgnc_symbol,"MARCH")] 
# just showing at some point they inserted an "F" in the MARCH* genes

# Using the gene.table to convert gene names
row.names(filteredCountsMat) <- gene.table$hgnc_symbol[match(
  row.names(filteredCountsMat), table = gene.table$ensembl_gene_id)]

## Normalise counts ----

# Upper-quartile normalise function
# This is a method to correct or differences in library depth between samples.
# The library-normalised counts can then be used to compare expression of a gene between the samples.
quartileNormalise <- function(data, q=4) {
  # Get the upper quartile (if q=4) value of each sample (samples as the rows of 'data')
  upperquartiles <- vector(length = nrow(data))
  for (i in 1:nrow(data)) {
    upperquartiles[i] <- quantile(data[i,])[q]
  }
  # now scale the data.. divide by uq and multiply by the average-uq
  normalised.data <- data
  for (i in 1:nrow(data)) {
    normalised.data[i,] <- normalised.data[i,] / upperquartiles[i] * median(upperquartiles)
  }
  # give back the normalised data result
  return(normalised.data)
}
filteredCountsMatUQ <- t(quartileNormalise(t(filteredCountsMat)))

# View distribution before normalising
boxplot(log2(1+filteredCountsMat), las=2, cex.axis=0.5, cex=0.5,
        ylab="log2(1+counts)")

# View distribution after normalising
boxplot(log2(1+filteredCountsMatUQ), las=2, cex.axis=0.5, cex=0.5,
        ylab="log2(1+counts)")

# Exploratory plots to view the amount of data collected in each sample:
# Plot the depth of sequencing per sample
plot(colSums(filteredCountsMat)/10^6, ylab="Counts (millions)", 
     ylim=c(0, max(colSums(filteredCountsMat)/10^6)),
     xlab = "Sample #")
# Add experiment groupings
plot(y = colSums(filteredCountsMat)/10^6, ylab="Counts (millions)", 
     ylim=c(0, max(colSums(filteredCountsMat)/10^6)),
     x = treatment.groups)

# Output: table of prepared gene expression data (counts x samples)
write.csv(filteredCountsMatUQ, file = paste0(save_dir, "filteredCountsMatUQ.csv"))

saveRDS(filteredCountsMatUQ, file = file.path(save_dir, "filteredCountsMatUQ.rds"))
saveRDS(filteredCountsMat, file = file.path(save_dir, "filteredCountsMat.rds"))
saveRDS(gene.table, file = file.path(save_dir, "geneTable.rds"))
