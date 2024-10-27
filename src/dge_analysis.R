# DEG analysis ----
# (Differential Gene Expression tests)

# load packages
library(limma)

# load data
treatment.groups <- readRDS(file.path(save_dir, "treatmentGroups.rds"))
raw.counts.matrix <- readRDS(file.path(save_dir, "rawCountsMatrix.rds"))
filteredCountsMat <- readRDS(file.path(save_dir, "filteredCountsMat.rds"))

## limma ----

# 'dgelist' object with groups
X.dgelist <- DGEList(filteredCountsMat, group = treatment.groups)
dim(X.dgelist)

# Normalise libraries
X.dgelist <- calcNormFactors(X.dgelist, method = "TMM")
X.dgelist$samples$norm.factors

# limma test design and contrast objects
design <- model.matrix(~ 0 + treatment.groups)
colnames(design) <- c("trt", "untrt")

contr <- makeContrasts(
  dex    = trt-untrt,
  levels = colnames(design)
)

v <- voom(counts=X.dgelist, design = design,plot = TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr)
efit <- eBayes(vfit)
plotSA(efit)

topTable.dex <- topTable(efit, coef = "dex", number = Inf, sort.by = "p")
head(topTable.dex)
hist(topTable.dex$P.Value)

saveRDS(topTable.dex, file = file.path(save_dir, "topTable.dex.rds"))
saveRDS(efit, file = file.path(save_dir, "efit.rds"))
