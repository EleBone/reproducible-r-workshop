# Plots

# load helpers ----
source('src/helpers.R')

# load packages
library(ggplot2)
library(ggrepel)

# load data
topTable.dex <- readRDS(here(save_dir, "topTable.dex.rds"))
efit <- readRDS(here(save_dir, "efit.rds"))
filteredCountsMatUQ <- readRDS(here(save_dir, "filteredCountsMatUQ.rds"))

# Check the published DEG: CRISPLD2 (ENSG00000103196)
gene.of.interest <- "ENSG00000103196" # (i.e. CRISPLD2, if not convering from the provided Ensembl ID)
gene.of.interest <- "CRISPLD2" # after converting to gene symbols

topTable.dex[which(topTable.dex$ID==gene.of.interest),]
plot(log2(1+filteredCountsMatUQ[gene.of.interest,]), x=treatment.groups,
     ylab="UQ-normalised mRNA counts", xlab="", main=gene.of.interest)

# Could continue with typical visualisations (heatmap, volcano plot), PCA analysis, etc.

## Volcano plot ----

# basic volcano plot (x-axis log2 fold change, y-axis Significance)
volcanoplot(efit, coef = 1, style = "p-value", highlight = 10, 
            names = row.names(efit), hl.col="blue",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
# issues seen (overlapping labels, unequal x-axis range)

# custom volcano plot
volcano.plot.data.frame <- data.frame(topTable.dex, EnsembleId = row.names(topTable.dex))

n.highlighted <- 30 # number of genes to be labelled, in order of p-value
volcano.plot.data.frame$label <- NA
volcano.plot.data.frame$label[1:n.highlighted] <- 
  volcano.plot.data.frame$ID[1:n.highlighted]

# basic test version
ggplot(volcano.plot.data.frame, aes(x=logFC, y=(-log10(adj.P.Val)))) + 
  geom_point() + 
  geom_text_repel(aes(label=label), max.overlaps = Inf, size=3, colour="blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  theme_light()

# (apply customisations from manual testing)

ggplot(volcano.plot.data.frame, aes(x=logFC, y=(-log10(adj.P.Val)))) + 
  xlim(c(-8,8)) + ylim(c(0,4.5)) +
  geom_point(size=1, alpha=0.4) + 
  geom_text_repel(aes(label=label), max.overlaps = Inf, size=3, colour="blue", 
                  force = 5, min.segment.length = 0, segment.alpha = 0.2) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, colour = "darkgreen", linetype = "dashed") + 
  annotate("label", y=-log10(0.05), x=-7, label="adj.p=0.05", colour="darkgreen") +
  theme_light()

ggsave(here(save_dir, "volcano_plot.png"))
