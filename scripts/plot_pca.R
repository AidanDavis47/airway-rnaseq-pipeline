##Script to make the pca plot


## route any errors to snakemake
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

##get rid of start up messages
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

## get the DESeq2 object from the snakemake rules
dds <- readRDS(snakemake@input[["dds"]])

## transform the counts so that most expressed gene is not dominant
if (nrow(dds) >= 1000){
    vsd <- vst(dds, blind = TRUE)
} else {
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
}

## get the PCA coordinates also get variance
pca <- plotPCA(vsd, intgroup = c("treatment", "donor"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))


## creat the actual plot
## treatment is identified via color and donor is identified via shape

plot <- ggplot(pca, aes(PC1, PC2, color = treatment, shape = donor)) + 
geom_point(size = 4) +
xlab(paste0("PC1 (", pv[1], "% of variance)")) +
ylab(paste0("PC2 (", pv[2], "% of variance)")) +
ggtitle("PCA of Airway Samples") +
theme_bw()

##save the plot
ggsave(snakemake@output[["pca"]], plot, width = 6, height = 4.5, dpi = 150)