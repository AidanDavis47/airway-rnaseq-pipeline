##script to make the heatmap plot


## route any errors to snakemake (pulled from pca plot)
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

##get rid of start up messages (pulled from pca plot)
suppressPackageStartupMessages({
    library(DESeq2)
    library(pheatmap)
})


##get the top genes so that we can show on the heatmap via snakefile
top_genes <- snakemake@params[["top_genes"]]

##use the results tbale to get the top genes and the names associated with them via snakefile
res <- read.csv(snakemake@input[["results"]])
res <- res[!is.na(res$padj), ]

## now get the expression values from the deseq object via snakefile (pulled from the pca plot, should work the same)
dds <- readRDS(snakemake@input[["dds"]])
if (nrow(dds) >= 1000){
    vsd <- vst(dds, blind = TRUE)
} else {
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
}

## get the top genes from the top row of results
top <- head(res$gene, top_genes)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat) ##make it so that color shows if a gene has gone up or down
rownames(mat) <- make.unique(res$symbol[match(top, res$gene)])

## have colors of bars be based of treatment and donor
annotation <- as.data.frame(colData(vsd)[, c("treatment", "donor")])

## now make the actual heat map
pheatmap(mat,
annotation_col = annotation,
show_colnames = FALSE,
main = paste("Top", top_genes, "Differentially Expressed Genes"),
filename = snakemake@output[["heatmap"]],
width = 7, height = 8)