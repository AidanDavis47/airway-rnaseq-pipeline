##this is the scripts file ran by the deseq2 rule in Snakefile line 95



#sends info to snakemakes logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
})

## gets sample information (the donor and the treatment that was used)
samples <- read.csv(snakemake@input[["samples"]])
rownames(samples) <- samples$sample
samples$donor <- factor(samples$donor)
samples$treatment <- factor(samples$treatment, levels = c("untreated", "dexamethasone"))

## Salmon output files
files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

## add transcript estimates 
tx2gene <- read.csv(snakemake@input[["tx2gene"]])
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx", "gene")])

##create deseq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ donor + treatment)

## filter genes, get rid of ones that have very little reads and keep ones with reads >= 10 in 4 samples
dds <- dds[rowSums(counts(dds) >= 10) >= 4,]

#Runs test comparing dexamethasone vs untreated donors
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", "dexamethasone", "untreated"))

#now add the names of the genes, sort by p-value and save
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$symbol <- tx2gene$symbol[match(res_df$gene, tx2gene$gene)]
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, snakemake@output[["results"]], row.names = FALSE)
saveRDS(dds, snakemake@output[["dds"]])