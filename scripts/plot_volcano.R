##script to generate a volcano plot of the results





## route any errors to snakemake (pulled from pca plot code)
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")



##get rid of start up messages (pulled and modified from pca plot code)
suppressPackageStartupMessages({
    ##dont need the deseq2 library since i am pulling directly from the results table
    library(ggplot2)
})


##get pramaters from the snakemake file for sig cutoff (padj) and fold change cutoff (lfc)
padj_cut <- snakemake@params[["padj"]]
lfc_cut <- snakemake@params[["lfc"]]

##now pull data from the results table
res <- read.csv(snakemake@input[["results"]]) ##getting the info from snake make input
res <- res[!is.na(res$padj), ] ##cut genes that DESeq2 did not or could not test

##now lable genes as not significant, went up, or went down
res$status <- "Not Significant"
res$status[res$padj < padj_cut & res$log2FoldChange > lfc_cut] <- "Up"
res$status[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Down"

#label the known steroid genes
known <- c("FKBP5", "TSC22D3", "PER1", "CRISPLD2", "DUSP1")
labels <- res[res$symbol %in% known, ]


##now creat the actual volcano plot
volcano <- ggplot(res, aes(log2FoldChange, -log10(padj), color = status)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", `Not Significant` = "grey70")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    geom_text(data = labels, aes(label = symbol), color = "black", vjust = -0.8, size = 3) +
    labs(title = "Dexamethasone Vs Untreated", x = "log2 Fold Change", y = "-log1- Adjusted P-Value", color = NULL) +
    theme_bw()

    ##save the plot to the plots folder in results via snakemake
    ggsave(snakemake@output[["volcano"]], volcano, width = 7, height = 5, dpi = 150)