# RNA-seq differential expression analysis
# GSE212353: Primary human CD8+ T cells
# Comparison: 24 h CD3/CD28 stimulation vs 0 h

library(tximport)
library(DESeq2)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(pheatmap)
# Sample metadata
samples <- c("SRR21354586", "SRR21354585",
             "SRR21354584", "SRR21354583")

condition <- factor(c("0h", "0h", "24h", "24h"),
                    levels=c("0h", "24h"))

donor <- factor(c("D1", "D2", "D1", "D2"))

coldata <- data.frame(
  row.names=samples,
  donor=donor,
  condition=condition
)
# Salmon quantification files
files <- file.path("../salmon_quant", samples, "quant.sf")
names(files) <- samples

# Check that all quantification files exist
stopifnot(all(file.exists(files)))

# Import GENCODE annotation and create transcript-to-gene mapping
gtf <- import("../reference/gencode.v50.annotation.gtf.gz")

tx2gene <- unique(
  as.data.frame(mcols(gtf[gtf$type == "transcript"]))[, c("transcript_id", "gene_id")]
)

# Import Salmon transcript-level estimates and summarize to genes
txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  ignoreAfterBar=TRUE
)

colnames(txi$counts) <- samples
# Differential expression analysis
# Paired design accounts for donor-to-donor variability
dds <- DESeqDataSetFromTximport(
  txi,
  colData=coldata,
  design=~ donor + condition
)

# Keep genes with at least 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# 24 h stimulation versus 0 h
res <- results(
  dds,
  contrast=c("condition", "24h", "0h")
)
# Add gene symbols to DESeq2 results
gene_info <- unique(
  as.data.frame(mcols(gtf[gtf$type == "gene"]))[, c("gene_id", "gene_name")]
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- merge(
  res_df,
  gene_info,
  by="gene_id",
  all.x=TRUE,
  sort=FALSE
)

# Save differential expression results
write.csv(
  res_df,
  "../results/DESeq2_24h_vs_0h_results.csv",
  row.names=FALSE
)
# PCA
vsd <- vst(dds, blind=FALSE)

pca_plot <- plotPCA(
  vsd,
  intgroup=c("condition", "donor")
)

pca_plot <- pca_plot +
  labs(
    title="CD3/CD28 Stimulation: 24 h vs 0 h"
  )

ggsave(
  "../results/PCA_0h_vs_24h.png",
  plot=pca_plot,
  width=7,
  height=5,
  dpi=300
)
# Volcano plot
res_df$significance <- "Not significant"

res_df$significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange > 1
] <- "Up"

res_df$significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange < -1
] <- "Down"

volcano_plot <- ggplot(
  res_df,
  aes(
    x=log2FoldChange,
    y=-log10(padj),
    color=significance
  )
) +
  geom_point(alpha=0.5, size=1.2) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  labs(
    title="CD3/CD28 Stimulation: 24 h vs 0 h",
    x="log2 Fold Change",
    y="-log10 adjusted p-value"
  )
label_genes <- c("IL2RA", "MIR155HG", "SLC7A5", "FABP5", "NAMPT")

volcano_plot <- volcano_plot +
  geom_text_repel(
    data=res_df[res_df$gene_name %in% label_genes, ],
    aes(label=gene_name),
    color="black",
    size=3,
    vjust=-0.8,
    show.legend=FALSE
  )
ggsave(
  "../results/Volcano_0h_vs_24h.png",
  plot=volcano_plot,
  width=7,
  height=6,
  dpi=300
)

# Heatmap of top 30 differentially expressed genes
top_genes <- head(
  res_df$gene_id[order(res_df$padj, na.last=NA)],
  30
)

heatmap_mat <- assay(vsd)[top_genes, ]

# Use gene symbols as row names
gene_labels <- gene_info$gene_name[
  match(rownames(heatmap_mat), gene_info$gene_id)
]

rownames(heatmap_mat) <- gene_labels

# Sample annotation
annotation_col <- data.frame(
  Condition=condition,
  Donor=donor
)

rownames(annotation_col) <- samples

pheatmap(
  heatmap_mat,
  scale="row",
  annotation_col=annotation_col,
  show_colnames=TRUE,
  fontsize_row=8,
  filename="../results/Heatmap_top30_DEGs.png",
  width=7,
  height=9
)

