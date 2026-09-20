# GSEA using GO Biological Process gene sets
# 24 h CD3/CD28 stimulation vs 0 h

library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(ggplot2)

# Load DESeq2 results
res_df <- read.csv(
  "../results/DESeq2_24h_vs_0h_results.csv",
  stringsAsFactors=FALSE
)

# Remove genes without test statistics
res_df <- res_df[
  !is.na(res_df$stat) &
  !is.na(res_df$gene_name) &
  res_df$gene_name != "",
]

# Create ranked gene list
# Positive values = higher expression at 24 h
# Negative values = lower expression at 24 h
ranks <- res_df$stat
names(ranks) <- res_df$gene_name

# Remove duplicated gene symbols
ranks <- ranks[!duplicated(names(ranks))]

# Sort from most positive to most negative
ranks <- sort(ranks, decreasing=TRUE)

# GO Biological Process annotations
go_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=names(ranks),
  keytype="SYMBOL",
  columns=c("SYMBOL", "GO", "ONTOLOGY")
)

go_bp <- go_map[
  !is.na(go_map$GO) &
  go_map$ONTOLOGY == "BP",
]

# Create GO gene sets
go_sets <- split(go_bp$SYMBOL, go_bp$GO)
go_sets <- lapply(go_sets, unique)

# Keep gene sets with enough genes
go_sets <- go_sets[
  sapply(go_sets, length) >= 10
]

# Run preranked GSEA
gsea_res <- fgsea(
  pathways=go_sets,
  stats=ranks,
  minSize=10,
  maxSize=500,
  eps=0
)

# Add GO term names
go_terms <- AnnotationDbi::select(
  GO.db,
  keys=gsea_res$pathway,
  keytype="GOID",
  columns="TERM"
)

gsea_res <- merge(
  gsea_res,
  go_terms,
  by.x="pathway",
  by.y="GOID",
  all.x=TRUE
)

# Sort by adjusted p-value
gsea_res <- gsea_res[
  order(gsea_res$padj),
]

# Remove list columns before saving
gsea_export <- gsea_res

gsea_export$leadingEdge <- NULL

write.csv(
  gsea_export,
  "../results/GSEA_GO_BP_24h_vs_0h.csv",
  row.names=FALSE
)

# Select top enriched pathways# Select top enriched pathways from both directions
gsea_valid <- gsea_res[!is.na(gsea_res$NES) & !is.na(gsea_res$padj), ]

top_positive <- head(
  gsea_valid[gsea_valid$NES > 0, ][order(gsea_valid[gsea_valid$NES > 0, ]$padj), ],
  8
)

top_negative <- head(
  gsea_valid[gsea_valid$NES < 0, ][order(gsea_valid[gsea_valid$NES < 0, ]$padj), ],
  8
)

top_gsea <- rbind(top_positive, top_negative)

# Save top results
top_gsea$leadingEdge <- NULL
write.csv(
  top_gsea,
  "../results/GSEA_GO_BP_top15_24h_vs_0h.csv",
  row.names=FALSE
)

# Portfolio-style GSEA plot

top_gsea$Direction <- ifelse(
  top_gsea$NES > 0,
  "24 h enriched",
  "0 h enriched"
)

top_gsea$TERM_plot <- stringr::str_wrap(top_gsea$TERM, width=45)

p_gsea <- ggplot(
  top_gsea,
  aes(
    x=NES,
    y=reorder(TERM_plot, NES),
    size=-log10(padj),
    shape=Direction
  )
) +
  geom_point() +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(
    title="GSEA of CD8+ T-cell Activation",
    subtitle="NES < 0: enriched at 0 h | NES > 0: enriched at 24 h",
    x="Normalized Enrichment Score (NES)",
    y=NULL,
    size="-log10(FDR)",
    shape="Direction"
  ) +
  theme_bw(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    axis.text.y=element_text(size=9)
  )

ggsave(
  "../results/GSEA_GO_BP_24h_vs_0h.png",
  p_gsea,
  width=14,
  height=9,
  dpi=300,
  bg="white"
)

