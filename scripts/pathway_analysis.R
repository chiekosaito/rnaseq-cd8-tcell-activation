# GO Biological Process enrichment analysis
# GSE212353: Primary human CD8+ T cells
# Comparison: 24 h CD3/CD28 stimulation vs 0 h

library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
# Load DESeq2 results
res_df <- read.csv(
  "../results/DESeq2_24h_vs_0h_results.csv",
  stringsAsFactors=FALSE
)

# Remove genes without adjusted p-values
res_df <- res_df[!is.na(res_df$padj), ]
# Select significantly upregulated genes
up_genes <- res_df[
  res_df$padj < 0.05 &
  res_df$log2FoldChange > 1,
]

up_symbols <- unique(up_genes$gene_name)
up_symbols <- up_symbols[!is.na(up_symbols) & up_symbols != ""]
# Define background genes from all genes tested by DESeq2
background_symbols <- unique(res_df$gene_name)
background_symbols <- background_symbols[
  !is.na(background_symbols) &
  background_symbols != ""
]# Retrieve Gene Ontology annotations
go_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=background_symbols,
  keytype="SYMBOL",
  columns=c("SYMBOL", "GO", "ONTOLOGY")
)

# Keep Biological Process annotations only
go_bp <- go_map[
  !is.na(go_map$GO) &
  go_map$ONTOLOGY == "BP",
]
# Prepare GO term gene sets
go_sets <- split(go_bp$SYMBOL, go_bp$GO)
go_sets <- lapply(go_sets, unique)

# Keep only upregulated genes present in the GO background
universe <- unique(go_bp$SYMBOL)
up_in_universe <- intersect(up_symbols, universe)

N <- length(universe)
n <- length(up_in_universe)

# Hypergeometric enrichment test
go_results <- lapply(names(go_sets), function(go_id) {

  term_genes <- intersect(go_sets[[go_id]], universe)

  K <- length(term_genes)
  k <- length(intersect(term_genes, up_in_universe))

  pvalue <- phyper(
    k - 1,
    K,
    N - K,
    n,
    lower.tail=FALSE
  )

  data.frame(
    GO=go_id,
    BackgroundGenes=K,
    UpGenes=k,
    GeneRatio=k/n,
    pvalue=pvalue
  )
})

go_results <- do.call(rbind, go_results)
go_results$padj <- p.adjust(go_results$pvalue, method="BH")
# Add GO term names
library(GO.db)

go_terms <- AnnotationDbi::select(
  GO.db,
  keys=go_results$GO,
  keytype="GOID",
  columns="TERM"
)

go_results <- merge(
  go_results,
  go_terms,
  by.x="GO",
  by.y="GOID",
  all.x=TRUE
)
# Add GO term names
library(GO.db)

go_terms <- AnnotationDbi::select(
  GO.db,
  keys=go_results$GO,
  keytype="GOID",
  columns="TERM"
)

go_results <- merge(
  go_results,
  go_terms,
  by.x="GO",
  by.y="GOID",
  all.x=TRUE
)

# Keep significant GO Biological Process terms
sig_go <- go_results[
  go_results$padj < 0.05 &
  go_results$UpGenes >= 5,
]

sig_go <- sig_go[order(sig_go$padj), ]

# Save enrichment results
write.csv(
  sig_go,
  "../results/GO_BP_upregulated_24h.csv",
  row.names=FALSE
)

# Select top 15 terms for visualization
top_go <- head(sig_go, 15)

# Create dot plot
p_go <- ggplot(
  top_go,
  aes(
    x=GeneRatio,
    y=reorder(stringr::str_wrap(TERM.x, width=45), GeneRatio),
    size=UpGenes,
    color=-log10(padj)
  )
) +
  geom_point() +
  labs(
    title="GO Enrichment: Upregulated Genes",
    subtitle="24 h CD3/CD28 stimulation vs 0 h",
    x="Gene Ratio",
    y=NULL,
    size="Upregulated genes",
    color="-log10(FDR)"
  ) +
  theme_bw(base_size=12)

ggsave(
  "../results/GO_BP_upregulated_24h.png",
  p_go,
  width=16,
  height=6,
  dpi=300
)

