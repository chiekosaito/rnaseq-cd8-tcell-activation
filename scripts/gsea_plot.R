library(ggplot2)
library(stringr)

top_gsea <- read.csv("../results/GSEA_GO_BP_top15_24h_vs_0h.csv")

top_gsea$Direction <- ifelse(
  top_gsea$NES > 0,
  "24 h enriched",
  "0 h enriched"
)

top_gsea$TERM_plot <- str_wrap(top_gsea$TERM, width=45)

p <- ggplot(
  top_gsea,
  aes(
    x=NES,
    y=reorder(TERM_plot, NES),
    size=-log10(padj),
    color=Direction
  )
) +
  geom_point() +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(
    limits=c(-2.5, 2.5),
    breaks=c(-2, -1, 0, 1, 2)
  ) +
  labs(
    title="GSEA of CD8+ T-cell Activation",
    subtitle="Negative NES: 0 h enriched | Positive NES: 24 h enriched",
    x="Normalized Enrichment Score (NES)",
    y=NULL,
    size="-log10(FDR)",
    color="Direction"
  ) +
  theme_bw(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    axis.text.y=element_text(size=9)
  )

ggsave(
  "../results/GSEA_GO_BP_bidirectional.png",
  p,
  width=14,
  height=9,
  dpi=300,
  bg="white"
)

