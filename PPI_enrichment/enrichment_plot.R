library(ggplot2)
library(dplyr)
library(data.table)

kegg.res <- fread("PPI_enrichment/NCAN_enrichment_GO.tsv",sep='\t',header=T)
dotplot <- ggplot(kegg.res, aes(
  x = strength,
  y = reorder(category, strength),
  size = Count,
  color = p.adjust
)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 10)) +
  scale_color_gradient(low = "#B3E2CD", high = "#006400", name = "-log10(p.adjust)") + 
  labs(
    x = "Count",
    y = "GO Term",
    title = "GO Enrichment Dotplot",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )
ggsave("NCAN_enrichment_GO_dotplot.pdf", width = 10, height = 6)
###############################################################################################
kegg.res <- fread("PPI_enrichment/SERPINA1_enrichment_GO.tsv",sep='\t',header=T)
dotplot <- ggplot(kegg.res, aes(
  x = strength,
  y = reorder(category, strength),
  size = Count,
  color = p.adjust
)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 10)) +
  scale_color_gradient(low = "#B3E2CD", high = "#006400", name = "-log10(p.adjust)") + 
  labs(
    x = "Count",
    y = "GO Term",
    title = "GO Enrichment Dotplot",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )
ggsave("SERPINA1_enrichment_GO_dotplot.pdf", width = 10, height = 8)
