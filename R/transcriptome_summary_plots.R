###### Transcriptome summary plots ######
# This script contains code used to make the PCA in figure 1B as well as the 
# Venn diagram in figure 4B.

library(ggplot2)
library(ggrepel)
library(VennDiagram)

###### Making the PCA plot ######
### Importing the data ###
pca_plot_data <- read.table(file = "./data/raw_data/PCA_plot_data.tsv",
                            header = TRUE,
                            sep = '\t')
pca_plot_data$group <- factor(pca_plot_data$group)


# Color pallet
cbPalette <- c("#E69F00", "#2557a8", "#009E73", "#a05d82")

# Plot (export size 2000x1400)
ggplot(pca_plot_data, aes(x = PC1, y = PC2, label = rownames(pca_plot_data), color = group)) + 
  geom_point(size = 5) +
  geom_text_repel(size = 12, segment.color = NA, point.padding = 0.4) +
  stat_ellipse(level = .95, size = 2) +
  labs(title = "Transcriptome PCA",
       x = "PC1 (29% of variance)",
       y = "PC2 (20.8% of variance)") +
  scale_color_manual(values = cbPalette) +
  theme_classic() +
  theme(axis.title = element_text(size = 52, face = "bold"),
        axis.text = element_text(size = 44, face = "bold"),
        axis.title.x = element_text(margin = margin(30, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 30, 0, 0)),
        plot.title = element_text(size = 72, face = "bold", margin = margin(0, 0, 40, 0)),
        axis.line = element_line(size = 2, color = "#4D4D4D"),
        axis.ticks = element_line(size = 2, color = "#4D4D4D"),
        axis.ticks.length = unit(10, "pt"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.position = "none")

###### Making the Venn diagram ######
dev.off() # Fixes a bug of plotting on top of the previous plot

# Making the Venn diagram (export size 1100x1100)
draw.quad.venn(area1 = 200, area2 = 200, area3 = 200, area4 = 200,
               n12 = 11, n13 = 66, n14 = 9, n23 = 13, n24 = 47, n34 = 19,
               n123 = 0, n124 = 4, n134 = 3, n234 = 9,
               n1234 = 0,
               fill = c("#9F5C82", "#019E73", "#2558A8", "#E69F00"),
               lwd = rep(3, 4),
               cex = 3.5,
               fontface = "bold",
               fontfamily = "sans",
               cat.cex = 2,
               cat.fontface = "bold",
               cat.fontfamily = "sans",
               margin = 0.24)