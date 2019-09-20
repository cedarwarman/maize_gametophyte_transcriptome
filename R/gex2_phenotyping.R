###### gex2 phenotyping ######
# This code examines gex2 undeveloped seed and seedless area phenotypes.

library(tidyverse)
library(googlesheets)

###### gex2 small seed counts #######
### Importing and tidying the data ###
seed_df <- read.table(file = "./data/raw_data/gex2_small_seeds.tsv",
                      header = TRUE,
                      sep = '\t')

# Rearranging the factors
seed_df$group <- factor(seed_df$group, 
                        levels = c("gfp-control1", 
                                   "gfp-control2", 
                                   "gfp-pollenMutant", 
                                   "gex2_R82A03", 
                                   "gex2_R84A12", 
                                   "homo"))

### Small seeds stats ###
# Pairwise t-test
pairwise.t.test(seed_df$small_seed_percent, seed_df$group, p.adjust.method = "none")

### Small seeds box plot ###
# Plotting (export 1900 x 1700)
ggplot(seed_df, aes(x = group, y = small_seed_percent, fill = group)) +
  geom_boxplot(size = 2, fill = "#dddddd", outlier.size = 5) +
  scale_x_discrete(labels = c("GFP line 1\n(tdsgR12H07)", 
                              "GFP line 2\n(tdsgR46C04)", 
                              "VC mutant\n(tdsgR96C12)", 
                              "gex2 het.\n(tdsgR82A03)", 
                              "gex2 het.\n(tdsgR84A12)", 
                              "gex2 homo.\n(tdsgR84A12)")) +
  labs(title = "Frequency of small or aborted seeds", y = "Percent small or aborted seeds") +
  scale_y_continuous(breaks = seq(0, 30, by = 5),
                     labels = seq(0, 30, by = 5),
                     limits = c(0, 30)) +
  theme_bw() +
  theme(axis.title = element_text(size = 52, face = "bold"),
        axis.text.x = element_text(size = 45, face = "bold", color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 50, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 30, 0, 0)),
        plot.title = element_text(size = 72, face = "bold", margin = margin(0, 0, 40, 0)), 
        axis.line = element_line(size = 2, color = "black"),
        axis.ticks = element_line(size = 2, color = "black"),
        axis.ticks.length = unit(10, "pt"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

###### gex2 seedless area ######
### Importing and tidying the data ###
area_df <- read.table(file = "./data/raw_data/gex2_seedless_area.tsv",
                            header = TRUE,
                            sep = '\t')

# Rearranging the factors
area_df$group <- factor(area_df$group, 
                        levels = c("gfp-control1", 
                                   "gfp-control2", 
                                   "gfp-pollenMutant", 
                                   "gex2_R82A03", 
                                   "gex2_R84A12", 
                                   "homo"))

### Seedless area stats ###
# Pairwise t-test for area in the different groups
pairwise.t.test(area_df$empty_area_percent, area_df$group, p.adjust.method = 'none')

### Seedless area box plot ###
# Plotting (export 1900 x 1700)
ggplot(area_df, aes(x = group, y = empty_area_percent, fill = group)) +
  geom_boxplot(size = 2, fill = "#dddddd", outlier.size = 5) +
  scale_x_discrete(labels = c("GFP line 1\n(tdsgR12H07)", 
                              "GFP line 2\n(tdsgR46C04)", 
                              "VC mutant\n(tdsgR96C12)", 
                              "gex2 het.\n(tdsgR82A03)", 
                              "gex2 het.\n(tdsgR84A12)", 
                              "gex2 homo.\n(tdsgR84A12)")) +
  labs(title = "Approximate seedless ear area", y = "Approx. percent seedless area") +
  theme_bw() +
  theme(axis.title = element_text(size = 52, face = "bold"),
        axis.text = element_text(size = 50, face = "bold"),
        axis.text.x = element_text(color = "black", angle = 40, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 30, 0, 0)),
        plot.title = element_text(size = 72, face = "bold", margin = margin(0, 0, 40, 0)), 
        axis.line = element_line(size = 2, color = "black"),
        axis.ticks = element_line(size = 2, color = "black"),
        axis.ticks.length = unit(10, "pt"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


