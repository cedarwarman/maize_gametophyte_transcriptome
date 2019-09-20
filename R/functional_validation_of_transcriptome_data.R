###### Functional validation of sequencing data ######
# Here we create generalized linear models describing transmission defects, 
# then plot the results (Figure 6)

library(tidyverse)

###### Importing the data ######
male_df <- read.table(file = "./data/raw_data/male_seed_counts.tsv",
                 header = TRUE,
                 sep = '\t')
male_df$female_family <- as.factor(male_df$female_family)

female_df <- read.table(file = "./data/raw_data/female_seed_counts.tsv",
                        header = TRUE,
                        sep = '\t')

# Splitting the male crosses into expression catagories
all_exp <- male_df
vc_high <- male_df[male_df$expression_category_a == "vegetative_cell_high", ]
sc_high <- male_df[male_df$expression_category_a == "sperm_cell_high", ]
seedling_only <- male_df[male_df$expression_category_a == "seedling_only", ]

###### Generalized linear models (GLM) ######
# Here we create generalized linear models with a quasibinomial distribution 
# and a logit link function

### GLM for Vegetative Cell category ###
glm_vc_quasi <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                          as.factor(male_family_number),
                          data = vc_high,
                          family = quasibinomial(link = "logit"))
summary(glm_vc_quasi)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_vc_quasi)
p_vc = coef(summary(glm_vc_quasi))[1,4]
for (i in 2:nrow(coef(summary(glm_vc_quasi)))){
  beta = coef(summary(glm_vc_quasi))[1,1] + coef(summary(glm_vc_quasi))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_vc = c(p_vc, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_vc) = sort(unique(vc_high$male_family_number))

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_vc_adj = p.adjust(p_vc,method = "BH")
p_results_vc=data.frame(raw_p = p_vc, adjusted_p = p_vc_adj)
p_results_vc


### GLM for Sperm Cell category ###
glm_sc_quasi <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                          as.factor(male_family_number),
                          data = sc_high,
                          family = quasibinomial(link = "logit"))
summary(glm_sc_quasi)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_sc_quasi)
p_sc = coef(summary(glm_sc_quasi))[1,4]
for (i in 2:nrow(coef(summary(glm_sc_quasi)))){
  beta = coef(summary(glm_sc_quasi))[1,1] + coef(summary(glm_sc_quasi))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_sc = c(p_sc, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_sc) = sort(unique(sc_high$male_family_number))

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_sc_adj = p.adjust(p_sc,method="BH")
p_results_sc=data.frame(raw_p = p_sc, adjusted_p = p_sc_adj)
p_results_sc


### GLM for Seedling Only category ###
glm_seed_quasi <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                               as.factor(male_family_number),
                               data = seedling_only,
                               family = quasibinomial(link = "logit"))
summary(glm_seed_quasi)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_seed_quasi)
p_seed = coef(summary(glm_seed_quasi))[1,4]
for (i in 2:nrow(coef(summary(glm_seed_quasi)))){
  beta = coef(summary(glm_seed_quasi))[1,1] + coef(summary(glm_seed_quasi))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_seed = c(p_seed, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_seed) = sort(unique(seedling_only$male_family_number))

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_seed_adj = p.adjust(p_seed,method="BH")
p_results_seed=data.frame(raw_p = p_seed, adjusted_p = p_seed_adj)
p_results_seed


### GLM for female crosses, all categories ###
glm_female_quasi <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                        as.factor(female_family),
                        data = female_df,
                        family = quasibinomial(link = "logit"))
summary(glm_female_quasi)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_female_quasi)
p_female = coef(summary(glm_female_quasi))[1,4]
for (i in 2:nrow(coef(summary(glm_female_quasi)))){
  beta = coef(summary(glm_female_quasi))[1,1] + coef(summary(glm_female_quasi))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_female = c(p_female, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_female) = sort(unique(female_df$female_family))

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_female_adj = p.adjust(p_female,method="BH")
p_results_female=data.frame(raw_p = p_female, adjusted_p = p_female_adj)
p_results_female


###### Proportion tests ######
# Here we check to see if the proportion of genes with a transmission defect
# is statistically different between each category (Seedling, Vegetative Cell, 
# Sperm Cell).

# Comparing the Vegetative Cell and Sperm Cell categories
df <- matrix(c(1, 7, 9, 25), nrow = 2, dimnames = list(c("sc", "vc"), c("defect", "no_defect")))
df
fisher.test(df, alternative = 'l')

# Comparing Seedling to Sperm Cell
df <- matrix(c(0, 1, 10, 9), nrow = 2, dimnames = list(c("seedling", "sc"), c("defect", "no_defect")))
df
fisher.test(df, alternative = 'l')

# Comparing Seedling to Vegetative Cell
df <- matrix(c(0, 7, 10, 25), nrow = 2, dimnames = list(c("seedling", "vc"), c("defect", "no_defect")))
df
fisher.test(df, alternative = 'l')

# Comparing low and high expression in the Vegetative Cell category (note, 
# relatively low, all genes in this portion of the study were highly expressed)
df <- matrix(c(1, 6, 19, 6), nrow = 2, dimnames = list(c("low_exp", "high_exp"), c("defect", "no_defect")))
df
fisher.test(df)


###### Analysis of Ac impact on transmission rate ######
# Here we test the effect of Ac presence by comparing transmission rates in 
# alleles that had both Ac present and absent.

### Importing the data ###
# Importing the data. This df was created from Supplemental Table 7. It 
# contains only families (8 of them) that have lines with both Ac and no Ac.
ac_df <- read.table(file = "./data/raw_data/families_with_both_ac_and_no_ac.tsv",
                    header = TRUE,
                    sep = '\t')

### Generalized linear model ###
# Both Ac and family as factors
glm_ac_no_ac <- glm(cbind(Marker_seeds, Non_marker_seeds) ~
                    as.factor(paternal_family) +
                    as.factor(Ac_present),
                    data = ac_df,
                    family = quasibinomial(link = "logit"))
summary(glm_ac_no_ac)

# Getting the proper p-values (p-values of individual parameters comparing 
# the transmission to Mendelian, 50%)
S = vcov(glm_ac_no_ac)
p_ac = coef(summary(glm_ac_no_ac))[1,4]
for (i in 2:nrow(coef(summary(glm_ac_no_ac)))){
  beta = coef(summary(glm_ac_no_ac))[1,1] + coef(summary(glm_ac_no_ac))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_ac = c(p_ac, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_ac) = sort(unique(ac_df$paternal_family))
names(p_ac)[[9]] <- 'ac_present'

# Multiple testing correction using Benjamini-Hochberg
p.adj.ac = p.adjust(p_ac,method="BH")
p_results_ac=data.frame(raw_p=p_ac, adjusted_p=p.adj.ac)

# Ac present is not significant when family is taken into account
p_results_ac


###### Building dataframes for plotting ######
### Setting up the female plot data ###
# Imported metadata includes line/allele information as well as the raw seed 
# counts and FPKM for the group that they belong to. 
plot_females <- read.table(file = "./data/raw_data/plot_female_metadata.tsv",
                           header = TRUE,
                           sep = '\t')
plot_females$female_family <- as.character(plot_females$female_family)

# Preparing the p-values from the GLMs
colnames(p_results_female) <- c("glm_ql_raw_p", "glm_ql_adj_p")
p_results_female <- rownames_to_column(p_results_female, "female_family")

# Joining the p-values to the metadata
plot_females <- full_join(plot_females, p_results_female)

# Adding a factor for p-value significance
p_value_threshold <- 0.05
plot_females$glm_ql_adj_p_sig <- NA
plot_females$glm_ql_adj_p_sig[plot_females$glm_ql_adj_p >= p_value_threshold] <- 0
plot_females$glm_ql_adj_p_sig[plot_females$glm_ql_adj_p < p_value_threshold] <- 1


### Setting up the Seedling category plot data ###
plot_seedling <- read.table(file = "./data/raw_data/plot_seedling_metadata.tsv",
                            header = TRUE,
                            sep = '\t')

# Preparing the p-values from the GLMs
colnames(p_results_seed) <- c("glm_ql_raw_p", "glm_ql_adj_p")
p_results_seed <- rownames_to_column(p_results_seed, "male_family_number")

# Joining the p-values to the metadata
plot_seedling$male_family_number <- as.character(plot_seedling$male_family_number)
plot_seedling <- full_join(plot_seedling, p_results_seed)

# Adding a factor for p-value significance
p_value_threshold <- 0.05
plot_seedling$glm_ql_adj_p_sig <- NA
plot_seedling$glm_ql_adj_p_sig[plot_seedling$glm_ql_adj_p >= p_value_threshold] <- 0
plot_seedling$glm_ql_adj_p_sig[plot_seedling$glm_ql_adj_p < p_value_threshold] <- 1


### Setting up the Sperm Cell category plot data ###
plot_sc <- read.table(file = "./data/raw_data/plot_sc_metadata.tsv",
                      header = TRUE,
                      sep = '\t')

# Preparing the p-values from the GLMs
colnames(p_results_sc) <- c("glm_ql_raw_p", "glm_ql_adj_p")
p_results_sc <- rownames_to_column(p_results_sc, "male_family_number")

# Joining the p-values to the metadata
plot_sc$male_family_number <- as.character(plot_sc$male_family_number)
plot_sc <- full_join(plot_sc, p_results_sc)

# Adding a factor for p-value significance
p_value_threshold <- 0.05
plot_sc$glm_ql_adj_p_sig <- NA
plot_sc$glm_ql_adj_p_sig[plot_sc$glm_ql_adj_p >= p_value_threshold] <- 0
plot_sc$glm_ql_adj_p_sig[plot_sc$glm_ql_adj_p < p_value_threshold] <- 1


### Setting up the Vegetative Cell category plot data ###
plot_vc <- read.table(file = "./data/raw_data/plot_vc_metadata.tsv",
                      header = TRUE,
                      sep = '\t')

# Preparing the p-values from the GLMs
colnames(p_results_vc) <- c("glm_ql_raw_p", "glm_ql_adj_p")
p_results_vc <- rownames_to_column(p_results_vc, "male_family_number")

# Joining the p-values to the metadata
plot_vc$male_family_number <- as.character(plot_vc$male_family_number)
plot_vc <- full_join(plot_vc, p_results_vc)

# Adding a factor for p-value significance
p_value_threshold <- 0.05
plot_vc$glm_ql_adj_p_sig <- NA
plot_vc$glm_ql_adj_p_sig[plot_vc$glm_ql_adj_p >= p_value_threshold] <- 0
plot_vc$glm_ql_adj_p_sig[plot_vc$glm_ql_adj_p < p_value_threshold] <- 1

###### Vegetative Cell category regressions ######
# Here we test if expression level is associated with percent transmission or 
# p-value in genes belonging to the Vegetative Cell category
vc_regressions <- plot_vc

### -log(p-value) by log(FPKM) ###
# Here we compare the -log10(p-value) to log2(FPKM)

# Due to a rounding error, on p-value is exactly 0. Because we take the -log10
# of the p-values, we convert this p-value to a very small number to avoid the 
# error.
vc_regressions$glm_ql_adj_p[vc_regressions$glm_ql_adj_p == 0] <- 0.000000000000000000001

fpkm_p_lm <- lm(-log10(glm_ql_adj_p) ~ log2(FPKM), data = vc_regressions)
summary(fpkm_p_lm)

### log(FPKM) by estimate ###
# Here we compare the log2(FPKM) to the estimaged transmission rate.
scatter.smooth(x = log2(vc_regressions$FPKM), y = vc_regressions$percent_GFP)

fpkm_est_lm <- lm(percent_GFP ~ log2(FPKM), data = vc_regressions)
summary(fpkm_est_lm)


###### Plotting females ######
color_vec <- c("#2d2d2d", "#e56e67") # Black and salmon

ggplot(data = plot_females, aes(x = log2(FPKM), y = 100 * percent_GFP, color = factor(glm_ql_adj_p_sig), fill = factor(glm_ql_adj_p_sig))) +
  geom_hline(yintercept = 50, size = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 21, stroke = 1, color = "black") +
  scale_y_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55),
                     labels = c("20", "25", "30", "35", "40", "45", "50", "55"),
                     limits = c(20, 56)) +
  scale_x_continuous(breaks = seq(5, 15, by = 2.5),
                     labels = as.character(seq(5, 15, by = 2.5)),
                     limits = c(4.8, 15)) +
  labs(title = "Female transmission rates (all alleles)",
       x = "log2(FPKM)",
       y = "% marker transmission") +
  scale_fill_manual(values = color_vec,
                    name = "Significance",
                    breaks = c("0", "1"),
                    labels = c(" p > 0.05", " p < 0.05")) +
  theme_bw() +
  theme(plot.title = element_text(size = 29, face = "bold", margin = margin(0, 0, 12, 0)),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.line = element_line(size = 0, color = "#4D4D4D"),
        axis.ticks = element_line(size = 1, color = "#4D4D4D"),
        axis.ticks.length = unit(4, "pt"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
        panel.border = element_rect(color = "#4D4D4D", size = 2, fill = NA),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

# Saving the plot
ggsave(filename = './data/plots/female_transmission.png', 
       device = 'png', 
       width = 9, 
       height = 6, 
       dpi = 400,
       units = 'in')


###### Plotting Seedling category ######
color_vec <- c("#2d2d2d", "#e56e67")

ggplot(data = plot_seedling, aes(x = log2(FPKM), y = 100 * percent_GFP, color = factor(glm_ql_adj_p_sig), fill = factor(glm_ql_adj_p_sig))) +
  geom_hline(yintercept = 50, size = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 21, stroke = 1, color = "black") +
  scale_y_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55),
                     labels = c("20", "25", "30", "35", "40", "45", "50", "55"),
                     limits = c(20, 55)) +
  scale_x_continuous(breaks = seq(5, 8, 1),
                     labels = seq(5, 8, 1), 
                     limits = c(5, 8)) +
  labs(title = "Seedling: male transmission rates",
       x = "log2(FPKM)",
       y = "% marker transmission") +
  scale_fill_manual(values = color_vec,
                    name = "Significance",
                    breaks = c("0", "1"),
                    labels = c(" p > 0.05", " p < 0.05")) +
  theme_bw() +
  theme(plot.title = element_text(size = 29, face = "bold", margin = margin(0, 0, 12, 0)),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.line = element_line(size = 0, color = "#4D4D4D"),
        axis.ticks = element_line(size = 1, color = "#4D4D4D"),
        axis.ticks.length = unit(4, "pt"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
        panel.border = element_rect(color = "#4D4D4D", size = 2, fill = NA),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

# Saving the plot
ggsave(filename = './data/plots/seedling_transmission.png', 
       device = 'png', 
       width = 9, 
       height = 6, 
       dpi = 400,
       units = 'in')


###### Plotting Sperm Cell category ######
color_vec <- c("#2d2d2d", "#e56e67")

ggplot(data = plot_sc, aes(x = log2(FPKM), y = 100 * percent_GFP, color = factor(glm_ql_adj_p_sig), fill = factor(glm_ql_adj_p_sig))) +
  geom_hline(yintercept = 50, size = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 21, stroke = 1, color = "black") +
  scale_y_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55),
                     labels = c("20", "25", "30", "35", "40", "45", "50", "55"),
                     limits = c(20, 55)) +
  scale_x_continuous(breaks = seq(5, 9, 1),
                     labels = seq(5, 9, 1),
                     limits = c(5, 9)) +
  labs(title = "Sperm cell: male transmission rates",
       x = "log2(FPKM)",
       y = "% marker transmission") +
  scale_fill_manual(values = color_vec,
                    name = "Significance",
                    breaks = c("0", "1"),
                    labels = c(" p > 0.05", " p < 0.05")) +
  theme_bw() +
  theme(plot.title = element_text(size = 29, face = "bold", margin = margin(0, 0, 12, 0)),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.line = element_line(size = 0, color = "#4D4D4D"),
        axis.ticks = element_line(size = 1, color = "#4D4D4D"),
        axis.ticks.length = unit(4, "pt"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
        panel.border = element_rect(color = "#4D4D4D", size = 2, fill = NA),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

# Saving the plot
ggsave(filename = './data/plots/sc_transmission.png', 
       device = 'png', 
       width = 9, 
       height = 6, 
       dpi = 400,
       units = 'in')


###### Plotting Vegetative Cell category #######


color_vec <- c("#2d2d2d", "#e56e67")

# (export 2200x1400, 2400x1300 for poster)
ggplot(data = plot_vc, aes(x = log2(FPKM), y = 100 * percent_GFP, color = factor(glm_ql_adj_p_sig), fill = factor(glm_ql_adj_p_sig))) +
  geom_hline(yintercept = 50, size = 1, linetype = "dashed") +
  geom_point(size = 5, shape = 21, stroke = 1, color = "black") +
  scale_y_continuous(breaks = c(20, 25, 30, 35, 40, 45, 50, 55),
                     labels = c("20", "25", "30", "35", "40", "45", "50", "55"),
                     limits = c(20, 55)) +
  scale_x_continuous(breaks = seq(5, 15, 2),
                     labels = seq(5, 15, 2),
                     limits = c(4.9, 15)) +
  labs(title = "Vegetative cell: male transmission rates",
       x = "log2(FPKM)",
       y = "% marker transmission") +
  scale_fill_manual(values = color_vec,
                    name = "Significance:",
                    breaks = c("0", "1"),
                    labels = c(" p > 0.05", " p < 0.05")) +
  theme_bw() +
  theme(plot.title = element_text(size = 29, face = "bold", margin = margin(0, 0, 12, 0)),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.line = element_line(size = 0, color = "#4D4D4D"),
        axis.ticks = element_line(size = 1, color = "#4D4D4D"),
        axis.ticks.length = unit(4, "pt"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
        panel.border = element_rect(color = "#4D4D4D", size = 2, fill = NA),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # legend.title = element_text(size = 22, face = "bold"),
        # legend.text = element_text(size = 22, face = "bold"),
        # legend.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
        # legend.position = "bottom",
        legend.position = "none")

# Saving the plot
ggsave(filename = './data/plots/vc_transmission.png', 
       device = 'png', 
       width = 9, 
       height = 6, 
       dpi = 400,
       units = 'in')

