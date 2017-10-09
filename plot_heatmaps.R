## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(cowplot)

## GENERATE HEATMAPS ----

# load data
heatmap.data <- read_csv("model_output/heatmap_data.csv")
raw.data <- read_csv("t_glmm_plus_connectance.csv")

# transform to original scale for plotting purposes
heatmap.data$sc.norm_deg_sp_1 <- heatmap.data$sc.norm_deg_sp_1*sd(raw.data$normalized_degree_sp_1) + mean(raw.data$normalized_degree_sp_1)
heatmap.data$sc.norm_deg_sp_2 <- heatmap.data$sc.norm_deg_sp_2*sd(raw.data$normalized_degree_sp_2) + mean(raw.data$normalized_degree_sp_2)
heatmap.data$type <- ifelse(heatmap.data$type == "M", "Mutualistic", "Antagonistic")

# color scheme for plot
viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
} # taken from https://github.com/paul-buerkner/brms/blob/master/R/plots.R

# heatmap of antagonistic vs. mutualistic interactions
type.norm_sp1.norm_sp2 <- ggplot(heatmap.data, aes(x = sc.norm_deg_sp_1, y = sc.norm_deg_sp_2, fill = estimate__)) +
  geom_raster() + 
  facet_wrap(~type) +
  scale_fill_gradientn(colors = viridis6(), name = "Probability of\n Interaction") +
  xlab("Resource Normalized Degree") +
  ylab("Consumer Normalized Degree") +
  bayesplot::theme_default()
ggsave("plots/type.norm_sp1.norm_sp2.pdf", width = 11, height = 8.5, units = "in")
  
## PLOT THE DIFFERENCE IN PREDICTIONS ----

# calculate the differences
ant <- heatmap.data %>%
  filter(type == "Antagonistic") %>%
  select(sc.norm_deg_sp_1, sc.norm_deg_sp_2, estimate.Ant = estimate__)

mut <- heatmap.data %>%
  filter(type == "Mutualistic") %>%
  select(sc.norm_deg_sp_1, sc.norm_deg_sp_2, estimate.Mut = estimate__)

diff.data <- left_join(mut, ant) %>%
  mutate(estimate.Diff = estimate.Mut - estimate.Ant)

# generate the plot
diff.plot <- ggplot(diff.data, aes(x = sc.norm_deg_sp_1, y = sc.norm_deg_sp_2, fill = estimate.Diff)) +
  geom_raster() +
  scale_fill_gradientn(colors = viridis6(), name = "Difference in \nInteraction Probability") + 
  xlab("Resource Normalized Degree") +
  ylab("Consumer Normalized Degree") +
  bayesplot::theme_default()
ggsave("plots/diff.plot.pdf")
