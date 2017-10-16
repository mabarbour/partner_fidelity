set.seed(34) # make sure simulations give same output

## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(cowplot)
library(visreg)

## LOAD AND MANAGE THE DATA ----
network.df <- read_csv("t_glmm_plus_connectance.csv") %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         c.norm_deg_sp_1 = normalized_degree_sp_1 - mean(normalized_degree_sp_1),
         c.norm_deg_sp_2 = normalized_degree_sp_2 - mean(normalized_degree_sp_2),
         sc.norm_deg_sp_1 = scale(normalized_degree_sp_1),
         sc.norm_deg_sp_2 = scale(normalized_degree_sp_2),
         expected_probability = (normalized_degree_sp_1 + normalized_degree_sp_2)/2)

plot(normalized_degree_sp_1 ~ type, network.df)
plot(normalized_degree_sp_2 ~ type, network.df)

network.df %>%
  select(species_1, subtype, type, normalized_degree_sp_1) %>%
  unique() %>%
  ggplot(., aes(x = type, y = normalized_degree_sp_1, group = subtype, fill = subtype)) + geom_boxplot()

network.df %>%
  select(species_2, subtype, type, normalized_degree_sp_2) %>%
  unique() %>%
  ggplot(., aes(x = type, y = normalized_degree_sp_2, group = subtype, fill = subtype)) + geom_boxplot()

aggregate.df <- network.df %>%
  group_by(id_pair, subtype, type) %>%
  summarise(co_occur = n(), interactions = sum(connected), 
            mean_norm_deg_sp_1 = mean(normalized_degree_sp_1), 
            mean_norm_deg_sp_2 = mean(normalized_degree_sp_2)) %>%
  ungroup() %>%
  mutate(sc_mean_deg_sp_1 = (mean_norm_deg_sp_1 - mean(mean_norm_deg_sp_1))/sd(mean_norm_deg_sp_1),
         sc_mean_deg_sp_2 = (mean_norm_deg_sp_2 - mean(mean_norm_deg_sp_2))/sd(mean_norm_deg_sp_2))

agg.glm <- glm(interactions/co_occur ~ subtype*sc_mean_deg_sp_1*sc_mean_deg_sp_2, data = aggregate.df, family = "binomial",  weights = co_occur)
summary(agg.glm)

exp(coef(agg.glm)) # when a species normalized degree increases by 1 SD, the odds of an interaction occuring increase by a factor of 1.4
exp(1.72) # when a species normalized degree increases by 1 SD, the odds of an interaction occuring increase by a factor of 5.6

mean(aggregate.df$mean_norm_deg_sp_1)
sd(aggregate.df$mean_norm_deg_sp_1)
# so when a species normalized degree increases from 0.24 to 0.43, the odds of an interaction occuring increase by a factor of ~1.4 in naive model, and 5.6 in a model that accounts for random effects.

visreg(agg.glm, xvar = "sc_mean_deg_sp_1", by = 'type', scale = 'response')
visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(type = "A"), scale = 'response')
visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(type = "M"), scale = 'response')

visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(subtype = "HostParasite"), scale = 'response')
visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(subtype = "PlantHerbivore"), scale = 'response')
visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(subtype = "PlantPollinator"), scale = 'response')
visreg2d(agg.glm, xvar = "sc_mean_deg_sp_1", yvar = "sc_mean_deg_sp_2", cond = list(subtype = "PlantSeedDisperser"), scale = 'response')


interaction.formula <- brmsformula(connected ~                                # probability of species interaction
                                     subtype*sc.norm_deg_sp_1*sc.norm_deg_sp_2 +   # 3-way interaction between type and normalized degree of each species
                                     (1 | network_id) +                       # allow probability of interaction to vary among networks
                                     #(1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                                     (1 | species_1) +                        # allow probability of interaction to vary among resource species
                                     (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                                     (1 | id_pair))                           # allow probability of interaction to vary among unique species interactions
interaction.priors <- c(set_prior("normal(0,2)", class = "b", coef = "subtypePlantHerbivore"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantPollinator"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantSeedDisperser"),
                        set_prior("normal(3.4,2)", class = "b", coef = "sc.norm_deg_sp_1"),
                        set_prior("normal(3.1,2)", class = "b", coef = "sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantHerbivore:sc.norm_deg_sp_1"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantHerbivore:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantHerbivore:sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantPollinator:sc.norm_deg_sp_1"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantPollinator:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantPollinator:sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantSeedDisperser:sc.norm_deg_sp_1"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantSeedDisperser:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "subtypePlantSeedDisperser:sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "sd"))

interaction.test <- brm(interaction.formula,
                        data = network.df, 
                        family = bernoulli(),
                        prior = interaction.priors,
                        algorithm = "sampling",
                        chains = 4,
                        control = list(adapt_delta = 0.9))

# plot(interaction.test) # everything looks good
summary(interaction.test)

subtype.data <- marginal_effects(interaction.test, 
                 effects = c("sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                 conditions = data.frame(subtype = c("HostParasite","PlantHerbivore","PlantPollinator","PlantSeedDisperser"), 
                                         row.names = c("Host-Parasite","Plant-Herbivore","Plant-Pollinator","Plant-Seed Disperser")), 
                 surface = TRUE)
plot(subtype.data, ncol = 2, stype = "raster")[[1]] +
  xlab("Normalized Degree Species 1 (scaled)") +
  ylab("Normalized Degree Species 2 (scaled)")
