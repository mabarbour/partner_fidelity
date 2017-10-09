
set.seed(34) # make sure simulations give same output

## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)
library(cowplot)

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


## SET NUMBER OF CORES FOR MODEL ----
options(mc.cores=parallel::detectCores ()) # Run on multiple cores

## EXPLORE THE DATA ----
table(network.df$subtype) # dataset dominated by host-parasite and plant-pollinator interactions

table(network.df$network_id) # 136 unique networks, substantial heterogeneity in which interactions come from each network
summary(c(table(network.df$network_id)))

table(network.df$id_pair) # 4075 unique interactions
summary(c(table(network.df$id_pair))) # most interactions co-occur 2-3 times, but the max is 19.

# both species 1 and 2 have similar range of normalized degree
summary(network.df$normalized_degree_sp_1)
summary(network.df$normalized_degree_sp_2)

mean.norm_deg_sp_1 <- mean(network.df$normalized_degree_sp_1)
mean.norm_deg_sp_2 <- mean(network.df$normalized_degree_sp_2)

## Explore the relationship between normalized degree and interaction type

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

ggplot(network.df, aes(x = normalized_degree_sp_1, y = connected, color = type)) +
  geom_point() + 
  binomial_smooth() 
ggplot(network.df, aes(x = normalized_degree_sp_2, y = connected, color = type)) +
  geom_point() + 
  binomial_smooth() 
ggplot(network.df %>% mutate(connected = as.factor(connected)), 
       aes(x = normalized_degree_sp_2, y = normalized_degree_sp_1, fill = connected, color = connected)) +
  geom_point(shape = 21, size = 2) + 
  facet_wrap(~type) + 
  scale_fill_manual(values = c("white", "black")) + 
  scale_color_manual(values = c("grey", "black")) 

# focus on mutualists
ggplot(network.df %>% mutate(connected = as.factor(connected), 
                             sp_1_type = factor(ifelse(normalized_degree_sp_1 > mean.norm_deg_sp_1, "sp_1_generalist", "sp_1_specialist"), levels = c("sp_1_generalist", "sp_1_specialist")),
                             sp_2_type = factor(ifelse(normalized_degree_sp_2 > mean.norm_deg_sp_2, "sp_2_generalist", "sp_2_specialist"), levels = c("sp_2_specialist", "sp_2_generalist"))) %>%
         filter(type == "M"), 
       aes(x = normalized_degree_sp_1, y = normalized_degree_sp_2, fill = connected, color = connected)) +
  geom_point(shape = 21, size = 4) + 
  facet_grid(sp_1_type ~ sp_2_type, scales = "free") + # facet_wrap()
  scale_fill_manual(values = c("white", "black")) + 
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw()


# focus on antagonists
ggplot(network.df %>% mutate(connected = as.factor(connected), 
                             sp_1_type = factor(ifelse(normalized_degree_sp_1 > mean.norm_deg_sp_1, "sp_1_generalist", "sp_1_specialist"), levels = c("sp_1_generalist", "sp_1_specialist")),
                             sp_2_type = factor(ifelse(normalized_degree_sp_2 > mean.norm_deg_sp_2, "sp_2_generalist", "sp_2_specialist"), levels = c("sp_2_specialist", "sp_2_generalist"))) %>%
         filter(type == "A"), 
       aes(x = normalized_degree_sp_1, y = normalized_degree_sp_2, fill = connected, color = connected)) +
  geom_point(shape = 21, size = 4) + 
  facet_grid(sp_1_type ~ sp_2_type, scales = "free") + # facet_wrap()
  scale_fill_manual(values = c("white", "black")) + 
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw()

# from species 1 perspective
ggplot(network.df %>% mutate(connected = connected, 
                             sp_1_type = factor(ifelse(normalized_degree_sp_1 > mean.norm_deg_sp_1, "sp_1_generalist", "sp_1_specialist"), levels = c("sp_1_generalist", "sp_1_specialist")),
                             sp_2_type = factor(ifelse(normalized_degree_sp_2 > mean.norm_deg_sp_2, "sp_2_generalist", "sp_2_specialist"), levels = c("sp_2_specialist", "sp_2_generalist"))),
       aes(x = normalized_degree_sp_1, y = connected, fill = sp_2_type, color = sp_2_type)) +
  geom_point(shape = 21, size = 4) + 
  binomial_smooth() +
  facet_wrap(~type, scales = "free") 

# from species 2 perspective
ggplot(network.df %>% mutate(connected = connected, 
                             sp_1_type = factor(ifelse(normalized_degree_sp_1 > mean.norm_deg_sp_1, "sp_1_generalist", "sp_1_specialist"), levels = c("sp_1_generalist", "sp_1_specialist")),
                             sp_2_type = factor(ifelse(normalized_degree_sp_2 > mean.norm_deg_sp_2, "sp_2_generalist", "sp_2_specialist"), levels = c("sp_2_specialist", "sp_2_generalist"))),
       aes(x = normalized_degree_sp_2, y = connected, fill = sp_1_type, color = sp_1_type)) +
  geom_point(shape = 21, size = 4) + 
  binomial_smooth() +
  facet_wrap(~type, scales = "free")

## BAYESIAN MODEL OF INTERACTION PROBABILITY ----
interaction.formula <- brmsformula(connected ~                                # probability of species interaction
                                     type*sc.norm_deg_sp_1*sc.norm_deg_sp_2 +   # 3-way interaction between type and normalized degree of each species
                                     (1 | network_id) +                       # allow probability of interaction to vary among networks
                                     (1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                                     (1 | species_1) +                        # allow probability of interaction to vary among resource species
                                     (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                                     (1 | id_pair))                           # allow probability of interaction to vary among unique species interactions


## RATIONALE FOR MODEL PRIORS ----

get_prior(interaction.formula, data = network.df, family = bernoulli()) # see which priors I need to set.

## INTERACTION TYPE

# playing around with different hypothetical probabilities suggests to me that the logistic regression coefficient will be 
# between -5 and 5. Therefore, a normal prior with mean 0 and sd = 2, seems reasonable to me.

# Let's take a hypothetical example, where the effect of interaction type is really large, ranging from
# a probability of 0.25 to 0.75. This corresponds to a logistic regression coefficient of ~2.2
pM <- 0.75 # probability of a mutualistic interaction when two species co-occur 
pA <- 0.25 # probability of an antagonistic interaction when two species co-occur 
mu_mutualist <- log(pM / (1-pM)) - log(pA / (1-pA)) # logistic regression coefficient, which is the difference in log-odds for 1 unit increase in the predictor

# I would expect such a large effect to be unlikely though. Also, whether mutualism or antagonism has a 
# positive or negative effect is unclear to me. I could imagine it either way. So let's specify a 
# regularizing prior where the prior for the coefficient is centered on zero, but the variance is large
# enough to allow for bigger effects if there is enough evidence to support them.
samp <- 1000
sd_type <- data.frame(density = c(rnorm(samp, 0, 0.5),
                                  rnorm(samp, 0, 1),
                                  rnorm(samp, 0, 2),
                                  rnorm(samp, 0, 4)),
                      standard_deviation = c(rep("sd_0.5", samp),
                                             rep("sd_1", samp),
                                             rep("sd_2", samp),
                                             rep("sd_4", samp)))
ggplot(sd_type, aes(x = density, fill = standard_deviation)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = mu_mutualist, linetype = "dotted") # big effect

# I think a normal prior with mean = 0 and sd = 2 would be appropriate.

## NORMALIZED DEGREES

# I expect for the normalized degree of a species to always have a positive relationship with the 
# probability of an interaction. Biologically, it does not make sense that this relationship could 
# ever be negative, which is why I will impose a lower bound of zero on this parameter. By default,
# I would expect this relationship to be 1:1 between the probability of observing an interaction and 
# the normalized degree of a species. Since the normalized degree of a species defines its probability
# of interacting with a co-occuring partner, regardless of its identity, my prior expectation is that 
# there is a 1:1 relationship between the normalized degree of a species and its probability of interacting
# with another species. In other words, the probability of observing an interaction is the same as
# the normalized degree. Let's get some intuition as to what this prior looks like when normalized degree
# is on a standardized scale (1 SD above and below the mean).
pSpec_sp1 <- mean(network.df$normalized_degree_sp_1) - sd(network.df$normalized_degree_sp_1) # probability of a relatively specialized species interacting
pGen_sp1 <- mean(network.df$normalized_degree_sp_1) + sd(network.df$normalized_degree_sp_1) # probability of a relatively generalized species interacting

mu_sc.norm_deg_sp1 <- log(pGen_sp1 / (1-pGen_sp1)) - log(pSpec_sp1 / (1-pSpec_sp1)) # logistic regression coefficient, which is the difference in log-odds for 1 standardized unit increase in a species normalized degree.

pSpec_sp2 <- mean(network.df$normalized_degree_sp_2) - sd(network.df$normalized_degree_sp_2) # probability of a relatively specialized species interacting
pGen_sp2 <- mean(network.df$normalized_degree_sp_2) + sd(network.df$normalized_degree_sp_2) # probability of a relatively generalized species interacting

mu_sc.norm_deg_sp2 <- log(pGen_sp2 / (1-pGen_sp2)) - log(pSpec_sp2 / (1-pSpec_sp2)) # logistic regression coefficient, which is the difference in log-odds for 1 standardized unit increase in a species normalized degree.

# Now, this is only a prior expectation, and to reflect some uncertainty around this prior, I'm going
# to specify a normal distribution centered on this mean and explore different variances
samp <- 1000
sd_sp1 <- data.frame(density = c(rnorm(samp, mu_sc.norm_deg_sp1, 0.5),
                                 rnorm(samp, mu_sc.norm_deg_sp1, 1),
                                 rnorm(samp, mu_sc.norm_deg_sp1, 2),
                                 rnorm(samp, mu_sc.norm_deg_sp1, 4)),
                     standard_deviation = c(rep("sd_0.5", samp),
                                            rep("sd_1", samp),
                                            rep("sd_2", samp),
                                            rep("sd_4", samp)))
ggplot(sd_sp1, aes(x = density, fill = standard_deviation)) + geom_density(alpha = 0.5)

# Most standardized logistic regression coefficients are less than 5, so I feel that specifying a standard
# deviation of 2 is reasonable since it more than covers the range from 0 to 5 (and I never expect this coefficient to be negative). It will
# also make the data "work" for larger values.

# For interactions involving interaction type, I'm going to specify the same regularizing prior
# as for "type", because I could see interaction type having positive or negative effects on this relationship
# and I want to make the data "work" for those strong effects.
# For the interaction with normalized degree, I expect this relationship to always be positive, but I don't know how much
# so I'm going to set a lower bound of zero, but I'll increase the standard deviation to 2 to allow more flexibility
# in having positive values.


## SET PRIORS ----
interaction.priors <- c(set_prior("normal(0,2)", class = "b", coef = "typeM"),
                        set_prior("normal(3.4,2)", class = "b", coef = "sc.norm_deg_sp_1"),
                        set_prior("normal(3.1,2)", class = "b", coef = "sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "typeM:sc.norm_deg_sp_1"),
                        set_prior("normal(0,2)", class = "b", coef = "typeM:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "typeM:sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "sd"))


## RUN MODEL (~45 MIN.) ----
interaction.test <- brm(interaction.formula,
                        data = network.df, 
                        family = bernoulli(),
                        prior = interaction.priors,
                        algorithm = "sampling",
                        chains = 4,
                        control = list(adapt_delta = 0.9),
                        save_model = "model_output/interaction_test")

# plot(interaction.test) # everything looks good
summary(interaction.test)

## SAVE MODEL SUMMARY INFORMATION ----

# fixed effects
interaction_fixef <- round(as.data.frame(summary(interaction.test)[[17]]),2)
interaction_fixef$coefficients <- rownames(interaction_fixef)
write_csv(select(interaction_fixef, coefficients, Estimate:Rhat), "model_output/interaction_summary_fixef.csv")

# random effect summary
plyr::ldply(summary(interaction.test)[[20]]) %>%
  write_csv(., "model_output/interaction_summary_ranef.csv")

# individual random effects
interaction_ranef_sp1 <- round(as.data.frame(ranef(interaction.test)$species_1),2)
interaction_ranef_sp1$species_1 <- rownames(interaction_ranef_sp1)
write_csv(interaction_ranef_sp1, "model_output/interaction_ranef_sp1.csv")

interaction_ranef_sp2 <- round(as.data.frame(ranef(interaction.test)$species_2),2)
interaction_ranef_sp2$species_2 <- rownames(interaction_ranef_sp2)
write_csv(interaction_ranef_sp2, "model_output/interaction_ranef_sp2.csv")

interaction_ranef_network_id <- round(as.data.frame(ranef(interaction.test)$network_id),2)
interaction_ranef_network_id$network_id <- rownames(interaction_ranef_network_id)
write_csv(interaction_ranef_network_id, "model_output/interaction_ranef_network_id.csv")

interaction_ranef_subtype <- round(as.data.frame(ranef(interaction.test)$subtype),2)
interaction_ranef_subtype$subtype <- rownames(interaction_ranef_subtype)
write_csv(interaction_ranef_subtype, "model_output/interaction_ranef_subtype.csv")

interaction_ranef_id_pair <- round(as.data.frame(ranef(interaction.test)$id_pair),2)
interaction_ranef_id_pair$id_pair <- rownames(interaction_ranef_id_pair)
write_csv(interaction_ranef_id_pair, "model_output/interaction_ranef_id_pair.csv")

## PLOT AND SAVE MODEL RESULTS ----

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
} # taken from https://github.com/paul-buerkner/brms/blob/master/R/plots.R

type_main <- plot(marginal_effects(interaction.test, effects = "type"))[[1]] +
  ggplot2::xlab("Interaction Type") +
  ggplot2::ylab("Probability of Interaction")
ggsave("plots/type_main.pdf")

norm_sp1 <- plot(marginal_effects(interaction.test, effects = "sc.norm_deg_sp_1"))[[1]] +
  ggplot2::xlab("Resource Normalized Degree") +
  ggplot2::ylab("Probability of Interaction")
ggsave("plots/norm_sp1.pdf")

norm_sp2 <- plot(marginal_effects(interaction.test, effects = "sc.norm_deg_sp_2"))[[1]] +
  ggplot2::xlab("Consumer Normalized Degree") +
  ggplot2::ylab("Probability of Interaction")
ggsave("plots/norm_sp2.pdf")

type.norm_sp1 <- plot(marginal_effects(interaction.test, effects = "type:sc.norm_deg_sp_1"))[[1]] +
  ggplot2::xlab("Interaction Type") +
  ggplot2::ylab("Probability of Interaction")
ggsave("plots/type.norm_sp1.pdf")

type.norm_sp2 <- plot(marginal_effects(interaction.test, effects = "type:sc.norm_deg_sp_2"))[[1]] +
  ggplot2::xlab("Interaction Type") +
  ggplot2::ylab("Probability of Interaction")
ggsave("plots/type.norm_sp2.pdf")

norm_sp1.norm_sp2 <- plot(marginal_effects(interaction.test, effects = "sc.norm_deg_sp_1:sc.norm_deg_sp_2", surface = TRUE), stype = "raster")[[1]] +
  ggplot2::xlab("Resource Normalized Degree (scaled)") +
  ggplot2::ylab("Consumer Normalized Degree (scaled)") +
  ggplot2::scale_fill_gradientn(colors = viridis6(), name = "Probability of\n Interaction")
ggsave("plots/norm_sp1.norm_sp2.pdf")

heatmap.data <- marginal_effects(interaction.test, 
                                 effects = c("sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                                 conditions = data.frame(type = c("M","A"), 
                                                         row.names = c("Mutualistic","Antagonistic")), 
                                 surface = TRUE)
write_csv(heatmap.data[[1]], "model_output/heatmap_data.csv") # export so I can tailor these plots later since we are more interested in these interactions

source("plot_heatmaps.R") # source code and plot 3-way interactions.

## EXAMINE ROBUSTNESS OF MODEL WITHOUT PRIORS ----

# take a really long time to run
interaction.test.nopriors <- brm(interaction.formula,
                        data = network.df, 
                        family = bernoulli(),
                        #prior = interaction.priors,
                        algorithm = "sampling",
                        chains = 4,
                        
                        # adjusted control parameters in order for model to converge.
                        control = list(adapt_delta = 0.999,
                                       max_treedepth = 20),
                        save_model = "model_output/no_priors/interaction_test_nopriors")

## Examine and save all summary output

# plot(interaction.test.nopriors) # everything looks good
summary(interaction.test.nopriors)

# fixed effects
noprior_fixef <- round(as.data.frame(summary(interaction.test.nopriors)[[17]]),2)
noprior_fixef$coefficients <- rownames(noprior_fixef)
write_csv(select(noprior_fixef, coefficients, Estimate:Rhat), "model_output/no_priors/nopriors_summary_fixef.csv")

# random effect summary
plyr::ldply(summary(interaction.test.nopriors)[[20]]) %>%
  write_csv(., "model_output/no_priors/nopriors_summary_randeff.csv")

# individual random effects
noprior_ranef_sp1 <- round(as.data.frame(ranef(interaction.test.nopriors)$species_1),2)
noprior_ranef_sp1$species_1 <- rownames(noprior_ranef_sp1)
write_csv(noprior_ranef_sp1, "model_output/no_priors/nopriors_ranef_sp1.csv")

noprior_ranef_sp2 <- round(as.data.frame(ranef(interaction.test.nopriors)$species_2),2)
noprior_ranef_sp2$species_2 <- rownames(noprior_ranef_sp2)
write_csv(noprior_ranef_sp2, "model_output/no_priors/nopriors_ranef_sp2.csv")

noprior_ranef_network_id <- round(as.data.frame(ranef(interaction.test.nopriors)$network_id),2)
noprior_ranef_network_id$network_id <- rownames(noprior_ranef_network_id)
write_csv(noprior_ranef_network_id, "model_output/no_priors/nopriors_ranef_network_id.csv")

noprior_ranef_subtype <- round(as.data.frame(ranef(interaction.test.nopriors)$subtype),2)
noprior_ranef_subtype$subtype <- rownames(noprior_ranef_subtype)
write_csv(noprior_ranef_subtype, "model_output/no_priors/nopriors_ranef_subtype.csv")

noprior_ranef_id_pair <- round(as.data.frame(ranef(interaction.test.nopriors)$id_pair),2)
noprior_ranef_id_pair$id_pair <- rownames(noprior_ranef_id_pair)
write_csv(noprior_ranef_id_pair, "model_output/no_priors/nopriors_ranef_id_pair.csv")
