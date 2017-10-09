## LOAD REQUIRED LIBRARIES ----
library(tidyverse)
library(brms)

## LOAD AND MANAGE THE DATA ----
network.df <- read_csv("t_glmm_plus_connectance.csv") %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         c.norm_deg_sp_1 = normalized_degree_sp_1 - mean(normalized_degree_sp_1),
         c.norm_deg_sp_2 = normalized_degree_sp_2 - mean(normalized_degree_sp_2),
         sc.norm_deg_sp_1 = scale(normalized_degree_sp_1),
         sc.norm_deg_sp_2 = scale(normalized_degree_sp_2),
         expected_probability = normalized_degree_sp_1*normalized_degree_sp_2, ## DOESN'T MAKE SENSE, NEEDS TO BE (DEGREE + DEGREE) / 2 WHICH CORRESPONDS TO JORDI'S NULL MODEL. OR WE DON'T USE THIS MEASURE AT ALL...
         c.expected_probability = scale(expected_probability, center = TRUE, scale = FALSE),
         asymmetry = abs(normalized_degree_sp_1 - normalized_degree_sp_2),
         c.asymmetry = scale(asymmetry, center = TRUE, scale = FALSE))#,
         #c.connectance = network_connectance - mean(network_connectance),
         #delta.norm_deg = normalized_degree_sp_1 - normalized_degree_sp_2)
network.df %>% group_by(id_pair, type, subtype) %>% summarise(interactions = sum(connected), cooccur = n()) %>% lm(interactions/cooccur ~ type + subtype, data=.) %>% summary()#ggplot(., aes(x = type, y = interactions/cooccur, color = subtype, fill = subtype)) + geom_boxplot()

network.df$resid_deg_sp_1  <- residuals(lm(normalized_degree_sp_1 ~ network_connectance, data = network.df))
network.df$resid_deg_sp_2  <- residuals(lm(normalized_degree_sp_2 ~ network_connectance, data = network.df))
network.df$type.SumContr <- network.df$type 
contrasts(network.df$type.SumContr) <- contr.sum

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

test <- glm(connected ~ type*sc.norm_deg_sp_1*sc.norm_deg_sp_2, data = network.df, family = "binomial")
summary(test)
plot(effect("type.SumContr", test))
plot(effect("type:sc.norm_deg_sp_1", test))
plot(effect("sc.norm_deg_sp_1*sc.norm_deg_sp_2", test))

network.df$cut.sp_2 <- ifelse(network.df$sc.norm_deg_sp_2 > 0, "Generalist", "Specialist")
ggplot(network.df, aes(x = sc.norm_deg_sp_1, y = connected, color = cut.sp_2)) + binomial_smooth()

ggplot(network.df, aes(x = expected_probability, y = connected, color = type)) + binomial_smooth()
ggplot(network.df, aes(x = asymmetry, y = connected, color = type)) + binomial_smooth()

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

## exploratory plots
ggplot(network.df, aes(x = log(number_sp_label_1 + number_sp_label_2), y = network_connectance)) +
  geom_point()
ggplot(network.df, aes(x = number_sp_label_1 + number_sp_label_2, y = normalized_degree_sp_1)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(network.df, aes(x = network_connectance, y = log(degree_1))) +
  geom_point() + geom_smooth(method = "lm") 
ggplot(network.df, aes(x = network_connectance, y = normalized_degree_sp_1)) +
  geom_point() + geom_smooth(method = "lm") +
  geom_hline(yintercept = mean.norm_deg_sp_1)
ggplot(network.df, aes(x = network_connectance, y = normalized_degree_sp_2)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(network.df, aes(x = network_connectance, y = resid_deg_sp_1)) +
  geom_point()# + geom_smooth(method = "lm")
ggplot(network.df, aes(x = network_connectance, y = resid_deg_sp_2)) +
  geom_point()# + geom_smooth(method = "lm")

ggplot(network.df, aes(x = network_connectance, y = connected, color = type)) +
  geom_point() + 
  binomial_smooth() 
ggplot(network.df, aes(x = normalized_degree_sp_1, y = connected)) +
  geom_point() + 
  binomial_smooth() 
ggplot(network.df, aes(x = normalized_degree_sp_2, y = connected)) +
  geom_point() + 
  binomial_smooth() 

ggplot(network.df %>% mutate(connected = as.factor(connected)), 
       aes(x = normalized_degree_sp_2, y = normalized_degree_sp_1, fill = connected, color = connected)) +
  geom_point(shape = 21, size = 1) + 
  facet_wrap(~type) + 
  scale_fill_manual(values = c("white", "black")) + 
  scale_color_manual(values = c("grey", "black")) + 
  theme_bw()

# mutualists
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


# antagonists
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
  facet_wrap(~type, scales = "free") + # facet_wrap()
  #scale_fill_manual(values = c("white", "black")) + 
  #scale_color_manual(values = c("grey", "black")) + 
  theme_bw()

# from species 2 perspective
ggplot(network.df %>% mutate(connected = connected, 
                             sp_1_type = factor(ifelse(normalized_degree_sp_1 > mean.norm_deg_sp_1, "sp_1_generalist", "sp_1_specialist"), levels = c("sp_1_generalist", "sp_1_specialist")),
                             sp_2_type = factor(ifelse(normalized_degree_sp_2 > mean.norm_deg_sp_2, "sp_2_generalist", "sp_2_specialist"), levels = c("sp_2_specialist", "sp_2_generalist"))),
       aes(x = normalized_degree_sp_2, y = connected, fill = sp_1_type, color = sp_1_type)) +
  geom_point(shape = 21, size = 4) + 
  binomial_smooth() +
  facet_wrap(~type, scales = "free") + # facet_wrap()
  #scale_fill_manual(values = c("white", "black")) + 
  #scale_color_manual(values = c("grey", "black")) + 
  theme_bw()

## BAYESIAN MODEL OF INTERACTION PROBABILITY
interaction.formula <- brmsformula(connected ~                                # probability of species interaction
                                     type*c.norm_deg_sp_1*c.norm_deg_sp_2 +   # 3-way interaction between type and normalized degree of each species
                                     (1 | network_id) +                       # allow probability of interaction to vary among networks
                                     (1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                                     (1 | species_1) +                        # allow probability of interaction to vary among resource species
                                     (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                                     (1 | id_pair))                           # allow probability of interaction to vary among unique species interactions
                                     #(1 | label_1) +
                                     #(1 | label_2))

## Get intuition for estimating priors
# playing around with different hypothetical probabilities suggests to me that the logistic regression coefficient will be 
# between -5 and 5. Therefore, a normal prior with mean 0 and sd = 2, seems reasonable to me.

# Let's take a hypothetical example, where the effect of interaction type is really large, ranging from
# a probability of 0.25 to 0.75. This corresponds to a logistic regression coefficient of ~2.2
pM <- 0.99 # probability of a mutualistic interaction when two species co-occur 
pA <- 0.01 # probability of an antagonistic interaction when two species co-occur 
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

# I think a normal prior with mean = 0 and sd = 1 would be appropriate.

get_prior(interaction.formula, data = network.df, family = bernoulli()) # see which priors I need to set.

## RATIONALE FOR MODEL PRIORS

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
pSpec_sp1 <- 0.01 #mean(network.df$normalized_degree_sp_1) - sd(network.df$normalized_degree_sp_1) # probability of a relatively specialized species interacting
pGen_sp1 <- 0.99 #mean(network.df$normalized_degree_sp_1) + sd(network.df$normalized_degree_sp_1) # probability of a relatively generalized species interacting

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
# deviation of 1 is reasonable since it covers the range from 0 to 5 (and I never expect this coefficient to be negative). It will
# also make the data "work" for larger values.

# For interactions involving interaction type, I'm going to specify the same regularizing prior
# as for "type", because I could see interaction type having positive or negative effects on this relationship
# and I want to make the data "work" for those strong effects.
# For the interaction with normalized degree, I expect this relationship to always be positive, but I don't know how much
# so I'm going to set a lower bound of zero, but I'll increase the standard deviation to 2 to allow more flexibility
# in having positive values.


# set the priors
# WHY CAN'T I SET A LOWER BOUND? 
interaction.priors <- c(set_prior("normal(0,1)", class = "b", coef = "type.SumContr1"),
                        set_prior("normal(3.4,1)", class = "b", coef = "sc.norm_deg_sp_1"),
                        set_prior("normal(3.1,1)", class = "b", coef = "sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "b", coef = "sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,1)", class = "b", coef = "type.SumContr1:sc.norm_deg_sp_1"),
                        set_prior("normal(0,1)", class = "b", coef = "type.SumContr1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,1)", class = "b", coef = "type.SumContr1:sc.norm_deg_sp_1:sc.norm_deg_sp_2"),
                        set_prior("normal(0,2)", class = "sd"))


## run the model.
sub <- network.df %>% filter(expected_probability > 0.5) # about 0.7 for normalized degree of each species.
table(sub$connected)

interaction.test <- brm(connected ~                                # probability of species interaction
                             c.expected_probability + type*c.asymmetry +   # 3-way interaction between type and normalized degree of each species
                             (1 | network_id) +                       # allow probability of interaction to vary among networks
                             (1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                             (1 | species_1) +                        # allow probability of interaction to vary among resource species
                             (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                             (1 | id_pair),
                           data = network.df, # filter(network.df, expected_probability < 0.5),
                           family = bernoulli(),
                           #prior = interaction.priors,
                           algorithm = "meanfield")
summary(interaction.test)
plot(marginal_effects(interaction.brm))

intercept <- -0.12
typeM <- 1.15

pA <- exp(intercept)/(1+exp(intercept)) # 0.47
pM <- exp(intercept + typeM)/(1 + exp(intercept + typeM)) # 0.74
log(pM / (1-pM)) - log(pA / (1-pA))

plot(marginal_effects(interaction.test, 
                      effects = c("c.expected_probability:c.asymmetry"),
                      conditions = data.frame(type = c("M","A"), 
                                              row.names = c("Mutualistic","Antagonistic"))))


# having trouble with "sampling" and current priors...
interaction.brm <- brm(interaction.formula,
                       data = network.df,
                       family = bernoulli(),
                       #prior = interaction.priors,
                       algorithm = "meanfield") #"sampling",
                       #control = list(adapt_delta = 0.99))

network.df <- network.df %>% mutate(expected.probability = normalized_degree_sp_1*normalized_degree_sp_2)
interaction.alt.brm <- brm(connected ~                                # probability of species interaction
                             type*expected.probability +   # 3-way interaction between type and normalized degree of each species
                             (1 | network_id) +                       # allow probability of interaction to vary among networks
                             (1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                             (1 | species_1) +                        # allow probability of interaction to vary among resource species
                             (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                             (1 | id_pair),
                       data = network.df,
                       family = bernoulli(),
                       #prior = interaction.priors,
                       algorithm = "meanfield")
summary(interaction.alt.brm)

launch_shiny(interaction.brm)
# analyze model output
plot(interaction.brm)
summary(interaction.brm)
plot(marginal_effects(interaction.brm))

plot(marginal_effects(interaction.brm, effects = "type"))

plot(marginal_effects(interaction.brm, effects = "type:c.norm_deg_sp_1"))

plot(marginal_effects(interaction.brm, effects = "type:c.norm_deg_sp_2"))

plot(marginal_effects(interaction.brm, effects = "c.norm_deg_sp_1:c.norm_deg_sp_2"))

# increasing the normalized degree of species 1 increases the probability of an interaction.
# this is also the case for species 2; however, mutualism seems to increase the probability
# of their being a connection when species 2 has a low normalized degree.
plot(marginal_effects(interaction.brm, 
                      effects = c("c.norm_deg_sp_1:c.norm_deg_sp_2"),
                      conditions = data.frame(type = c("M","A"), 
                                              row.names = c("Mutualistic","Antagonistic"))))

##
summary(network.df$sc.norm_deg_sp_1)
interaction.newdata <- expand.grid(c.norm_deg_sp_1 = seq(min(network.df$c.norm_deg_sp_1), max(network.df$c.norm_deg_sp_1), by = 0.01),
                                   c.norm_deg_sp_2 = seq(min(network.df$c.norm_deg_sp_2), max(network.df$c.norm_deg_sp_2), by = 0.01),
                                   type = c("M","A"),
                                   network_id = NA,
                                   id_pair = NA,
                                   species_1 = NA,
                                   species_2 = NA,
                                   subtype = NA)
interaction.fits <- fitted(interaction.brm, newdata = interaction.newdata)
plot.data <- cbind.data.frame(interaction.newdata, interaction.fits, rbind(expect, expect)) %>%
  select(c.norm_deg_sp_1, c.norm_deg_sp_2, norm_deg_sp_1, norm_deg_sp_2, type, Estimate, lowerCI = `2.5%ile`, upperCI = `97.5%ile`) %>%
  mutate(uncertainty = upperCI - lowerCI,
         expected.probability = norm_deg_sp_1*norm_deg_sp_2,
         deviation = Estimate - expected.probability)
expect <- data.frame(expand.grid(norm_deg_sp_1 = seq(min(network.df$normalized_degree_sp_1), max(network.df$normalized_degree_sp_1), by = 0.01),
                                   norm_deg_sp_2 = seq(min(network.df$normalized_degree_sp_2), max(network.df$normalized_degree_sp_2), by = 0.01)))
ggplot(expand.grid(norm_deg_sp_1 = seq(min(network.df$normalized_degree_sp_1), max(network.df$normalized_degree_sp_1), by = 0.01),
                   norm_deg_sp_2 = seq(min(network.df$normalized_degree_sp_2), max(network.df$normalized_degree_sp_2), by = 0.01)) %>%
         mutate(expected.probability = norm_deg_sp_1*norm_deg_sp_2),
       aes(x = norm_deg_sp_1, y = norm_deg_sp_2)) +
  geom_tile(aes(fill = expected.probability)) +
  scale_fill_gradient(low = "white", high = "red", limit = c(0,1))

# default expected probability. It actually doesn't make sense to contstrain this to one, because 
# both antagonistic and mutualistic interactions have a higher probability of interacting than we would
# expect from simply their normalized degree.
ggplot(plot.data, aes(x = c.norm_deg_sp_1, y = c.norm_deg_sp_2)) + 
  geom_tile(aes(fill = Estimate)) +
  facet_wrap(~type)+
  scale_fill_gradient(low = "white", high = "red", limit = c(0,1))

ggplot(plot.data, aes(x = norm_deg_sp_1, y = norm_deg_sp_2)) + 
  geom_tile(aes(fill = deviation)) +
  facet_wrap(~type)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0)

ggplot(plot.data, aes(x = c.norm_deg_sp_1, y = c.norm_deg_sp_2)) + 
  geom_tile(aes(fill = uncertainty)) +
  facet_wrap(~type.SumContr) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)

plot.data.diff <- plot.data %>%
  select(c.norm_deg_sp_1, c.norm_deg_sp_2, type, Estimate) %>%
  spread(type, Estimate) %>%
  mutate(predicted.probability.diff = M - A,
         std.probability.diff = M/A,
         std.percent.probability.diff = (M-A)/A)
summary(plot.data.diff$predicted.probability.diff)

ggplot(plot.data.diff, aes(x = c.norm_deg_sp_1, y = c.norm_deg_sp_2)) + 
  geom_tile(aes(fill = predicted.probability.diff)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)

ggplot(plot.data.diff, aes(x = c.norm_deg_sp_1, y = c.norm_deg_sp_2)) + 
  geom_tile(aes(fill = std.probability.diff)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1) +
  scale_x_continuous(labels = c("0","0.25","0.5","0.75","1"), breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_y_continuous(labels = c("0","0.25","0.5","0.75","1"), breaks = c(-0.25, 0, 0.25, 0.5, 0.75))

#ggplot(plot.data.diff, aes(x = sc.norm_deg_sp_1, y = sc.norm_deg_sp_2)) + 
#  geom_tile(aes(fill = uncertainty)) +
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)


## SPLIT MUTALISTIC AND ANTAGONISTIC DATASETS TO ENSURE PROPER INTERPRETATION
subset.formula <- brmsformula(connected ~                                     # probability of species interaction
                                     c.norm_deg_sp_1*c.norm_deg_sp_2 +        # 2-way interaction between normalized degree of each species
                                     (1 | network_id) +                       # allow probability of interaction to vary among networks
                                     (1 | subtype) +                          # allow probability of interaction to vary among interaction subtypes.
                                     (1 | species_1) +                        # allow probability of interaction to vary among resource species
                                     (1 | species_2) +                        # allow probability of interaction to vary among consumer species
                                     (1 | id_pair))        

mutualism.brm <- brm(subset.formula,
                       data = filter(network.df, type == "M"),
                       family = bernoulli(),
                       algorithm = "meanfield")

# analyze model output
plot(mutualism.brm)
summary(mutualism.brm)

# increasing the normalized degree of species 1 increases the probability of an interaction with species 2.
# this occurs regardless of the normalized degree of species 2.
plot(marginal_effects(mutualism.brm))


antagonistic.brm <- brm(subset.formula,
                     data = filter(network.df, type == "A"),
                     family = bernoulli(),
                     algorithm = "meanfield")

# analyze model output
plot(antagonistic.brm)
summary(antagonistic.brm)

# increasing the normalized degree of species 1 increases the probability of an interaction with species 2.
# this is enhanced by the normalized degree of species 2.
plot(marginal_effects(antagonistic.brm))
