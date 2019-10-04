## Address Reviewer 2's comments

## Prep ----

# copied text from reproduce_analyses.Rmd

set.seed(34) # make sure simulations give same output

# load required libraries
library(tufte)
library(tidyverse)
library(brms)
library(cowplot)
library(broom)
library(kableExtra)
library(gmodels)        # make custom contrasts
library(codingMatrices) # visualize contrast matrix

# invalidate cache when the tufte version changes
knitr::opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tufte'))
options(htmltools.dir.version = FALSE)

# load and manage data
full.df <- read_csv("data/dataset.csv") %>%
  rename(resource_sp = species_1,
         consumer_sp = species_2,
         r_ND = normalized_degree_sp_1,
         c_ND = normalized_degree_sp_2) %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         subtype = as.factor(subtype)) 

weighted_subset.df <- read_csv("data/weighted_subset.csv") %>%
  rename(resource_sp = species_1,
         consumer_sp = species_2,
         r_ND = normalized_degree_sp_1,
         c_ND = normalized_degree_sp_2) %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         subtype = as.factor(subtype), 
         log.sum_r = log(sum_1),
         log.sum_c = log(sum_2))

# set the number of cores for model
options(mc.cores=parallel::detectCores ()) 

# scale data
full.df <- full.df %>%
  mutate(sc.r_ND = scale(r_ND),
         sc.c_ND = scale(c_ND))

weighted_subset.df <- weighted_subset.df %>%
  mutate(sc.r_ND = scale(r_ND),
         sc.c_ND = scale(c_ND),
         sc.log.sum_r = scale(log.sum_r),
         sc.log.sum_c = scale(log.sum_c))

# create contrast matrix
type_matrix <- rbind("_M.vs.A" = c(-1, 1))

# make contrasts for R
type_contrast <- make.contrasts(type_matrix)

# clarify description of contrasts
show_type_contrast <- mean_contrasts(type_contrast)

# add column names for clarity
attr(show_type_contrast, "dimnames")[[2]] <- levels(full.df$type)

# show contrasts
show_type_contrast %>%
  kable(., caption = "Contrasts for \\textbf{type}.", booktabs = T) %>%
  kable_styling(latex_options = c("hold_position"))

# set contrasts for both datasets
contrasts(full.df$type) <- type_contrast
contrasts(weighted_subset.df$type) <- type_contrast

# test effect of subtype within each type
full.df <- full.df %>%
  mutate(subtype_Herb.vs.Para = ifelse(type == "M", 0,
                                       ifelse(subtype == "PlantHerbivore", 1/2, -1/2)),
         subtype_Poll.vs.Disp = ifelse(type == "A", 0,
                                       ifelse(subtype == "PlantPollinator", 1/2, -1/2)))

# create contrast matrix
subset.type_matrix <- rbind("M.vs.A" = c(-1, 1/2, 1/2),
                            "Poll.vs.Disp" = c(0, 1, -1))
# "Herb.vs.Para" = c(-1, 1, 0, 0)

# make contrasts for R
subset.type_contrast <- make.contrasts(subset.type_matrix)

# clarify description of contrasts
show_subset.type_contrast <- mean_contrasts(subset.type_contrast)

# add column names for clarity
attr(show_subset.type_contrast, "dimnames")[[2]] <- levels(weighted_subset.df$subtype) # original "type"

# show contrasts
show_subset.type_contrast

# set contrasts
contrasts(weighted_subset.df$subtype) <- subset.type_contrast # changed from "type"

weighted_subset.df <- weighted_subset.df %>%
  mutate(subtype_Poll.vs.Disp = ifelse(type == "A", 0,
                                       ifelse(subtype == "PlantPollinator", 1/2, -1/2)))


## Do antagonistic networks have on average fewer total interactions tallied than mutualistic? ----

# group data at network level
weighted_subset.df.network_level <- weighted_subset.df %>%
  group_by(network_id, type, subtype) %>%
  # we look at both the qualitative (connected) and quantitative (sum_shared) interaction counts
  summarise_at(vars(connected, sum_shared), list(sum)) %>%
  ungroup()

# number of unique interactions is higher in antagonistic vs mutualistic networks
ggplot(weighted_subset.df.network_level, aes(x = type, y = connected)) +
  geom_boxplot()

ggplot(weighted_subset.df.network_level, aes(x = subtype, y = connected)) +
  geom_boxplot()

# number of total interactions sampled is higher in antagonistic vs. mutualistic networks
ggplot(weighted_subset.df.network_level, aes(x = type, y = log(sum_shared))) +
  geom_boxplot()

ggplot(weighted_subset.df.network_level, aes(x = subtype, y = log(sum_shared))) +
  geom_boxplot()

# mutualistic networks have lower interaction counts (qualitative) compared to antagonistic ones.
summary(glm(connected ~ type, data = weighted_subset.df.network_level, family = quasipoisson(link = "log")))

# mutualistic networks have lower interaction counts (quantitative) compared to antagonistic ones,
# but not significantly so, which is surprising
summary(glm(sum_shared ~ type, data = weighted_subset.df.network_level, family = quasipoisson(link = "log")))

## Rarefaction ----

library(lme4) # using as a substitute for the Bayesian models

# explore distribution of observed interactions
summary(weighted_subset.df$sum_shared) # huge variation in number of observed interactions
hist(log1p(weighted_subset.df$sum_shared)) # log distribution is still quite skewed

## Make really long dataset where each row corresponds to a sample interaction

# focus dataset
weighted_long <- select(weighted_subset.df, network_id, type, subtype, resource_sp, consumer_sp, connected, sum_shared)

# for each interaction, replicate it as many times as it was observed
weighted_long.list <- list()
for(i in 1:nrow(weighted_long)){ 
  if(weighted_long$sum_shared[i] > 0){
    weighted_long.list[[i]] <- slice(weighted_long[i, ], rep(1:n(), each = weighted_long$sum_shared[i]))
  } else {
    weighted_long.list[[i]] <- weighted_long[i, ]
  }
}
weighted_long.df <- plyr::ldply(weighted_long.list) %>% select(-sum_shared) # tidy into a dataframe
head(weighted_long.df)

# now I need to randomly subsample the antagonistic networks to the mutualistic ones
table(distinct(weighted_long.df, network_id, subtype)$subtype) 

n_sims <- 1000

three_way.FE <- c() # vector of fixed-effect estimate for three-way interaction term
three_way.P <- c() # vector of P-values for the three-way term
for(i in 1:n_sims){
  # subsample Antagonistic networks to match the number of mutualistic ones
  subsample_A_networks <- distinct(filter(weighted_long.df, type == "A"), network_id) %>% sample_n(size = 19, replace = F) 
  
  M_networks <- distinct(filter(weighted_long.df, type == "M"), network_id) # no need to subsample because there are fewer
  
  # randomly pair them together (note randomization was already done with subsample_A_networks)
  random_pairs <- data.frame(network_pairs = LETTERS[1:19],
                             M = as.character(M_networks$network_id),
                             A = as.character(subsample_A_networks$network_id))
  
  # create data frame to determine subsample size for each network
  subsample_df <- left_join(random_pairs, select(weighted_subset.df.network_level, M = network_id, M_sum = sum_shared)) %>%
    left_join(., select(weighted_subset.df.network_level, A = network_id, A_sum = sum_shared)) %>%
    # subsample is based on smallest number of observed interactions in a network
    mutate(subsample_size = ifelse(M_sum < A_sum, M_sum, A_sum)) %>%
    select(network_pairs, M, A, subsample_size) %>%
    # organize for joining to other datasets
    gather(type, network_id, -subsample_size, -network_pairs)
  
  # now I need to make the new data frame based on the random subsample
  new_data <- weighted_long.df %>%
    # filter to randomly sampled networks
    filter(network_id %in% c(as.character(M_networks$network_id), as.character(subsample_A_networks$network_id))) %>%
    # following example from: https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
    # to get different sample sizes for each group
    group_by(network_id) %>%
    nest() %>%
    left_join(., subsample_df) %>%
    mutate(samp = map2(data, subsample_size, sample_n)) %>%
    select(-data) %>%
    unnest(samp) %>%
    # now that I've subsampled, I can reduce the dataset to the normal size for analysis
    # i.e. only 1 interaction per pair per network
    distinct %>%
    # add back in some data for the analysis
    # note that sc.r_ND and sc.c_ND are calculated from the original data and not the subsampled interactions
    # I think this is okay, because their normalized degree is based on the entire interaction network
    # rather than their partner fidelity across networks
    left_join(., select(weighted_subset.df, network_id, resource_sp, consumer_sp, id_pair, sc.r_ND, sc.c_ND))
  
  # run test with glmer
  new_glmer <- glmer(connected ~ subtype + type*sc.r_ND*sc.c_ND + (1|consumer_sp) + (1|resource_sp) + (1|id_pair) + (1|network_id), 
                     data = new_data, family = binomial(link = "logit"), 
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  
  # get model effects for 3-way interaction term
  new_glmer_effects <- broom::tidy(new_glmer, effects = "fixed")
  three_way.FE[i] <- filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$estimate
  three_way.P[i] <- filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$p.value
  print(paste(i/n_sims*100,"% done!")) # track progress
}

# write data from simulation 
write_csv(x = data.frame(FixedEffect = three_way.FE, P = three_way.P), path = "data/three_way_permutation_data.csv")

# reference statistic based on three-way term in Table 7 of supplementary material
ref.FE <- -0.68

# looks to me like equalizing the sample size doesn't alter the distribution of the test statistic,
# suggesting that our analysis is robust
ggplot(data.frame(FE = three_way.FE), aes(x = FE)) +
  geom_histogram(alpha = 0.5, color = "grey") +
  geom_vline(xintercept = ref.FE, linetype = "dotted") +
  xlab("Coefficient estimate") +
  ylab("Count")

# proportion of times that equalizing sample size reduces effect size
sum(three_way.FE > ref.FE) / length(three_way.FE) # 55% of the time.
# i.e. since the reference statistic is not larger than 90% or more of the subsampled
# distribution values, then sampling intensity does not matter for our results.

# proportion of times that equalizing sample size flips the sign of the effect
sum(three_way.FE > 0) / length(three_way.FE) # 2%

# proportion of times the 3-way interaction is statistically significant
sum(three_way.P < 0.05) / length(three_way.P) # 35%

