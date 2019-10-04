## Rarefaction

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

#### 10.4.Fri ----
summary(weighted_subset.df$sum_shared) # huge variation in number of observed interactions
hist(log1p(weighted_subset.df$sum_shared)) # log distribution is still quite skewed

# make really long with number of observed interactions
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

subsample_A_networks <- distinct(filter(weighted_long.df, type == "A"), network_id) %>% sample_n(size = 19, replace = F) 
subsample_A_networks

#assign_pollinator <- subsample_A_networks$network_id[1:10]
#assign_disperser <- subsample_A_networks$network_id[11:19]

#pollinator_networks <- distinct(filter(weighted_long.df, type == "M", subtype == "PlantPollinator"), network_id)
#disperser_networks <- distinct(filter(weighted_long.df, type == "M", subtype == "PlantSeedDisperser"), network_id)

M_networks <- distinct(filter(weighted_long.df, type == "M"), network_id) # no need to subsample because there are fewer
# also the randomization is already done with the subsample_A_networks

random_pairs <- data.frame(network_pairs = LETTERS[1:19],
           M = as.character(M_networks$network_id),
           A = as.character(subsample_A_networks$network_id))

# for all but one occassion, the mutualistic networks have fewer observed interactions
subsample_df <- left_join(random_pairs, select(weighted_subset.df.network_level, M = network_id, M_sum = sum_shared)) %>%
  left_join(., select(weighted_subset.df.network_level, A = network_id, A_sum = sum_shared)) %>%
  mutate(subsample_size = ifelse(M_sum < A_sum, M_sum, A_sum)) %>%
  select(network_pairs, M, A, subsample_size) %>%
  gather(type, network_id, -subsample_size, -network_pairs)

# now I need to make the new data frame based on the random subsample
new_data <- weighted_long.df %>%
  filter(network_id %in% c(as.character(M_networks$network_id), as.character(subsample_A_networks$network_id))) %>%
  # distinct(network_id) # checked and it worked
  group_by(network_id) %>%
  nest() %>%
  left_join(., subsample_df) %>%
  # following example from: https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
  mutate(samp = map2(data, subsample_size, sample_n)) %>%
  select(-data) %>%
  unnest(samp) %>%
  # now that I've subsampled, I can reduce the dataset to the normal size for analysis
  # i.e. only 1 interaction per pair per network
  distinct %>%
  # add back in some data for the analysis
  left_join(., select(weighted_subset.df, network_id, resource_sp, consumer_sp, id_pair, sc.r_ND, sc.c_ND))
  
# run test with glmer
library(lme4)
new_glmer <- glmer(connected ~ subtype + type*sc.r_ND*sc.c_ND + (1|consumer_sp) + (1|resource_sp) + (1|id_pair) + (1|network_id), 
                    data = new_data, family = binomial(link = "logit"), 
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

# get model effects for 3-way interaction term
new_glmer_effects <- broom::tidy(new_glmer, effects = "fixed")
filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$estimate
filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$p.value
  
#### Address Reviewer #2 comments ----

# Do antagonistic networks have on average fewer total interactions tallied than mutualistic?
weighted_subset.df.network_level <- weighted_subset.df %>%
  group_by(network_id, type, subtype) %>%
  # we look at both the qualitative (connected) and quantitative (sum_shared) interaction counts
  summarise_at(vars(connected, sum_shared), list(sum)) %>%
  ungroup()

left_join(subsample_A_networks, weighted_subset.df.network_level) %>% left_join(., select(random_pairs, network_pairs, network_id = A))
left_join(M_networks, weighted_subset.df.network_level) %>% left_join(., select(random_pairs, network_pairs, network_id = M))

# NA values represent pairs that did not co-occur in more than 1 network
# 0 indicates no interaction in this network
# 1 indicates an interaction
filter(weighted_subset.df, network_id == "A_HP_001") %>%
  select(resource_sp, consumer_sp, connected) %>%
  spread(consumer_sp, connected)

# as above, except now the value corresponds to the number of observed interactions
filter(weighted_subset.df, network_id == "A_HP_001") %>%
  select(resource_sp, consumer_sp, sum_shared) %>%
  spread(consumer_sp, sum_shared)

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

## Get reference statistic
# didn't converge # ref.glmer <- lme4::glmer(connected ~ subtype + type*sc.r_ND*sc.c_ND + (1|consumer_sp) + (1|resource_sp) + (1|id_pair) + (1|network_id), data = weighted_subset.df, family = binomial(link = "logit"))
table(distinct(weighted_subset.df, subtype, network_id)$subtype)

# max number of interactions sampled in mutualistic network
max_M_sum_shared <- max(filter(weighted_subset.df.network_level, type == "M")$sum_shared) 

# 19 networks 
filter(weighted_subset.df.network_level, type == "A", sum_shared < max_M_sum_shared)

subsample_A_networks <- filter(weighted_subset.df.network_level, type == "A") %>% sample_n(size = 19, replace = F) %>% mutate(network_id = as.character(network_id)) %>% .$network_id

assign_pollinator <- subsample_A_networks[1:10]
assign_disperser <- subsample_A_networks[11:19]

pollinator_networks <- as.character(filter(weighted_subset.df.network_level, type == "M", subtype == "PlantPollinator")$network_id)
disperser_networks <- as.character(filter(weighted_subset.df.network_level, type == "M", subtype == "PlantSeedDisperser")$network_id)

data.frame(network_pairs = as.character(1:19),
           M = c(pollinator_networks, disperser_networks),
           A = c(assign_pollinator, assign_disperser))


weighted_subset.df


pollinator_networks
## 
library(weboflife)
pollination_networks <- get_networks(interaction_type = "pollination", data_type = "weighted")

test <- pollination_networks[pollinator_networks]

test1 <- test[[1]] %>%
  gather(Consumer, value = InteractionCount, -V1) %>%
  rename(Resource = V1) %>%
  filter(InteractionCount > 0) %>%
  unite(col = "ConsumerResource", Consumer, Resource)

apply(test1, 1, function(x) rep(x[,1], x[,2]))

rep(test1[1,1], test1[1,2]) # this is what I want

## Rarefaction analysis code chunk ----
min_sample <- min(rowSums(select(df, euro:platy))) # get minimum number of samples

# get individual samples for each site
individual_samples <- apply(select(df, euro:platy), 1, function(x) rep(names(x), x))

# convert to a list with site information
samples_list <- list()
for(i in 1:20){
  samples_list[[i]] <- data.frame(site = i, species = individual_samples[[i]])
}
samples_df <- plyr::ldply(samples_list) # convert to data frame

# perform rarefaction, keeping track of species identities
n_sims <- 10000 # number of rarefied samples for each site

sim_data <- list()
for(i in 1:n_sims){
  # for each site, randomly sample 'min_sample' from each site (without replacement) and store in a data frame
  sim_data[[i]] <- group_by(samples_df, site) %>% sample_n(., size = min_sample, replace = FALSE) %>% mutate(sim = i)
}
sim_df <- plyr::ldply(sim_data) 

# aggregate simulation data
sum_sims <- sim_df %>%
  group_by(site, species, sim) %>%
  summarise(abund = n()) %>%
  arrange(site, sim, species) %>%
  as.data.frame()

spread_sims <- spread(sum_sims, species, abund, fill=0) %>%
  # check that total samples for each site is equal to 'min_sample'
  mutate(check = rowSums(select(., -site, -sim)) - min_sample) # all should be zero
summary(spread_sims$check) # and they are
