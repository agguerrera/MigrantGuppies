# Load required packages
library(dplyr) # Make sure MASS package is not loaded or it will mask dplyr functions
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(car)
library(here)

# Necessary input files:
# rgb_uv_pool5DT1PR.rds : contains Colormesh data
# pool_info.csv : contain trial- and pool-specific information for each male

#######################################
# Guppy color-pattern rarity analysis #
#######################################
# Each trial has RC1, RC2, RR, MU; repeated migrants are identifiable by ID + trial
# Trial rarity is calculated based on the Euclidean distance to the mean trial centroid
# Pool rarity is calculated based on the Euclidean distance to the mean pool centroid. 
# We exclude migrant males from the mean pool centroid as they were sourced from different pools
# The pool rarity for residents is thus a measure of how rare each resident male (RC1, RC2, RR) 
# is from the mean color pattern of all males collected in that focal pool.
# Migrant male pool-level rarity is then computed separately by adding an individual migrant's 
# color vector to the color matrix of the focal resident pool.

#WARNING: This clears your global working environment to start fresh between scripts
rm(list = ls())

colors <- readRDS(here::here("data", "rgb_uv_pool5DT1PR.rds"))
pool_info <- read.csv(here("data", "pool_info.csv")) %>% mutate(trial = as.integer(trial))

# We first make a variable that stores all RGB and UV color measurements for all males
first_color <- match("r_interior1", names(colors)) # finds where first color columm starts, excludes metadata
color_vars <- names(colors)[seq.int(first_color, ncol(colors))]

# We need to determine what resident focal pool is associated with each trial male
# Then create a dataframe containing each trial's four males, their color measurements, and the recipient pool against which their rarity will be calculated.
resident_members <- pool_info %>%
  dplyr::filter(type %in% c("RC1", "RC2", "RR")) %>%
  dplyr::left_join(
    colors %>% dplyr::select(ID, PoolStream),
    by = "ID",
    relationship = "many-to-one")

recipient_pools <- resident_members %>%
  group_by(trial) %>%
  summarise(n_residents = n(), n_pools = n_distinct(PoolStream, na.rm = TRUE),
            RecipientPool = first(PoolStream), .groups = "drop")

trial_males <- pool_info %>%
  dplyr::left_join(
    recipient_pools %>%
      dplyr::select(trial, RecipientPool),
    by = "trial",
    relationship = "many-to-one"
  ) %>%
  dplyr::left_join(
    colors %>%
      dplyr::select(ID, dplyr::all_of(color_vars)),
    by = "ID",
    relationship = "many-to-one"
  )

# Euclidean distance from a reference centroid
# The sweep function is very important! It tells R to subtract each color value from the appropriate column mean!)
centroid_distance <- function(x, center = colMeans(x)) {
  sqrt(rowSums(sweep(x, 2, center, "-")^2))
}

# Trial rarity: distance from trial centroid of the four males presented together.
trial_groups <- split(trial_males, trial_males$trial)
TrialRarity_df <- bind_rows(lapply(trial_groups, function(dat) {
  x <- as.matrix(dat[, color_vars, drop = FALSE])
  tibble(ID = dat$ID, trial = dat$trial, TrED = centroid_distance(x))
}))

# Let's look at pairwise color-pattern distances based on in the field male type assignments
comparison_types <- c("MU:RC1", "MU:RC2", "RC1:RC2", "RC1:RR", "RC2:RR", "MU:RR")
pairwise_distances <- bind_rows(lapply(trial_groups, function(dat) {
  x <- as.matrix(dat[, color_vars, drop = FALSE])
  dm <- as.matrix(dist(x))
  dimnames(dm) <- list(dat$type, dat$type)
  tibble(trial = dat$trial[1], comparison = comparison_types,
         distance = vapply(strsplit(comparison_types, ":"),
                           function(z) dm[z[1], z[2]], numeric(1)))
}))
pairwise_summary <- pairwise_distances %>%
  group_by(comparison) %>%
  summarise(mean = mean(distance), SE = sd(distance) / sqrt(n()), .groups = "drop")

print(pairwise_summary)

# Supplemental Figure S5: Pairwise color pattern distances
# Now we can plot it
pairwise_plot_data <- pairwise_distances %>%
  mutate(comparison = case_when(
    comparison == "RC1:RC2" ~ "RC:RC",
    comparison %in% c("RC1:RR", "RC2:RR") ~ "RC:RR",
    comparison %in% c("MU:RC1", "MU:RC2") ~ "RC:MU",
    comparison == "MU:RR" ~ "RR:MU")) %>%
  group_by(trial, comparison) %>%
  summarise(distance = mean(distance), .groups = "drop") %>%
  mutate(comparison = factor(comparison, levels = c("RC:RC", "RC:RR", "RC:MU", "RR:MU")))

FigS5 <- ggplot(pairwise_plot_data, aes(comparison, distance)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.5), width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color="red",position = position_dodge(width = 0.5)) + #adds mean
  labs(x = "Male type comparison", y = "Pairwise color pattern distances") + 
  theme(axis.text = element_text(size = 15, color = "black"), axis.title = element_text(size = 15)) + 
  theme_classic() 

ggsave(filename = here("figures", "FigS5.png"),
       plot = FigS5,
       width = 8,
       height = 6,
       dpi = 300)

# Pool rarity: all four trial males are compared to centroid of focal resident pool
# Migrants are not used to calculate the centroid.Their rarity is calculated separately. 
trial_pool_map <- resident_members %>%
  dplyr::group_by(trial) %>%
  dplyr::summarise(
    PoolStream = dplyr::first(PoolStream),
    .groups = "drop"
  )

focal_pools <- unique(trial_pool_map$PoolStream)

pool_centers <- lapply(
  split(
    colors[colors$PoolStream %in% focal_pools, ],
    colors$PoolStream[colors$PoolStream %in% focal_pools]
  ),
  function(dat) {
    x <- as.matrix(dat[, color_vars, drop = FALSE])
    colMeans(x)
  }
)

PoolRarity_df <- bind_rows(lapply(trial_groups, function(dat) {
  pool <- unique(dat$RecipientPool)
  x <- as.matrix(dat[, color_vars, drop = FALSE])
  tibble(ID = dat$ID, trial = dat$trial,
         Pool.ED = centroid_distance(x, pool_centers[[pool]]))
}))

# Combine raw trial and pool rarity measurements
combined_rarity <- trial_males %>%
  dplyr::select(
    dplyr::all_of(setdiff(
      names(pool_info),
      c("Trial", "TrialRarity", "PoolRarity", "TrED", "Pool.ED")
    )),
    RecipientPool
  ) %>%
  dplyr::left_join(
    TrialRarity_df %>%
      dplyr::select(ID, trial, TrED),
    by = c("ID", "trial"),
    relationship = "one-to-one"
  ) %>%
  dplyr::left_join(
    PoolRarity_df %>%
      dplyr::select(ID, trial, Pool.ED),
    by = c("ID", "trial"),
    relationship = "one-to-one"
  ) %>%
  dplyr::rename(
    Trial = trial,
    TrialRarity = TrED,
    PoolRarity = Pool.ED
  )

# Export dataset

write.csv(combined_rarity, "~/Desktop/Trinidad-Migrant-MS/AmNat/Dryad/combined_rarity.csv", row.names = FALSE)

# Comparisons across male types.
combined_rarity <- combined_rarity %>%
  mutate(male_type = factor(male_type, levels = c("common", "rare", "migrant")))


# Number of trials where mean common-resident pool rarity exceeds migrant rarity.
trial_means <- combined_rarity %>%
  group_by(Trial, male_type) %>%
  summarise(mean_pool_rarity = mean(PoolRarity), .groups = "drop") %>%
  pivot_wider(names_from = male_type, values_from = mean_pool_rarity)
cat("Trials with common > migrant pool rarity:",
    sum(trial_means$common > trial_means$migrant), "of", nrow(trial_means), "\n")

# Prepare data
rarity <- combined_rarity %>%
  dplyr::mutate(
    Trial = factor(Trial),
    male_type = factor(
      male_type,
      levels = c("common", "rare", "migrant")
    ),
    PoolRarity_scaled = as.numeric(scale(PoolRarity, scale=TRUE, center=TRUE)), 
    TrialRarity_scaled = as.numeric(scale(TrialRarity, scale=TRUE, center=TRUE))
  )


# Linear mixed-effects model mixed-effects model: Trial rarity
TrialRarity_model <- glmmTMB(
  TrialRarity_scaled ~ male_type,
  family = gaussian(),
  data = rarity
)

# Overall test of male type
car::Anova(TrialRarity_model, type = "II")

# Estimated marginal means
TrialRarity_emm <- emmeans::emmeans(
  TrialRarity_model,
  ~ male_type
)

# Pairwise comparisons with Bonferoni adjustment
pairs(TrialRarity_emm, adjust = "holm", type="response")

# Extract estimates for plotting
TrialRarity_emm <- as.data.frame(
  summary(TrialRarity_emm)
)

# Model diagnostics
res <- DHARMa::simulateResiduals(TrialRarity_model)
plot(res)

# Plotting
FigS4.1 <- ggplot(rarity,
       aes(x = male_type, y = TrialRarity_scaled)) +
  
  # Individual male observations
  geom_jitter(
    width = 0.12,
    height = 0,
    color = "grey60",
    alpha = 0.5,
    size = 2
  ) +
  
  # Model estimated means
  geom_point(
    data = TrialRarity_emm,
    aes(x = male_type, y = emmean, color = male_type),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # Standard errors
  geom_errorbar(
    data = TrialRarity_emm,
    aes(
      x = male_type,
      ymin = emmean - SE,
      ymax = emmean + SE,
      color = male_type
    ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "common" = "orchid",
      "rare" = "orange",
      "migrant" = "#009E73"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "common" = "Resident Common",
      "rare" = "Resident Rare",
      "migrant" = "Migrant Unfamiliar"
    )
  ) +
  
  labs(
    x = "Male Type",
    y = "Male Color Pattern Rarity in Trial (scaled)"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15),
    legend.position = "none"
  )

ggsave(filename = here("figures", "FigS4.1.png"),
       plot = FigS4.1,
       width = 8,
       height = 6,
       dpi = 300)

# Linear mixed-effects model mixed-effects model: Pool rarity
PoolRarity_model <- glmmTMB(
  PoolRarity_scaled ~ male_type,
  family = gaussian(),
  data = rarity
)

summary(PoolRarity_model)

# Overall test of male type
car::Anova(PoolRarity_model, type = "II")

# Estimated marginal means
PoolRarity_emm <- emmeans::emmeans(
  PoolRarity_model,
  ~ male_type
)

# Pairwise comparisons with Bonferoni adjustment
pairs(PoolRarity_emm, adjust = "holm", type="response")

# Extract estimates for plotting
PoolRarity_emm <- as.data.frame(
  summary(PoolRarity_emm)
)

# Model diagnostics
res <- DHARMa::simulateResiduals(PoolRarity_model)
plot(res)

# Plotting S4: pool rarity vs male type
FigS4.2 <- ggplot(rarity,
       aes(x = male_type, y = PoolRarity_scaled)) +
  
  # Individual male observations
  geom_jitter(
    width = 0.12,
    height = 0,
    color = "grey60",
    alpha = 0.5,
    size = 2
  ) +
  
  # Model estimated means
  geom_point(
    data = PoolRarity_emm,
    aes(x = male_type, y = emmean, color = male_type),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # Standard errors
  geom_errorbar(
    data = PoolRarity_emm,
    aes(
      x = male_type,
      ymin = emmean - SE,
      ymax = emmean + SE,
      color = male_type
    ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "common" = "orchid",
      "rare" = "orange",
      "migrant" = "#009E73"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "common" = "Resident Common",
      "rare" = "Resident Rare",
      "migrant" = "Migrant Unfamiliar"
    )
  ) +
  
  labs(
    x = "Male Type",
    y = "Male Color Pattern Rarity in Pool (scaled)"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15),
    legend.position = "none"
  )


ggsave(filename = here("figures", "FigS4.2.png"),
       plot = FigS4.2,
       width = 8,
       height = 6,
       dpi = 300)

# Calculate mean rarity for each male type within trials
trial_comparison <- rarity %>%
  dplyr::filter(male_type %in% c("common", "migrant")) %>%
  dplyr::group_by(Trial, male_type) %>%
  dplyr::summarise(
    MeanPoolRarity = mean(PoolRarity),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = male_type,
    values_from = MeanPoolRarity
  ) %>%
  dplyr::mutate(
    MigrantLessRare = migrant < common,
    MigrantMoreRare = migrant > common
  )

# Number of trials with migrants less rare than common residents
sum(trial_comparison$MigrantLessRare)

# Number of trials with migrants more rare than common residents
sum(trial_comparison$MigrantMoreRare)


# Compare migrant rarity against both common males in each trial
trial_comparison <- combined_rarity %>%
  dplyr::filter(type %in% c("RC1", "RC2", "MU")) %>%
  dplyr::select(Trial, type, TrialRarity) %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = TrialRarity
  ) %>%
  dplyr::mutate(
    n_common_more_rare = (MU < RC1) + (MU < RC2),
    category = dplyr::case_when(
      n_common_more_rare == 0 ~ "Neither common male",
      n_common_more_rare == 1 ~ "One common male",
      n_common_more_rare == 2 ~ "Both common males"
    )
  )

# Number of trials in each category
trial_comparison %>%
  dplyr::count(category)

# A nice way to look at individual trial pairwise distances
get_trial_distances <- function(trial_number) {
  
  trial_males <- pool_info %>%
    dplyr::filter(trial == trial_number) %>%
    dplyr::select(ID, type)
  
  trial_colors <- trial_males %>%
    dplyr::left_join(
      colors %>%
        dplyr::select(ID, dplyr::all_of(color_vars)),
      by = "ID"
    ) %>%
    dplyr::arrange(match(type, c("RC1", "RC2", "RR", "MU")))
  
  stopifnot(
    nrow(trial_colors) == 4L,
    !anyNA(trial_colors[, color_vars, drop = FALSE])
  )
  
  color_matrix <- as.matrix(
    trial_colors[, color_vars, drop = FALSE]
  )
  
  rownames(color_matrix) <- trial_colors$type
  
  round(
    as.matrix(dist(color_matrix, method = "euclidean")),
    2
  )
}

get_trial_distances(11)
