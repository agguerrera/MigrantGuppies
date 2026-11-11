# Load required packages
library(tidyverse)
library(patchwork)
library(zoo)
library(glmmTMB)
library(car)
library(emmeans)
library(psych)
library(broom.mixed)
library(sjPlot)
library(performance)
library(DHARMa)
library(MASS)
library(buildmer)
library(here)


# Necessary input files:
# combined_rarity.csv : color and trial-specific info
# focal_behaviors.csv : behavior information
# annotated-behaviors.csv : time-stamped behaviors for supplemental

#WARNING: This clears your global working environment to start fresh between scripts
rm(list = ls())

# ============================================================
# BEHAVIORAL STATISTICS
# ============================================================
behavior_df <- read.csv(here("data", "annotated-behaviors.csv"))
focal_behaviors <- read.csv(here("data", "focal_behaviors.csv"))
combined_rarity <- read.csv(here("data", "combined_rarity.csv"))


# Merge behavioral data with pool and trial rarity
trin_df <- combined_rarity %>%
  dplyr::left_join(
    focal_behaviors,
    by = c("ID", "Trial"),
    relationship = "many-to-one"
  )

#Calling variables
trin_df$trial <- as.factor(trin_df$Trial)
trin_df$unique_id <- as.factor(trin_df$unique_id)
trin_df$mig_res <- as.factor(trin_df$mig_res)
trin_df$time_of_day <- as.numeric(trin_df$time_of_day)
trin_df$female_id <- as.factor(trin_df$female_id)
trin_df$male_id <- as.factor(trin_df$male_id)
trin_df$focal <-as.factor(trin_df$focal.x)
trin_df$male_type <- as.factor(trin_df$male_type)
trin_df$pool_rarity_f <- as.numeric(trin_df$PoolRarity) # Rarity from female perspective
trin_df$pool_rarity_f_scaled <- scale(trin_df$pool_rarity_f, center=TRUE, scale=TRUE) # Rarity from female perspective
trin_df$trial_rarity <- as.numeric(trin_df$TrialRarity) # Rarity from female perspective
trin_df$trial_rarity_scaled <- scale(trin_df$trial_rarity, scale=TRUE, center=TRUE) # Rarity from female perspective
trin_df$TotalOr <- as.numeric(trin_df$TotalOr)
trin_df$uv_sum <- as.numeric(trin_df$uv_sum)
trin_df$TotalOr_scaled <- scale(trin_df$TotalOr,center=TRUE, scale=TRUE)
trin_df$TotalBlk_scaled <- scale(trin_df$TotalBlk,center=TRUE, scale=TRUE)
trin_df$uv_sum_scaled <- scale(trin_df$uv_sum, center=TRUE, scale=TRUE)
trin_df$sigmoids_comp <- as.numeric(trin_df$sigmoids_comp)
trin_df$total_sigmoids <- as.numeric(trin_df$sigmoids_alone + trin_df$sigmoids_comp)
trin_df$SL <- as.numeric(trin_df$sl)
trin_df$female_pool <- as.factor(trin_df$female_pool)
trin_df$stream <- as.factor(trin_df$stream)
trin_df$SL_scaled <- scale(trin_df$SL, center=TRUE, scale=TRUE)
trin_df$total_sigmoids_scaled <- as.numeric(scale(trin_df$total_sigmoids, center=TRUE, scale=TRUE))
trin_df$Interest <- trin_df$Orient_to_Active_Follow + trin_df$Female_Approach_While_Following  + trin_df$Orient_Response_to_Sigmoid + trin_df$Approach_Response_to_Sigmoid + trin_df$Orient_to_No_Follow + trin_df$Female_Approach_No_Follow #positive response to a male 

#################################
####Female interest analysis####
################################
# Female interest model
interest.p <- glmmTMB(Interest 
                       ~ male_type
                       + total_sigmoids_scaled
                       + stream 
                       + stream*total_sigmoids_scaled 
                       + (1|unique_id)
                       + (1|trial), data = trin_df, family = poisson())
summary(interest.p)
Anova(interest.p,type=c("III")) 
VarCorr(interest.p) #Pull out random effects values
tidy(interest.p, effects = "fixed", conf.int = TRUE) #Pull out CIs

#Test model assumptions
sim_interest <- simulateResiduals(fittedModel = interest.p)
testZeroInflation(sim_interest) 
testDispersion(sim_interest) 
testCategorical(sim_interest, trin_df$male_type, plot=T)
testCategorical(sim_interest, trin_df$stream, plot=T)
testResiduals(sim_interest)
outliers(sim_interest)
testOutliers(sim_interest)
check_collinearity(interest.p)

#Estimate marginal means of categorical predictors
#Male type
emm_interest <- emmeans(interest.p, ~ male_type, type = "response")
pairs(emm_interest, adjust="holm", type="response")
confint(emm_interest)

#Stream
courtship_trends <- emtrends(interest.p,~ stream, var = "total_sigmoids_scaled")
confint(courtship_trends)
summary(courtship_trends, infer = TRUE)
pairs(courtship_trends, adjust = "holm")

#Figure 3
#Plotting predicted results for male type
emm_interest <- data.frame(emm_interest)
emm_interest$male_type <- factor(emm_interest$male_type, levels = c("common", "rare", "migrant"))
interest_df <- trin_df
interest_df$male_type <- factor(interest_df$male_type,levels = c("common", "rare", "migrant"))

Fig3a <- ggplot() +
  #geom_point(data = interest_df, aes(x = male_type, y = Interest), color="grey", position = position_jitter(width = 0.1, height = 0), alpha = 0.6, size = 2) + # include this line for Supplemental Figure
  geom_point(data = emm_interest, 
             aes(x = male_type, y = rate, color = male_type), 
             position = position_dodge(width = 0.4), size = 2, inherit.aes = FALSE) +
  geom_errorbar(data = emm_interest, aes(x = male_type, ymin = rate - SE, ymax = rate + SE, color = male_type),
               width = 0.3, position = position_dodge(width = 0.4), inherit.aes = FALSE) +
  labs(x = "Male Type", y = "Female interest") +
  scale_x_discrete(limits = c("common", "rare", "migrant"),
    labels = c(
      "common" = "Resident\nCommon",
      "rare" = "Resident\nRare",
      "migrant" = "Migrant\nUnfamiliar")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1.5)) + #change to 8 for full raw distribution
  theme(axis.text = element_text(size = 15, color = "black"),axis.title = element_text(size = 15)) + 
  scale_color_manual(values = c("orchid", "orange", "#009E73"))

Fig3a

ggsave(filename = here("figures", "Fig3a.png"),
  plot = Fig3a,
  width = 8,
  height = 6,
  dpi = 300)

# Plotting predicted results for stream x courtship effort
# Define representative range of courtship effort
sigmoid_range <- quantile(
  trin_df$total_sigmoids_scaled,
  probs = c(0.05, 0.95),
  na.rm = TRUE)

# Generate predictions
sigmoid_seq <- seq(
  sigmoid_range[1],
  sigmoid_range[2],
  length.out = 100)

# Estimated marginal means across courtship effort
emm_interaction <- emmeans(interest.p,~ stream | total_sigmoids_scaled,at = list(total_sigmoids_scaled = sigmoid_seq),type = "response")
emm_interaction_df <- as.data.frame(summary(emm_interaction, type = "response"))

# Set stream order
emm_interaction_df$stream <- factor(emm_interaction_df$stream,levels = c("ARIPO", "NAR", "TEMP"))

# Plot interaction
Fig3b <- ggplot(emm_interaction_df,aes(x = total_sigmoids_scaled,y = rate,linetype = stream)) +
  geom_ribbon(aes(ymin = asymp.LCL,ymax = asymp.UCL,group = stream),fill = "grey70",alpha = 0.3,color = NA) + 
  geom_line(color = "black",linewidth = 1.2) +
  scale_linetype_manual(values = c("ARIPO" = "solid","NAR" = "dashed","TEMP" = "dotted"),
    labels = c("ARIPO" = "Aripo Intro","NAR" = "Naranjo","TEMP" = "Mainstream")) +
  labs(x = "Courtship Effort (Scaled Total Sigmoids)",y = "Predicted Female Interest",linetype = "Stream") +
  theme_classic() +
  theme(axis.text = element_text(size = 15, color = "black"),axis.title = element_text(size = 15),legend.text = element_text(size = 13),legend.title = element_text(size = 13),legend.position = "right")

ggsave(filename = here("figures", "Fig3b.png"),
       plot = Fig3b,
       width = 8,
       height = 6,
       dpi = 300)

########################
# Complimentary analysis of female interest in migrants vs residents with continuous rarity metrics
########################
# Continuous trial-level rarity
trarity <- glmmTMB(Interest 
                    ~ mig_res
                    + trial_rarity_scaled
                    + stream
                    + total_sigmoids_scaled
                    + stream*total_sigmoids_scaled 
                    + (1|unique_id) 
                    + (1|trial), data = trin_df, family = poisson())
summary(trarity)
Anova(trarity, type="III")
VarCorr(trarity) #Pull out random effects values
tidy(trarity, effects = "fixed", conf.int = TRUE) #Pull out CIs

#Check model assumptions: 
rarity_comp <- simulateResiduals(fittedModel = trarity)
testZeroInflation(rarity_comp) 
testDispersion(rarity_comp) 
testCategorical(rarity_comp, trin_df$mig_res, plot=T)
testCategorical(rarity_comp, trin_df$stream, plot=T)
testResiduals(rarity_comp)
outliers(rarity_comp)
testOutliers(rarity_comp)
check_collinearity(trarity)

# Migrant status effect
#Estimate marginal means 
emm_trarity <- emmeans(trarity, ~ mig_res,type = "response")
pairs(emm_trarity, adjust = "holm", type = "response")

# Stream x courtship effort effect
# Estimated slopes 
tr.sig.trends <- emtrends(trarity, ~ stream, var = "total_sigmoids_scaled")
summary(tr.sig.trends)
confint(tr.sig.trends)
# Pairwise slope contrasts with Holm correction
summary(pairs(tr.sig.trends, adjust = "holm"))


# Continuous pool-level rarity
# Final pool rarity model
prarity <- glmmTMB(Interest
                     ~ mig_res
                     + pool_rarity_f_scaled
                     + total_sigmoids_scaled
                     + stream
                     + stream*total_sigmoids_scaled
                     + (1|unique_id)
                     + (1|trial), data = trin_df, family = poisson())
summary(prarity)
Anova(prarity, type="III")
VarCorr(prarity) #Pull out random effects values
tidy(prarity, effects = "fixed", conf.int = TRUE) #Pull out CIs

#Check model assumptions: 
prarity_comp <- simulateResiduals(fittedModel = prarity)
testZeroInflation(prarity_comp) 
testDispersion(prarity_comp) 
testCategorical(prarity_comp, trin_df$mig_res, plot=T)
testCategorical(prarity_comp, trin_df$stream, plot=T)
testResiduals(prarity_comp)
outliers(prarity_comp)
testOutliers(prarity_comp)
check_collinearity(prarity)

# Migrant status effect
#Estimate marginal means 
emm_prarity <- emmeans(prarity, ~ mig_res,type = "response")
pairs(emm_prarity, adjust = "holm", type = "response")

# Stream x courtship effort effect
# Estimated slopes 
p.sig.trends <- emtrends(prarity, ~ stream, var = "total_sigmoids_scaled")
summary(p.sig.trends)
confint(p.sig.trends)
# Pairwise slope contrasts with Holm correction
summary(pairs(p.sig.trends, adjust = "holm"))


#Figure 4: female interest in migrant type (combining both models)
rarity_df <- data.frame(trin_df)
rarity_df$mig_res <- factor(rarity_df$mig_res, levels = c("Res", "Mig"))
emm_prarity <- data.frame(emm_prarity)
emm_prarity$mig_res <- factor(emm_prarity$mig_res,levels = c("Res", "Mig"))

# Combine estimated marginal means
emm_combined <- bind_rows(
  as.data.frame(emm_prarity) %>%
    mutate(Rarity = "Pool rarity"),
  as.data.frame(emm_trarity) %>%
    mutate(Rarity = "Trial rarity"))

# Set factor levels
emm_combined$mig_res <- factor(
  emm_combined$mig_res,
  levels = c("Res", "Mig"))

emm_combined$Rarity <- factor(
  emm_combined$Rarity,
  levels = c("Pool rarity", "Trial rarity"))

# Position dodge
pd <- position_dodge(width = 0.25)

# Plot
Fig4 <- ggplot(
  emm_combined,
  aes(x = mig_res, y = rate, color = mig_res,shape = Rarity, group = Rarity)) +
  geom_point(size = 3,position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = rate - SE, ymax = rate + SE), width = 0.12, position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c("Res" = "darkgrey", "Mig" = "#009E73")) +
  scale_shape_manual(values = c(
      "Pool rarity" = 16,  # Circle
      "Trial rarity" = 15)) +  # Square
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(
    x = "Migrant Status",
    y = "Predicted Female Interest",
    color = "Migrant Status",
    shape = "Rarity metric") +
  theme_classic() +
  theme(axis.text = element_text(size = 15, color = "black"),axis.title = element_text(size = 15),legend.text = element_text(size = 12),legend.title = element_text(size = 12))

ggsave(filename = here("figures", "Fig4.png"),
       plot = Fig4,
       width = 8,
       height = 6,
       dpi = 300)


########################################
#####Male-male competition analysis ####
########################################
compfollows <-glmmTMB(comp_follows
                     ~ male_type 
                     + stream
                     + TotalBlk_scaled
                     + uv_sum_scaled   
                     + male_type*TotalBlk_scaled 
                     + stream*uv_sum_scaled
                     + (1|female_pool)
                     + (1|unique_id)  
                     + (1|trial), zi=~1, family=ziGamma(link="log"), data=trin_df)
summary(compfollows)
Anova(compfollows,type=c("III")) 
VarCorr(compfollows)
tidy(compfollows, effects = "fixed", conf.int = TRUE) #Pull out CIs

#Check model assumptions: 
sim_comp <- simulateResiduals(fittedModel = compfollows)
testZeroInflation(sim_comp) 
testDispersion(sim_comp) 
testCategorical(sim_comp, trin_df$male_type, plot=T)
testCategorical(sim_comp, trin_df$stream, plot=T)
testResiduals(sim_comp)
outliers(sim_comp)
testOutliers(sim_comp)
check_collinearity(compfollows)

# Estimated slopes for each male type
comp_trends <- emtrends(compfollows,~ male_type,var = "TotalBlk_scaled")
summary(comp_trends)
confint(comp_trends)

# Pairwise contrasts with Holm correction
summary(pairs(comp_trends, adjust = "holm"))


#Figure 4
#Plotting predicted results interaction between male type x black
blk_seq <- seq(min(trin_df$TotalBlk_scaled, na.rm = TRUE),
               max(trin_df$TotalBlk_scaled, na.rm = TRUE),
               length.out = 100)

emm_cont <- emmeans(
  compfollows,
  ~ male_type | TotalBlk_scaled,
  at = list(TotalBlk_scaled = blk_seq),
  type = "response")

pairs(emm_cont)
emm_cont_df <- as.data.frame(emm_cont)

Fig5 <- ggplot(emm_cont_df,aes(x = TotalBlk_scaled, y = response, color = male_type, fill = male_type)) +
  #geom_point(data=trin_df, aes(x=TotalBlk_scaled, y=comp_follows, color=male_type, alpha=.4)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL),alpha = 0.2, color = NA) +
  geom_line(aes(y = asymp.LCL),linetype = "dotted", linewidth = 0.5) +
  geom_line(aes(y = asymp.UCL),linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 1) + #predicted relationship
  theme_classic() +
  labs(x = "Total Black (scaled)",
       y = "Predicted Competitive Following",
       color = "Male type",
       fill = "Male type") +
  scale_color_manual(values = c(
    "common" = "orchid",
    "rare" = "orange",
    "migrant" = "#009E73")) +
  scale_fill_manual(values = c(
    "common" = "orchid",
    "rare" = "orange",
    "migrant" = "#009E73"))

ggsave(filename = here("figures", "Fig5.png"),
       plot = Fig5,
       width = 8,
       height = 6,
       dpi = 300)


####Supplemental Figure Plotting####
# Supplemental Figure S3: Lagged female interest
behavior_df <- behavior_df %>%
  mutate(
    trial = as.character(trial),
    focal = as.character(focal)
  )

# Identify male-following bouts
annotated_df <- behavior_df %>%
  group_by(trial, focal) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    # Time at which male following ends
    follow_stop_time = if_else(
      behavior == "focal_follows_female_alone" &
        status == "STOP",
      as.numeric(time),
      NA_real_
    )
  ) %>%
  # Carry most recent STOP timestamp forward
  tidyr::fill(follow_stop_time, .direction = "down") %>%
  mutate(
    # Is the male currently following?
    male_following =
      cumsum(
        behavior == "focal_follows_female_alone" &
          status == "START"
      ) >
      cumsum(
        behavior == "focal_follows_female_alone" &
          status == "STOP"
      ),
    
    # Seconds since following ended
    time_since_follow = if_else(
      !male_following,
      as.numeric(time) - follow_stop_time,
      NA_real_
    )
  ) %>%
  ungroup()

#Identify lagged female interest
lagged_interest <- annotated_df %>%
  filter(
    behavior %in% c(
      "female_orients_towards_focal",
      "female_swims_towards_male"
    ),
    !male_following,
    !is.na(time_since_follow),
    time_since_follow >= 0
  )

# Plot time since male pursuit and female orientation
FigS3<- ggplot(lagged_interest,
       aes(x = time_since_follow)) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    closed = "left",
    color = "grey",
    fill = "grey50"
  ) +
  scale_x_continuous(
    breaks = scales::breaks_width(10)
  ) +
  labs(
    x = "Time since last male follow (s)",
    y = "Count of female interest events"
  ) +
  theme_classic()

ggsave(filename = here("figures", "FigS3.png"),
       plot = FigS3,
       width = 8,
       height = 6,
       dpi = 300)

#Supplemental Figure S6: Trial-level rarity of migrant males compared to resident common males
rarity_diff <- trin_df %>%
  group_by(trial) %>%
  mutate(migrant_rarity = trial_rarity_scaled[male_type == "migrant"]) %>%
  filter(male_type == "common") %>%
  mutate(rarity_difference = migrant_rarity - trial_rarity_scaled, type = c("RC1", "RC2")) %>%
  ungroup()

rarity_diff_plot <- rarity_diff %>%
  group_by(trial) %>%
  mutate(rarity_category = case_when(all(rarity_difference > 0) ~ "Migrant more rare than both", all(rarity_difference < 0) ~ "Migrant less rare than both", TRUE ~ "Migrant intermediate")) %>%
  ungroup()

# Set category order
rarity_diff_plot <- rarity_diff_plot %>%
  mutate(rarity_category = factor(rarity_category, levels = c(
        "Migrant more rare than both",
        "Migrant intermediate",
        "Migrant less rare than both")))

# Get trial order: category first, then numerical trial number
trial_order <- rarity_diff_plot %>%
  distinct(trial, rarity_category) %>%
  arrange(rarity_category, trial) %>%
  pull(trial)

# Apply that order to the y-axis
rarity_diff_plot <- rarity_diff_plot %>%
  mutate(trial_plot = factor(trial, levels = trial_order))

FigS6 <- ggplot(rarity_diff_plot,
       aes(x = rarity_difference,
           y = trial_plot,
           group = trial)) +
  geom_vline(xintercept = 0,linetype = "dashed",linewidth = 0.6) +
  geom_line(aes(color = rarity_category),alpha = 0.4) +
  geom_point(aes(shape = type,color = rarity_category),size = 2.5) +
  scale_color_manual(values = c(
      "Migrant more rare than both" = "deeppink",
      "Migrant intermediate" = "grey50",
      "Migrant less rare than both" = "forestgreen")) +
  labs(x = "Migrant rarity − common male rarity",y = "Trial",shape = "Common male",color = NULL ) + theme_classic()

ggsave(filename = here("figures", "FigS6.png"),
       plot = FigS6,
       width = 8,
       height = 6,
       dpi = 300)

# Supplemental Figure S7: Color distributions across streams
# UV plot
p_uv <- ggplot(trin_df, aes(stream, uv_sum)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white",median.color = "purple") +
  geom_point(position = position_jitter(width = 0.1, height = 0),color = "purple", alpha = 0.6) +
  theme_classic() + labs(x = "Tributary", y = "Total UV Reflectance")

# Black plot
p_blk <- ggplot(trin_df, aes(stream, TotalBlk)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white",median.color = "black") +
  geom_point(position = position_jitter(width = 0.1, height = 0),color = "black", alpha = 0.6) +
  theme_classic() + labs(x = "Tributary", y = "Proportion of Black")

# Orange plot
p_or <- ggplot(trin_df, aes(stream, TotalOr)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", median.color = "darkorange") +
  geom_point(position = position_jitter(width = 0.1, height = 0), color = "darkorange", alpha = 0.6) +
  theme_classic() + labs(x = "Tributary", y = "Proportion of Orange")

# Put orange and black (proportions) on same axis
y_max <- max(c(trin_df$TotalBlk, trin_df$TotalOr), na.rm = TRUE)
y_min <- min(c(trin_df$TotalBlk, trin_df$TotalOr), na.rm = TRUE)
y_breaks <- pretty(c(y_min, y_max), n = 5)

p_blk <- p_blk + scale_y_continuous(breaks = y_breaks) + coord_cartesian(ylim = c(y_min, y_max))
p_or <- p_or + scale_y_continuous(breaks = y_breaks) + coord_cartesian(ylim = c(y_min, y_max))

# Combine plots  
FigS7 <- p_uv | p_blk | p_or

ggsave(filename = here("figures", "FigS7.png"),
       plot = FigS7,
       width = 8,
       height = 6,
       dpi = 300)

# Supplemental Figure S9: Sigmoid and following time across male type
trin_df$male_type <- factor(trin_df$male_type, levels = c("common", "rare", "migrant"))

FigS9.1 <- ggplot(trin_df, aes(x = male_type, y = follows_fem)) +
  stat_summary(fun = mean, geom = "col", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.6) +
  theme_classic() +
  labs(x = "Male type", y = "Male follows female (s) (Mean ± SE)") +
  theme(legend.position = "none")

ggsave(filename = here("figures", "FigS9.1.png"),
       plot = FigS9.1,
       width = 8,
       height = 6,
       dpi = 300)

FigS9.2 <- ggplot(trin_df, aes(x = male_type, y = total_sigmoids)) +
  stat_summary(fun = mean, geom = "col", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.6) +
  theme_classic() +
  labs(x = "Male type", y = "Male sigmoid displays (Mean ± SE)") +
  theme(legend.position = "none")

ggsave(filename = here("figures", "FigS9.2.png"),
       plot = FigS9.2,
       width = 8,
       height = 6,
       dpi = 300)

#Supplemental Figure S10: Predictor correlelogram
predictor_data <- trin_df[, c("total_sigmoids_scaled", "TotalOr_scaled",
                              "TotalBlk_scaled", "uv_sum_scaled",
                              "SL_scaled", "pool_rarity_f_scaled",
                              "trial_rarity_scaled")]

colnames(predictor_data) <- c(
  "Sigmoids (scaled)",
  "Proportion of orange (scaled)",
  "Proportion of black (scaled)",
  "UV reflectance (scaled)",
  "Male size (scaled)",
  "Pool rarity",
  "Trial rarity"
)

# Export figure
png(filename = here("figures", "FigS10.png"),
    width = 8, height = 6, units = "in", res = 300)

psych::pairs.panels(
  predictor_data,
  method = "pearson",
  hist.col = "grey",
  density = TRUE,
  ellipses = FALSE)

dev.off()

#Supplemental Analysis & Figure S11: male type x black on number of interruptive bouts 
int <-glmmTMB(interrupts_no
              ~ male_type 
              + stream
              + TotalBlk_scaled
              + uv_sum_scaled   
              + male_type*TotalBlk_scaled 
              + stream*uv_sum_scaled
              + (1|unique_id) 
              + (1|trial), zi=~1, family=nbinom2(), data=trin_df)
summary(int)
Anova(int,type=c("III")) 

#Estimate marginal means
emtrends(int, pairwise ~ male_type, var = "TotalBlk_scaled", infer=T, adjust="holm")

FigS11 <- ggplot(trin_df,
       aes(x = TotalBlk_scaled, y = interrupts_no,
           color = male_type, fill = male_type)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(
    method = "glm",
    method.args = list(family = negative.binomial(theta = 1)),
    se = TRUE,
    linewidth = 0.9,
    alpha = 0.18
  ) +
  scale_color_manual(
    values = c(
      "common" = "orchid",
      "rare" = "orange",
      "migrant" = "#009E73"
    ),
    labels = c(
      "common" = "Resident Common",
      "rare" = "Resident Rare",
      "migrant" = "Migrant Unfamiliar"
    )
  ) +
  scale_fill_manual(
    values = c(
      "common" = "orchid",
      "rare" = "orange",
      "migrant" = "#009E73"
    ),
    guide = "none"
  ) +
  labs(
    x = "Proportion of black",
    y = "No. of interruptions",
    color = "Male type"
  ) +
  theme_classic()

ggsave(filename = here("figures", "FigS11.png"),
       plot = FigS11,
       width = 8,
       height = 6,
       dpi = 300)



