# Load required packages
library(psych)
library(dplyr)
library(tidyr)

# Necessary input files:
# rgb_uv_pool5DT1PR.rds

#WARNING: This clears your global working environment to start fresh between scripts
rm(list = ls())

#Read in color data 
rgb_uv_pool5DT1PR <- readRDS(here::here("data", "rgb_uv_pool5DT1PR.rds"))

# Define color columns
color_cols <- 9:28832

# Extract full color matrix
dat5DT1PR_raw <- rgb_uv_pool5DT1PR[, color_cols]

# Globally standardize each color variable across all individuals
# This gives each color variable mean = 0 and SD = 1 across the full dataset
dat5DT1PR_scaled <- scale(dat5DT1PR_raw, center = TRUE, scale = TRUE)


dat5DT1PR_scaled <- data.frame(
  Stream = rgb_uv_pool5DT1PR$Stream,
  dat5DT1PR_scaled
)


# if you get a working memory error between stream calculations, remove previous stream matrix to make more storage space
#Aripo Intro
Aripo_df <- subset(dat5DT1PR_scaled, Stream=="Aripo")
# Color variables only
Aripo.colors <- Aripo_df[, 2:28825]
Aripo_var <- cov(Aripo.colors)
tr(Aripo_var)

#Nar
Nar_df <- subset(dat5DT1PR_scaled, Stream=="Nar")
# Color variables only
Nar.colors <- Nar_df[, 2:28825]
Nar_var <- cov(Nar.colors)
tr(Nar_var)

#Main
Main_df <- subset(dat5DT1PR_scaled, Stream=="Temp")
# Color variables only
Main.colors <- Main_df[, 2:28825]
Main_var <- cov(Main.colors)
tr(Main_var)
