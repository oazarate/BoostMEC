library(tidyverse)
library(lightgbm)
library(markovchain)
library(TmCalculator)
library(janitor)

set.seed(50)

#Required functions
#Sequence manipulation, Markov features, data preparation
function_scripts <- list.files("functions", full.names = TRUE)
sapply(function_scripts, source)

#Model objects
ratio_mat_2nd <- readRDS("model/ratio-matrices-order-2.rds")
lgbm_mod <- lgb.load("model/model.txt")
lgbm_dat_rules <- readRDS("model/training_rules.rds")
features <- read_csv("model/selected_features.csv") %>% pull(features)

#Sequences
dat <- read_csv("data/sample_dataset.csv")
#Computed free energy
dat_fe <- read_csv("data/min_free_energy.csv")

#Format data
data_processed <- prepare_data(df = dat,
                               free_energy_df = dat_fe,
                               ratio_list = ratio_mat_2nd,
                               selected_features = features,
                               encoding_rules = lgbm_dat_rules)

boostmec_predictions <- predict(lgbm_mod, data_processed)

write_csv(tibble(predicted_efficiency = boostmec_predictions), "predictions/sample_dataset_predictions.csv")

#If observed efficiency is provided in sample_dataset.csv, run the lines below to see agreement via Spearman correlation:

#cor(dat$efficiency, boostmec_predictions, method = "spearman") %>%
#  round(3)
