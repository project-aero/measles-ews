# ews-re-corrs-across-bandwidths.R
#  Script to calculate the Spearman rank correlation between candidate EWS
#  and effective reproductive ratio through time at different window sizes
#  for the EWS.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(spaero)
library(DescTools)


# Load data and model results ---------------------------------------------

file_name <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
measles_data <- readRDS(file_name) %>%
  filter(year > 1994) %>%
  dplyr::select(-year, -week_of_year, -obs_week, -time) %>%
  mutate(
    region = str_sub(region, 1, str_length(region)-7)  # remove (City) from name
  )

focal_states <- c(
  "effective_r_nonseasonal",
  "effective_r_seasonal",
  "transmission_rate"
)

filtered_states <- {}
for(do_city in unique(measles_data$region)){
  tmp_states <- readRDS(paste0("../results/filtered-states-", do_city, ".RDS")) %>%
    filter(state %in% focal_states) %>%
    unnest() %>%
    dplyr::select(date, med, state) %>%
    rename(state_value = med) %>%
    mutate(
      region = do_city
    )
  
  filtered_states <- bind_rows(filtered_states, tmp_states)
}
