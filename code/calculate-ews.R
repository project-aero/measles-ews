# calculate-ews.R:
#  Script to calculate 10 early warning signals for time series incidence
#  data from four cities in Niger.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(spaero)


# Load data ---------------------------------------------------------------

file_name <- "../data/clean-data/weekly-measles-and-demog-niger-cities-clean.RDS"
measles_data <- readRDS(file_name) %>%
  dplyr::select(-population_smooth, -births_per_week_smooth)


# Calculate early warning signals -----------------------------------------

all_stats <- {}
for(do_city in unique(measles_data$region)){
  city_data <- measles_data %>%
    filter(region == do_city)
  
  city_stats <- spaero::get_stats(
    x = city_data$cases,
    center_trend = "local_constant", 
    center_kernel = "uniform", 
    center_bandwidth = 35, 
    stat_trend = "local_constant", 
    stat_kernel = "uniform", 
    stat_bandwidth = 35, 
    lag = 1, 
    backward_only = TRUE
  )$stats
  
  city_stats_tb <- as_tibble(city_stats) %>%
    mutate(time = 1:n()) %>%
    gather(key = ews, value = value, -time) %>%
    mutate(city = do_city)
  
  all_stats <- bind_rows(all_stats, city_stats_tb)
}


