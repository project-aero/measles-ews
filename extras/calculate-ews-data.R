# calculate-ews-data.R
#  Script to calculate early warning signals over a moving window based
#  on weekly case incidence data from four Nigerien cities. EWS are
#  calculated using the spaero::get_stats() function. We use a bandwidth
#  of 35 weeks, meaning each window is 35*2 = 70 weeks in length.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(spaero)


# Load data ---------------------------------------------------------------

fname <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
measles_data <- readRDS(fname) %>%
  filter(year > 1994)  # drop first NA year, only used for modeling


# Calculate EWS -----------------------------------------------------------

all_stats <- tibble()  # empty tibble for storage

for(do_region in unique(measles_data$region)){
  
  cases <- measles_data %>%
    filter(region == do_region) %>%
    pull(cases)
  
  city_stats <- spaero::get_stats(
    x = cases,
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
    mutate(
      time_iter = 1:n(),
      date = unique(measles_data$date)
    ) %>%
    gather(key = ews, value = value, -time_iter, -date) %>%
    mutate(region = do_region)
  
  all_stats <- bind_rows(all_stats, city_stats_tb)
}


# Plot data and EWS -------------------------------------------------------

plot_it <- function(dates, cases, ews, lab, title){
  par(mar = c(5,5,2,5))
  plot(dates, cases, type = "l", main = title)
  par(new = T)
  plot(dates, ews, type = "l", col = "red", axes=F, xlab=NA, ylab=NA)
  axis(side = 4)
  mtext(side = 4, line = 3, lab)
}

dates <- unique(measles_data$date) 
for(do_region in unique(measles_data$region)){
  par(mfrow = c(2,5))
  for(do_ews in unique(all_stats$ews)){
    cases <- filter(measles_data, region == do_region) %>% pull(cases)
    ews <- filter(all_stats, region == do_region & ews == do_ews) %>% pull(value)
    plot_it(dates, cases, ews, lab = do_ews, title = do_region)
  }
}




