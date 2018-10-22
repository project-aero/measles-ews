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
library(ggthemes)


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
    # filter(state %in% focal_states) %>%
    unnest() %>%
    dplyr::select(date, med, state) %>%
    rename(state_value = med) %>%
    mutate(
      region = do_city
    )
  
  filtered_states <- bind_rows(filtered_states, tmp_states)
}


z <- filtered_states %>%
  filter(state == "effective_r_seasonal" & region == "Niamey") %>%
  pull(state_value)

mvz <- spaero::get_stats(
  x = log(z),
  center_trend = "local_constant", 
  center_kernel = "uniform", 
  center_bandwidth = 26, 
  stat_trend = "local_constant", 
  stat_kernel = "uniform", 
  stat_bandwidth = 26, 
  lag = 1, 
  backward_only = TRUE
)$stats$mean

# plot(z, type = "l", col = "blue")
# lines(exp(mvz), col = "red")

y <- filtered_states %>%
  filter(state == "cases" & region == "Niamey") %>%
  pull(state_value)
# 
# plot(y, type = "p")
# lines(exp(mvz)*500, col = "red")

# Calculate EWS and correlations at different bandwidths ------------------

bandwidth_vector <- c(13, 26, 52, 78, 104)
correlations <- tibble()
for(do_bw in bandwidth_vector){
  
  window_bandwidth <- do_bw
  all_stats <- {}
  
  for(do_city in unique(measles_data$region)){
    city_reff <- filtered_states %>%
      filter(state == "effective_r_seasonal" & region == do_city) %>%
      pull(state_value)
    
    smooth_reff <- spaero::get_stats(
      x = log(city_reff),
      center_trend = "local_constant", 
      center_kernel = "uniform", 
      center_bandwidth = window_bandwidth, 
      stat_trend = "local_constant", 
      stat_kernel = "uniform", 
      stat_bandwidth = window_bandwidth, 
      lag = 1, 
      backward_only = TRUE
    )$stats$mean
    
    lag <- 4
    smooth_reff2 <- numeric(length(smooth_reff)+lag)
    smooth_reff2[] <- NA
    smooth_reff2[(lag+1):length(smooth_reff2)] <- smooth_reff
    
    city_data <- measles_data %>%
      filter(region == do_city)
    
    city_stats <- spaero::get_stats(
      x = city_data$cases,
      center_trend = "local_constant", 
      center_kernel = "uniform", 
      center_bandwidth = window_bandwidth, 
      stat_trend = "local_constant", 
      stat_kernel = "uniform", 
      stat_bandwidth = window_bandwidth, 
      lag = 1, 
      backward_only = TRUE
    )$stats
    
    city_stats_tb <- as_tibble(city_stats) %>%
      mutate(
        time_iter = 1:n(),
        date = city_data$date,
        state_value = smooth_reff2[1:length(smooth_reff)]
      ) %>%
      gather(key = ews, value = value, -time_iter, -date, -state_value) %>%
      mutate(region = do_city)
    
    all_stats <- bind_rows(all_stats, city_stats_tb)
  }
  
  all_stats <- all_stats %>%
    rename(ews_value = value)
  
  # ews_states <- all_stats %>%
  #   left_join(filtered_states, by = c("date", "region"))
  
  ews_state_corrs <- all_stats %>%
    group_by(ews, region) %>%
    summarise(
      spearman_value = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[1],
      spearman_lwr = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[2],
      spearman_upr = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[3], 
      spearman_pvalue = cor.test(ews_value, state_value, use = "pairwise.complete.obs", method = "spearman")[["p.value"]]
    ) %>%
    mutate(
      pos = spearman_value > 0,
      sig = spearman_pvalue < 0.05,
      color_id_final = "cnull",
      color_id_final = ifelse(pos == TRUE & sig == TRUE, "apos", color_id_final),
      color_id_final = ifelse(pos == FALSE & sig == TRUE, "bneg", color_id_final)
    ) %>%
    mutate(
      bandwidth = do_bw
    )
  
  correlations <- bind_rows(correlations, ews_state_corrs)
  
}

re_seasonal_corrs <- correlations

ggplot(re_seasonal_corrs, aes(x = region, y = ews, fill = color_id_final)) +
  geom_tile(color = "white") +
  facet_wrap(~as.factor(bandwidth), nrow = 1) +
  scale_fill_manual(values = c(ptol_pal()(2)[1], ptol_pal()(2)[2], "grey35"))

