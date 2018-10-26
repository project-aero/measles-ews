# analyze-emergence-sims.R
#  Script calculate EWS over two windows from simulated series for each
#  city. Simulations start at low number of susceptibles to mimic emergence.
#
# Author: Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(pROC)
library(spaero)


# Load simulations --------------------------------------------------------

all_sims <- readRDS("../simulations/emergence-simulations.RDS") %>%
  unnest()


# Calculate yearly average RE ---------------------------------------------

re_time_avg <- all_sims %>%
  group_by(city, time) %>%
  summarise(mean_re = mean(RE_seas)) %>%
  mutate(
    year = round(time)
  ) %>%
  ungroup() %>%
  group_by(city, year) %>%
  summarise(time_mean_re = mean(mean_re)) %>%
  mutate(
    diff_one = 1 - time_mean_re
  ) 

re_one_year <- re_time_avg %>%
  ungroup() %>%
  group_by(city) %>%
  filter(diff_one == min(diff_one)) %>%
  dplyr::select(city, year) %>%
  ungroup()


# Plot RE series ----------------------------------------------------------

re_series <- ggplot(re_time_avg, aes(x = year, y = time_mean_re)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_segment(data = re_one_year, aes(x = year, xend = year, y = 0, yend = 1), color = ptol_pal()(2)[2]) +
  geom_line(data = all_sims, aes(x = time, y = RE_seas, group = sim), alpha = 0.05, color = "grey") +
  geom_line(size = 0.4, color = ptol_pal()(2)[1]) +
  geom_point(size = 1, color = ptol_pal()(2)[1]) +
  scale_x_continuous(breaks = 1995:2006) +
  labs(x = "Date", y = expression(R[E](t))) +
  facet_wrap(~city, ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "../figures/re-simulated-series.pdf", plot = re_series, width = 3, height = 5, units = "in")




# Format data for 2-window EWS --------------------------------------------

data_for_ews <- tibble()
bws <- tibble()

for(do_city in unique(all_sims$city)){
  tmpsims <- all_sims %>%
    filter(city == do_city)
  
  tmpyear <- re_one_year %>%
    filter(city == do_city) %>%
    pull(year)
  
  tmptimes <- tmpsims %>%
    filter(time > 1994.99) %>%
    filter(round(time) < tmpyear) %>%
    pull(time) %>%
    unique()
  
  if((length(tmptimes) %% 2) != 0){
    tmptimes <- tmptimes[2:length(tmptimes)]
  }
  
  ews_time_ids <- tibble(
    time = tmptimes,
    half = c(rep("first", length(tmptimes)/2), rep("second", length(tmptimes)/2))
  )
  
  window_bandwidth <- length(tmptimes)/2
  
  tmp_ews_data <- tmpsims %>%
    dplyr::select(city, time, reports, sim) %>%
    left_join(ews_time_ids, by = "time") %>%
    filter(is.na(half) == FALSE)
  
  data_for_ews <- bind_rows(tmp_ews_data, data_for_ews)
  bws <- bind_rows(
    bws, 
    tibble(
      city = do_city,
      bandwidth = window_bandwidth
    )
  )
}


# Calculate EWS -----------------------------------------------------------

ews_out <- {}

for(do_city in unique(data_for_ews$city)){
  
  city_data <- data_for_ews %>%
    filter(city == do_city)
  
  window_bandwidth <- bws %>%
    filter(city == do_city) %>%
    pull(bandwidth)
  
  for(i in unique(data_for_ews$sim)){
    
    tmp_data <- city_data %>%
      filter(sim == i)
    
    tmp_first <- spaero::get_stats(
      x = pull(filter(tmp_data, half == "first"), reports),
      center_trend = "local_constant",
      center_kernel = "uniform",
      center_bandwidth = window_bandwidth,
      stat_trend = "local_constant",
      stat_kernel = "uniform",
      stat_bandwidth = window_bandwidth,
      lag = 1,
      backward_only = FALSE)$stats
    
    tmp_second <- spaero::get_stats(
      x = pull(filter(tmp_data, half == "second"), reports),
      center_trend = "local_constant",
      center_kernel = "uniform",
      center_bandwidth = window_bandwidth,
      stat_trend = "local_constant",
      stat_kernel = "uniform",
      stat_bandwidth = window_bandwidth,
      lag = 1,
      backward_only = FALSE)$stats
    
    tmp_out1 <- as_tibble(tmp_first) %>%
      summarise_all(funs(mean, .args = list(na.rm = TRUE))) %>%
      mutate(
        half = "first",
        sim = i
      )
    
    tmp_out2<- as_tibble(tmp_second) %>%
      summarise_all(funs(mean, .args = list(na.rm = TRUE))) %>%
      mutate(
        half = "second",
        sim = i
      )
    
    tmp_out <- bind_rows(tmp_out1, tmp_out2) %>%
      mutate(city = do_city)
    
    ews_out <- bind_rows(ews_out, tmp_out)
    
  }  # end simulation loop
  
}  # end city loop


# Plot the results --------------------------------------------------------

scale_it <- function(x){
  (x-min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}

ews_long <- ews_out %>%
  gather(key = metric, value = value, -half, -sim, -city) %>%
  group_by(city, metric) %>%
  mutate(
    scaled_value = scale_it(value)
  ) %>%
  filter(metric != "variance_first_diff") %>%
  filter(value < 100)

gout <- list()
for(do_city in sort(unique(ews_long$city))){
  if(do_city != "Ziner"){
    gout[[do_city]] <- ggplot(filter(ews_long, city == do_city), aes(fill = half, x = value)) +
      geom_histogram(bins = 20) +
      facet_wrap(~metric, scales = "free", nrow = 1) +
      scale_fill_manual(values = ptol_pal()(2)) +
      theme_minimal(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = FALSE) +
      labs(y = "Count", x = "") +
      ggtitle(do_city)
  }
  
  if(do_city == "Zinder"){
    gout[[do_city]] <- ggplot(filter(ews_long, city == do_city), aes(fill = half, x = value)) +
      geom_histogram(bins = 20) +
      facet_wrap(~metric, scales = "free", nrow = 1) +
      scale_fill_manual(values = ptol_pal()(2)) +
      theme_minimal(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = FALSE) +
      labs(y = "Count", x = "EWS value") +
      ggtitle(do_city)
  }
  
}

ews_hists <- cowplot::plot_grid(plotlist = gout, align = "v", nrow = 4, labels = "AUTO")
ggsave(filename = "../figures/ews-histograms-simulation.pdf", plot = ews_hists, height = 5, width = 9, units = "in")

