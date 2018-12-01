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

all_files <- list.files("../simulations/")
sim_file_ids <- grep("emergence-simulations-grid", all_files)
sim_files <- all_files[sim_file_ids]
all_sims <- {}
for(do_file in sim_files){
  tmp_file <- paste0("../simulations/", do_file)
  tmp <- readRDS(tmp_file) %>%
    filter(time > 0 & susc_discount < 0.6)
  all_sims <- bind_rows(all_sims, tmp)
}


# Calculate yearly average RE ---------------------------------------------

re_time_avg <- all_sims %>%
  group_by(city, time, susc_discount) %>%
  summarise(mean_re = mean(RE_seas)) %>%
  mutate(
    year = round(time)
  ) %>%
  ungroup() %>%
  group_by(city, year, susc_discount) %>%
  summarise(time_mean_re = mean(mean_re)) %>%
  filter(year <= 20)

re_one_year <- re_time_avg %>%
  ungroup() %>%
  group_by(city, susc_discount) %>%
  filter(round(time_mean_re,1) >= 1) %>%
  dplyr::select(city, year) %>%
  filter(year == min(year)) %>%
  ungroup()

plot_sims <- all_sims %>% 
  filter(sim < 21 & time <= 20)


# Plot RE series ----------------------------------------------------------

re_series <- ggplot(re_time_avg, aes(x = year, y = time_mean_re)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_segment(
    data = re_one_year,
    aes(x = year, xend = year, y = 0, yend = 1),
    color = ptol_pal()(2)[2]
  ) +
  geom_line(
    data = plot_sims,
    aes(x = time, y = RE_seas, group = sim),
    alpha = 0.1,
    color = "grey"
  ) +
  geom_line(size = 0.4, color = ptol_pal()(2)[1]) +
  geom_point(size = 1, color = ptol_pal()(2)[1]) +
  labs(x = "Simulation year", y = expression(R[E](t))) +
  facet_grid(susc_discount~city) +
  theme_minimal() 

ggsave(
  filename = "../figures/re-simulated-series.pdf", 
  plot = re_series, 
  width = 3, 
  height = 5, 
  units = "in"
)




# Format data for 2-window EWS --------------------------------------------

data_for_ews <- tibble()
bws <- tibble()

for(do_city in unique(all_sims$city)){
  for(do_discount in unique(all_sims$susc_discount)){
    tmpsims <- all_sims %>%
      filter(city == do_city & susc_discount == do_discount)
    
    tmpyear <- re_one_year %>%
      filter(city == do_city & susc_discount == do_discount) %>%
      pull(year)
    
    tmptimes <- tmpsims %>%
      filter(round(time) < tmpyear) %>%
      pull(time) %>%
      unique()
    
    if((length(tmptimes) %% 2) != 0){
      tmptimes <- tmptimes[2:length(tmptimes)]
    }
    
    ews_time_ids <- tibble(
      time = tmptimes,
      half = c(
        rep("first", length(tmptimes)/2), 
        rep("second", length(tmptimes)/2)
      )
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
        susc_discount = do_discount,
        bandwidth = window_bandwidth
      )
    )
  }
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
  # filter(value < 100) %>%
  ungroup() %>%
  mutate(
    metric = ifelse(metric == "variance", "Variance", metric),
    metric = ifelse(metric == "variance_first_diff", "Var. 1st Diff.", metric),
    metric = ifelse(metric == "autocovariance", "Autocovar.", metric),
    metric = ifelse(metric == "autocorrelation", "Autocorr.", metric),
    metric = ifelse(metric == "decay_time", "Decay time", metric),
    metric = ifelse(metric == "mean", "Mean", metric),
    metric = ifelse(metric == "index_of_dispersion", "Index of dis.", metric),
    metric = ifelse(metric == "coefficient_of_variation", "Coeff. var.", metric),
    metric = ifelse(metric == "skewness", "Skewness", metric),
    metric = ifelse(metric == "kurtosis", "Kurtosis", metric)
  )



gout <- list()
for(do_city in sort(unique(ews_long$city))){
  if(do_city != "Zinder"){
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

ews_hists <- cowplot::plot_grid(
  plotlist = gout, 
  align = "v", 
  nrow = 4, 
  labels = "AUTO"
)
ggsave(
  filename = "../figures/ews-histograms-simulation.pdf", 
  plot = ews_hists, 
  height = 5,
  width = 9, 
  units = "in"
)


# Calculate AUC -----------------------------------------------------------

cats <- tibble(
  half = c("first", "second"),
  cat = c(0, 1)
)

ews_long <- ews_long %>%
  left_join(cats, by = "half")

auc_tbl <- {}
for(do_city in unique(ews_long$city)){
  for(do_metric in unique(ews_long$metric)){
    tmp <- filter(ews_long, metric == do_metric & city == do_city)
    roc_obj <- roc(tmp$cat, tmp$value)
    tmp_auc <- auc(roc_obj)
    plot(roc_obj, main = do_metric)
    tmp_tbl <- tibble(
      city = do_city,
      metric = do_metric,
      AUC = as.numeric(tmp_auc)
    )
    
    auc_tbl <- bind_rows(auc_tbl, tmp_tbl)
  }
}

write.csv(x = auc_tbl, "../results/emergence-aucs.csv")


# plt_tbl <- auc_tbl %>%
#   filter(metric %in% c("Variance", "Autocovar.", "Autocorr.",
#                        "Decay time", "Mean"))

# auc_plot <- ggplot(auc_tbl, aes(x = metric, y = AUC-0.5, fill = AUC)) +
#   geom_col(position = position_dodge()) +
#   scale_y_continuous(limits = c(0,0.5)) +
#   scale_fill_viridis_c(limits = c(0,1), direction = -1, option = "C") +
#   facet_wrap(~city, nrow = 1) +
#   theme_minimal() +
#   labs(x = NULL)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(panel.spacing = unit(1, "lines"))
# ggsave(filename = "../figures/sim-emergence-aucs.pdf", plot = auc_plot, width = 8.5, height = 2.8, units = "in")
