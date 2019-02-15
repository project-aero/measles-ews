# analyze-emergence-sims.R
#  Script calculate EWS over two windows from simulated series for each
#  city. Simulations start at low number of susceptibles to mimic emergence.
#
# Author: Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(spaero)


# Load simulations --------------------------------------------------------

# This chunck loops over all the relevant emergence simulation files
# and combines them into one long tibble for futher processing and
# subsetting. We simulated emergence at lower levels of susceptible
# discounting, but focus this analysis on discount factors that are less
# than 0.6. 

all_files <- list.files("../simulations/")
sim_file_ids <- grep("emergence-simulations-grid", all_files)
sim_files <- all_files[sim_file_ids]
ignore_id <- grep("0.6|0.7|0.8|0.9|-1.RDS", sim_files)
sim_files_reduced <- sim_files[-ignore_id]
all_sims_list <- list()
counter <- 1
for(do_file in sim_files_reduced){
  tmp_file <- paste0("../simulations/", do_file)
  tmp <- readRDS(tmp_file) %>%
    filter(time > 0)
  all_sims_list[[counter]] <- tmp
  counter <- counter+1
}

all_sims <- as_tibble(data.table::rbindlist(all_sims_list))
rm(all_sims_list)  # remove the large list from memory


# Find year of critical transition ----------------------------------------

# The year of the critical transition is defined as the year just after
# effective reproduction number (R_E) reaches or exceeds the critical
# value of 1. For example, if R_E reaches or exceeds 1 at some point during
# the fifth year of the simulation, then the critical transition year is
# defined as the sixth year of the simulation. Thus, the window for 
# calculating early warning signals ends at the end of the fifth year. The
# dplyr chain below finds these critical years.

re_one_year <- all_sims %>%
  group_by(city, time, susc_discount) %>%
  summarise(mean_re = mean(RE_seas)) %>%  
  mutate(year = round(time)) %>%  # create 'year' variable
  ungroup() %>%
  group_by(city, susc_discount) %>%
  filter(mean_re >= 1) %>%  # drop times where Re less than 1
  filter(year == min(year) + 1) %>%  # filter to critical transition year
  distinct(city, susc_discount, year, .keep_all = TRUE) %>%  # drop duplicates
  dplyr::select(city, year, susc_discount) %>%
  ungroup()


# re_time_avg <- all_sims %>%
#   group_by(city, time, susc_discount) %>%
#   summarise(mean_re = mean(RE_seas)) %>%
#   mutate(
#     year = round(time)
#   ) %>%
#   ungroup() %>%
#   group_by(city, year, susc_discount) %>%
#   summarise(time_mean_re = mean(mean_re)) %>%
#   filter(year <= 25)
# 
# re_one_year <- re_time_avg %>%
#   ungroup() %>%
#   group_by(city, susc_discount) %>%
#   filter(round(time_mean_re,1) >= 1) %>%
#   dplyr::select(city, year, susc_discount) %>%
#   filter(year == min(year)) %>%
#   ungroup()


# Plot RE series ----------------------------------------------------------

plot_sims <- all_sims %>% 
  filter(sim < 21 & time <= 25)  # just plot 20 reps for 25 years

re_series <- ggplot() +
  geom_hline(data = re_one_year, aes(yintercept = 1), linetype = 2) +
  geom_segment(
    data = re_one_year,
    aes(x = year, xend = year, y = 0, yend = 1),
    color = ptol_pal()(2)[2]
  ) +
  geom_line(
    data = plot_sims,
    aes(x = time, y = RE_seas, group = sim),
    alpha = 0.1,
    color = "grey55",
    size = 0.3
  ) +
  labs(x = "Simulation year", y = expression(R[E](t))) +
  facet_grid(susc_discount~city) +
  theme_minimal() 

# ggsave(
#   filename = "../figures/effective-r-emergence-grid.pdf", 
#   plot = re_series, 
#   width = 8.5, 
#   height = 5, 
#   units = "in"
# )


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
      dplyr::select(city, time, reports, sim, susc_discount) %>%
      left_join(ews_time_ids, by = "time") %>%
      filter(is.na(half) == FALSE)
    
    data_for_ews <- bind_rows(data_for_ews, tmp_ews_data)
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

bws %>%
  arrange(city, susc_discount) %>%
  write.csv(file = "../results/emergence-bandwidths.csv")

# Calculate EWS -----------------------------------------------------------

my_get_stats <- function(x, bw){
  suppressWarnings(
    ews_list <- spaero::get_stats(
      x = x,
      center_trend = "local_constant",
      center_kernel = "uniform",
      center_bandwidth = bw,
      stat_trend = "local_constant",
      stat_kernel = "uniform",
      stat_bandwidth = bw,
      lag = 1,
      backward_only = FALSE)$stats
  )
  
  stats <- sapply(ews_list, mean, na.rm = TRUE)
  return(stats)
}

ews_out <- {}

for(do_city in unique(data_for_ews$city)){
  
  city_data <- data_for_ews %>%
    filter(city == do_city)
  
  for(do_discount in unique(data_for_ews$susc_discount)){
   
    discount_data <- city_data %>%
      filter(susc_discount == do_discount) %>%
      dplyr::select(time, reports, sim, half)
    
    first_half <- discount_data %>%
      filter(half == "first") %>%
      dplyr::select(-half) %>%
      spread(key = sim, value = reports) %>%
      dplyr::select(-time) %>%
      as.matrix()
    
    second_half <- discount_data %>%
      filter(half == "second") %>%
      dplyr::select(-half) %>%
      spread(key = sim, value = reports) %>%
      dplyr::select(-time) %>%
      as.matrix()
    
    window_bandwidth <- bws %>%
      filter(city == do_city & susc_discount == do_discount) %>%
      pull(bandwidth)
    
    ews_first <- apply(X = first_half, MARGIN = 2, FUN = my_get_stats, 
                       bw = window_bandwidth)
    
    ews_second <- apply(X = second_half, MARGIN = 2, FUN = my_get_stats, 
                        bw = window_bandwidth)
    
    first_tbl <- tibble(
      metric = row.names(ews_first)
    ) %>%
      bind_cols(as_tibble(ews_first)) %>%
      gather(key = sim, value = value, - metric) %>%
      mutate(half = "first")
    
    second_tbl <- tibble(
      metric = row.names(ews_second)
    ) %>%
      bind_cols(as_tibble(ews_second)) %>%
      gather(key = sim, value = value, - metric) %>%
      mutate(half = "second")
    
    tmp_out <- bind_rows(first_tbl, second_tbl) %>%
      mutate(
        city = do_city,
        susc_discount = do_discount
      )
    
    ews_out <- bind_rows(ews_out, tmp_out)
    
  }  # end susceptible discount loop
  
}  # end city loop

# Save the results
write.csv(x = ews_out, file = "../results/ews-emergence.csv", row.names = FALSE)



# Reformat results --------------------------------------------------------

col_spec <- cols(
  metric = col_character(),
  sim = col_integer(),
  value = col_double(),
  half = col_character(),
  city = col_character(),
  susc_discount = col_double()
)

ews_long <- read_csv("../results/ews-emergence.csv", col_types = col_spec) %>%
  filter(metric != "variance_first_diff") %>%
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


# Calculate AUC -----------------------------------------------------------

calc_auc <- function(predictions, is_null){
  # Function to calculate AUC, from Eamon.
  # Returns exact same results as those from popular R packages.
  
  r <- rank(predictions)
  r1 <- sum(r[!is_null])
  n1 <- sum(!is_null)
  n2 <- sum(is_null)
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}

cats <- tibble(
  half = c("first", "second"),
  cat = c(TRUE, FALSE)  # null (TRUE), test (FALSE)
)

ews_long <- ews_long %>%
  left_join(cats, by = "half")

auc_tbl <- {}
for(do_city in unique(ews_long$city)){
  for(do_discount in unique(ews_long$susc_discount)){
    for(do_metric in unique(ews_long$metric)){
      tmp <- filter(ews_long, metric == do_metric & city == do_city & susc_discount == do_discount)
      tmp_auc <- calc_auc(predictions = tmp$value, is_null = tmp$cat)
      tmp_tbl <- tibble(
        city = do_city,
        susc_discount = do_discount,
        metric = do_metric,
        AUC = as.numeric(tmp_auc)
      )
      
      auc_tbl <- bind_rows(auc_tbl, tmp_tbl)
    }
  }
}

write.csv(x = auc_tbl, "../results/emergence-grid-aucs.csv")




# OLD CODE ----------------------------------------------------------------

# gout <- list()
# for(do_city in sort(unique(ews_long$city))){
#   if(do_city != "Zinder"){
#     gout[[do_city]] <- ggplot(filter(ews_long, city == do_city), aes(fill = half, x = value)) +
#       geom_histogram(bins = 20) +
#       facet_wrap(~metric, scales = "free", nrow = 1) +
#       scale_fill_manual(values = ptol_pal()(2)) +
#       theme_minimal(base_size = 8) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       guides(fill = FALSE) +
#       labs(y = "Count", x = "") +
#       ggtitle(do_city)
#   }
#   
#   if(do_city == "Zinder"){
#     gout[[do_city]] <- ggplot(filter(ews_long, city == do_city), aes(fill = half, x = value)) +
#       geom_histogram(bins = 20) +
#       facet_wrap(~metric, scales = "free", nrow = 1) +
#       scale_fill_manual(values = ptol_pal()(2)) +
#       theme_minimal(base_size = 8) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       guides(fill = FALSE) +
#       labs(y = "Count", x = "EWS value") +
#       ggtitle(do_city)
#   }
#   
# }
# 
# ews_hists <- cowplot::plot_grid(
#   plotlist = gout, 
#   align = "v", 
#   nrow = 4, 
#   labels = "AUTO"
# )
# ggsave(
#   filename = "../figures/ews-histograms-simulation.pdf", 
#   plot = ews_hists, 
#   height = 5,
#   width = 9, 
#   units = "in"
# )




# ggplot(auc_tbl, aes(x = as.factor(susc_discount), y = metric, fill = abs(AUC-0.5))) +
#   geom_tile() +
#   scale_fill_viridis_c(limits = c(0,0.5), direction = -1, option = "C", name = "| AUC - 0.5 |") +
#   facet_wrap(~city, nrow = 1) +
#   labs(x = "Level of susceptible depletion", y = NULL) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(panel.spacing = unit(1, "lines"))



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
