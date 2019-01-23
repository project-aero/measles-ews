# analyze-elimination-sims.R
#  Script to calculate EWS over two windows from simulated series for each
#  city. Simulations start at low vaccine coverage and go to high, 
#  mimicking elimination.
#
# Author: Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(pROC)
library(spaero)


# Define function to calculate R0 from seasonal params --------------------

calc_R0 <- function(beta, qis, season, eta = (365/8), 
                    mu = 0.05, nu = 0.05, gamma = (365/5), p = 0.7){
  B <- as.numeric((1 + exp(season %*% qis)) * beta)
  R0 <- (eta*B*mu*(1-p)) / (nu*(eta + nu)*(gamma + nu))
  return(R0)
}


# Calculate vaccination thresholds for each city --------------------------

vacc_thresholds <- {}

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model 
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  pomp_file <- paste0("./measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  # Calculte R0 
  qis <- mles %>%
    dplyr::select(b1, b2, b3, b4, b5, b6) %>%
    as.numeric()
  
  beta <- mles %>%
    pull(beta_mu)
  
  bases <- as_tibble(fitted_pomp@covar) %>%
    dplyr::select(starts_with("x")) %>%
    dplyr::slice(1:365) %>%
    mutate(
      day = 1:365
    ) %>%
    gather(key = base, value = value, -day)
  
  season <- bases %>%
    spread(key = base, value = value) %>%
    dplyr::select(-day) %>%
    as.matrix()
  
  N <- round(mean(fitted_pomp@covar[, "N"]))
  
  R0 <- calc_R0(beta = beta, qis = qis, season = season)
  crit_vacc_cover <- 1 - (1/max(R0))
  vacc_thresholds <- bind_rows(
    vacc_thresholds, 
    tibble(
      city = do_city,
      threshold = crit_vacc_cover
    )
  )
}


# Load simulations --------------------------------------------------------

all_files <- list.files("../simulations/")
sim_file_ids <- grep("elimination-simulations-grid", all_files)
sim_files <- all_files[sim_file_ids]
all_sims_list <- list()
counter <- 1
for(do_file in sim_files){
  tmp_file <- paste0("../simulations/", do_file)
  tmp <- readRDS(tmp_file) %>%
    filter(time > 0) %>%
    left_join(vacc_thresholds, by = "city") %>%
    filter(vacc_coverage <= threshold)
  all_sims_list[[counter]] <- tmp
  counter <- counter+1
}

all_sims <- as_tibble(data.table::rbindlist(all_sims_list))
rm(all_sims_list)

# Format data for 2-window EWS --------------------------------------------

data_for_ews <- tibble()
bws <- tibble()

for(do_city in unique(all_sims$city)){
  for(do_speed in unique(all_sims$vacc_speed)){
    
    tmpsims <- all_sims %>%
      filter(city == do_city) %>%
      filter(vacc_speed == do_speed)
    
    window_size <- tmpsims %>%
      filter(vacc_coverage > 0.7 & sim == 1) %>%
      nrow()  # returns number to times in window leading to tcritical
    
    tmptimes <- tmpsims %>%
      pull(time) %>%
      unique()  # returns unique observation times across simulations
    
    # if((length(tmptimes) %% 2) != 0){
    #   tmptimes <- tmptimes[2:length(tmptimes)]
    # }
    # 
    # ews_time_ids <- tibble(
    #   time = tmptimes,
    #   half = c(
    #     rep("first", length(tmptimes)/2), 
    #     rep("second", length(tmptimes)/2)
    #   )
    # )
    # 
    # window_bandwidth <- length(tmptimes)/2
    
    times_before_vaccine <- tmptimes[tmptimes < 50]
    start_ts <- length(times_before_vaccine) - window_size + 1
    times_before_vaccine <- times_before_vaccine[start_ts:length(times_before_vaccine)]
    times_after_vaccine <- tmptimes[tmptimes >= 50]

    simtimes <- c(times_before_vaccine, times_after_vaccine)

    ews_time_ids <- tibble(
      time = simtimes,
      half = ifelse(
        time < 50,
        "first",
        "second"
      )
    )
    
    # Bandwidth is the full window for each half
    window_bandwidth <- length(simtimes)/2
    
    # Merge in `time` and `half` information
    tmp_ews_data <- tmpsims %>%
      dplyr::select(city, time, reports, sim, vacc_speed) %>%
      left_join(ews_time_ids, by = "time") %>%
      filter(is.na(half) == FALSE)
    
    # Bind for time series storage
    data_for_ews <- bind_rows(data_for_ews, tmp_ews_data)
    
    # Store bandwidths for each city and speed
    bws <- bind_rows(
      bws, 
      tibble(
        city = do_city,
        vacc_speed = do_speed,
        bandwidth = window_bandwidth
      )
    )
    
  }  # end vaccine coverage speed loop
}  # end city loop


# Plot an example series and window ---------------------------------------

ews_data_example <- data_for_ews %>%
  filter(city == "Maradi" & sim == 1 & vacc_speed == 0.000015)

min_time <- min(ews_data_example$time)

example_file <- "../simulations/elimination-simulations-grid-Maradi-1.5e-05.RDS"
example_data <- readRDS(example_file) %>%
  filter(sim == 2) %>%
  left_join(vacc_thresholds, by = "city") %>%
  mutate(
    use_it = ifelse(vacc_coverage > threshold | time < min_time, FALSE, TRUE)
  ) %>%
  filter(use_it == FALSE) %>%
  mutate(group = ifelse(time < 50, 1, 2))

saveRDS(list(example_data, ews_data_example), "../simulations/single-elimination-example.RDS")
bws %>%
  arrange(city, vacc_speed) %>%
  write.csv(file = "../results/elimination-bandwidths.csv")

# ggplot() +
#   geom_line(data = example_data, aes(x = time, y = reports, group = group), color = "tan", alpha = 0.4, size = 0.3) +
#   geom_vline(data = ews_data_example, aes(xintercept = max(time)), color = "dodgerblue4",  linetype = 2) +
#   geom_vline(data = ews_data_example, aes(xintercept = min(time)), color = "dodgerblue4",  linetype = 2) +
#   geom_vline(data = ews_data_example, aes(xintercept = 50), color = "grey45", linetype = 1) +
#   geom_line(data = ews_data_example, aes(x = time, y = reports), color = "tan", size = 0.3) +
#   annotate(geom = "text", x = 28, y = 1450, label = "Null window") +
#   annotate(geom = "text", x = 70, y = 1450, label = "Test window") +
#   annotate(geom = "text", x = 38, y = 900, label = "start of\nvaccination campaign", size = 3, color = "grey25") +
#   annotate(geom = "text", x = 81, y = 900, label = "coverage reaches\nherd immunity", size = 3, color = "dodgerblue4") +
#   labs(x = "Simulation time (year)", y = "Reported cases") +
#   theme_classic(base_size = 14)
# ggsave(filename = "../figures/elimination-series-example.pdf", height = 3, width = 7, units = "in")


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
  
  for(do_speed in unique(data_for_ews$vacc_speed)){
    
    discount_data <- city_data %>%
      filter(vacc_speed == do_speed) %>%
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
      filter(city == do_city & vacc_speed == do_speed) %>%
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
        vacc_speed = do_speed
      )
    
    ews_out <- bind_rows(ews_out, tmp_out)
    
  }  # end susceptible discount loop
  
}  # end city loop

# Save the results
write.csv(x = ews_out, file = "../results/ews-elimination.csv", row.names = FALSE)


# Format results ----------------------------------------------------------

ews_long <- ews_out %>%
  group_by(city, metric) %>%
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


# Calculate AUC -----------------------------------------------------------

cats <- tibble(
  half = c("first", "second"),
  cat = c(0, 1)
)

ews_long <- ews_long %>%
  left_join(cats, by = "half")

auc_tbl <- {}
for(do_city in unique(ews_long$city)){
  for(do_speed in unique(ews_long$vacc_speed)){
    for(do_metric in unique(ews_long$metric)){
      tmp <- filter(ews_long, metric == do_metric & city == do_city & vacc_speed == do_speed)
      roc_obj <- roc(tmp$cat, tmp$value)
      tmp_auc <- auc(roc_obj)
      tmp_tbl <- tibble(
        city = do_city,
        vacc_speed = do_speed,
        metric = do_metric,
        AUC = as.numeric(tmp_auc)
      )
      
      auc_tbl <- bind_rows(auc_tbl, tmp_tbl)
    }
  }
}

write.csv(x = auc_tbl, "../results/elimination-grid-aucs.csv")




# EXTRAS ------------------------------------------------------------------

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
#   filename = "../figures/ews-histograms-elimination.pdf",
#   plot = ews_hists,
#   height = 5,
#   width = 9,
#   units = "in"
# )
# 
# 
# ggplot(auc_tbl, aes(x = as.factor(vacc_speed), y = metric, fill = abs(AUC-0.5))) +
#   geom_tile() +
#   viridis::scale_fill_viridis(limits = c(0,0.5), direction = -1, option = "C", name = "| AUC - 0.5 |") +
#   facet_wrap(~city, nrow = 1) +
#   labs(x = "Speed of vaccination campaign", y = NULL) +
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
