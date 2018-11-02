# ews-re-simulation-correlations.R
#  Script to calculate the correlation between R(E) and EWS across many
#  simulations from the fitted models.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(spaero)
library(pROC)


# Define function for grouping time series based on Re --------------------

get_groups <- function(x){
  x <- trunc(x)
  group <- 1
  group_vec <- numeric(length(x))
  group_vec[1] <- group
  for(i in 2:length(x)){
    test <- x[i] < x[i-1]
    if(test == FALSE){group_vec[i] <- group}
    if(test == TRUE){group <- group+1; group_vec[i] <- group}
  }
  return(group_vec)
}


# Load and stack simulations ----------------------------------------------

all_files <- list.files(path = "../simulations/")
sim_files <- all_files[grep("mle", all_files)]
all_sims <- tibble()
for(do_file in sim_files){
  city <- str_sub(string = do_file, start = 10, end = 15)
  
  tmp <- readRDS(paste0("../simulations/", do_file)) %>%
    unnest() %>%
    mutate(
      city = city
    )
  
  all_sims <- bind_rows(all_sims, tmp)
}

all_sims <- all_sims %>%
  nest(-city)


# Calculate EWS and smoothed effective R ----------------------------------

window_bandwidth <- 52
cities <- all_sims$city

ews_out <- tibble()

for(do_city in cities){
  tmp_city <- all_sims %>%
    filter(city == do_city) %>%
    unnest() %>%
    ungroup() %>%
    filter(time != 0)
  
  for(i in 1:1){
    tmp_sim <- tmp_city %>%
      filter(sim == i)
    
    tmp_mean <- tmp_sim %>%
      mutate(year = trunc(time)) %>%
      mutate(
        above = ifelse(RE_seas >= 1, 1, 0)
      ) %>%
      group_by(year) %>%
      summarise(
        sum_re = sum(above)
      ) %>%
      mutate(
        trunc_re = ifelse(sum_re > (52/2), 1, 0)
      )
    
    group_vec <- get_groups(tmp_mean$trunc_re)
    
    year_groups <- tibble(
      year = tmp_mean$year,
      group_id = as.factor(group_vec),
      above = tmp_mean$trunc_re
    )
    
    tmp_sim <- tmp_sim %>%
      mutate(year = trunc(time)) %>%
      left_join(year_groups, by = "year") %>%
      filter(above == 0) %>%
      ungroup() %>%
      group_by(group_id) %>%
      slice(105:n()) %>%
      ungroup()
    
    # ggplot(filter(tmp_sim, time < 301), aes(x = time, y = reports, color = group_id)) +
      # geom_line()
    
    
    for(do_group in unique(tmp_sim$group_id)){
      tmp_grp <- tmp_sim %>%
        filter(group_id == do_group)
      
      if(nrow(tmp_grp) > 104){
        tmp_times <- tmp_grp %>%
          pull(time) %>%
          unique()
        
        if((length(tmp_times) %% 2) != 0){
          tmp_times <- tmp_times[2:length(tmp_times)]
        }
        
        ews_time_ids <- tibble(
          time = tmp_times,
          half = c(
            rep("first", length(tmp_times)/2), 
            rep("second", length(tmp_times)/2)
          )
        )
        
        window_bandwidth <- length(tmp_times)/2
        
        tmp_data <- tmp_grp %>%
          dplyr::select(city, time, RE_seas, reports, sim) %>%
          left_join(ews_time_ids, by = "time") %>%
          filter(is.na(half) == FALSE)
        
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
            sim = i,
            group = do_group
          )
        
        tmp_out2 <- as_tibble(tmp_second) %>%
          summarise_all(funs(mean, .args = list(na.rm = TRUE))) %>%
          mutate(
            half = "second",
            sim = i,
            group = do_group
          )
        
        tmp_out <- bind_rows(tmp_out1, tmp_out2) %>%
          mutate(city = do_city)
        
        ews_out <- bind_rows(ews_out, tmp_out)
      }
      
      # plot(tmp_grp$reports, type = "l")
      # lines(tmp_grp$RE_seas*100, col = "red")
      
      
    }  # next group
    print(i)
  }  # next sim, i
}  # next city, do_city



ews_long <- ews_out %>%
  gather(key = metric, value = value, -half, -sim, -city, -group) %>%
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

ggplot(filter(ews_long, city == "Niamey"), aes(fill = half, x = value)) +
  geom_histogram(bins = 20) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  scale_fill_manual(values = ggthemes::ptol_pal()(2)) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Count", x = "") +
  ggtitle(do_city)



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

write.csv(x = auc_tbl, "../results/endemic-aucs.csv")


# ggplot(auc_tbl, aes(x = metric, y = abs(AUC-0.5), fill = AUC)) +
#   geom_col(position = position_dodge()) +
#   scale_y_continuous(limits = c(0,0.5)) +
#   scale_fill_viridis_c(limits = c(0,1), direction = -1, option = "C") +
#   facet_wrap(~city, nrow = 1) +
#   theme_minimal() +
#   labs(x = NULL, y = "|AUC - 0.5|")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(panel.spacing = unit(1, "lines"))





  

