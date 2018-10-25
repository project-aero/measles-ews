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
library(DescTools)


# Load and stack simulations ----------------------------------------------

sim_files <- list.files(path = "../simulations/")
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

window_bandwidth <- 26
cities <- all_sims$city

all_corrs <- tibble()

for(do_city in cities){
  tmp_city <- all_sims %>%
    filter(city == do_city) %>%
    unnest() %>%
    ungroup()
  
  for(i in 1:max(tmp_city$sim)){
    tmp_sim <- tmp_city %>%
      filter(sim == i)
    
    tmp_ews <- spaero::get_stats(
      x = tmp_sim$reports,
      center_trend = "local_constant", 
      center_kernel = "uniform", 
      center_bandwidth = window_bandwidth, 
      stat_trend = "local_constant", 
      stat_kernel = "uniform", 
      stat_bandwidth = window_bandwidth, 
      lag = 1, 
      backward_only = TRUE
    )$stats %>%
      as_tibble()
    
    tmp_re <- spaero::get_stats(
      x = log(tmp_sim$RE_seas),
      center_trend = "local_constant", 
      center_kernel = "uniform", 
      center_bandwidth = window_bandwidth, 
      stat_trend = "local_constant", 
      stat_kernel = "uniform", 
      stat_bandwidth = window_bandwidth, 
      lag = 1, 
      backward_only = TRUE
    )$stats$mean
    
    tmp_corr <- cor(
      tmp_re, tmp_ews, 
      use = "pairwise.complete.obs",
      method = "spearman"
    )
    
    tmp_out <- tibble(
      city = do_city,
      sim = i,
      ews = names(tmp_ews),
      spearman_value = as.numeric(tmp_corr)
    )
    
    all_corrs <- bind_rows(all_corrs, tmp_out)
    print(i)
  }  # next sim, i
}  # next city, do_city


outfile <- paste0("../results/sim-corrs-ews-re.RDS")
saveRDS(object = all_corrs, file = outfile)

# ggplot(all_corrs, aes(x = ews, y = spearman_value)) +
#   geom_boxplot(outlier.size = 0.5) +
#   geom_hline(aes(yintercept = 0), color = "red") +
#   facet_wrap(~city, ncol = 1) +
#   labs(y = expression(paste("Spearman's ",rho)), x = "Early warning signal") +
#   coord_flip()

