# calculate-boot-cis.R
#  Script to calculate approximate 95% confidence intervals for all parameters
#  from bootstrapped fits of model simulations.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(dplyr)


# Define computation grid -------------------------------------------------

nboots <- 100
nmifs <- 50
comp_grid <- expand.grid(1:nboots, 1:50)
colnames(comp_grid) <- c("boot_series", "param_set")
comp_grid$do_grid <- 1:nrow(comp_grid)


# Load loglikelihood results and combine ----------------------------------

loglik_files <- list.files("../results/", pattern = "bootstrap-mif-lls")
param_summaries <- tibble()

for(do_file in loglik_files){
  do_city <- str_sub(do_file, 19, -5)
  tmp_file <- read.csv(paste0("../results/", do_file)) %>%
    slice(2:n()) %>%  # ignore first row of NAs
    left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
    group_by(boot_series) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    dplyr::select(-do_grid, -loglik, -loglik_se, -boot_series, -param_set,
                  -b1, -b2, -b3, -b4, -b5, -b6) %>%
    gather(key = parameter, value = mle_value) %>%
    group_by(parameter) %>%
    summarise(
      mean_value = mean(mle_value),
      median_value = median(mle_value),
      std_dev = sd(mle_value),
      upper_95_ci = quantile(mle_value, 0.975),
      lower_95_ci = quantile(mle_value, 0.025)
    ) %>%
    mutate(city = do_city)
  
  param_summaries <- bind_rows(param_summaries, tmp_file)
}

# Load example pomp file for basis function
pomp_file <- paste0("measles-pomp-object-", do_city, ".RDS")
measles_pomp <- readRDS(pomp_file)  # exemplar bases


# Calculate seasonal transmission functions -------------------------------

bases <- as_tibble(measles_pomp@covar) %>%
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

seasonal_functions <- tibble()
for(do_file in loglik_files){
  do_city <- str_sub(do_file, 19, -5)
  
  b_splines <- read.csv(paste0("../results/", do_file)) %>%
    slice(2:n()) %>%  # ignore first row of NAs
    left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
    group_by(boot_series) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    dplyr::select(b1, b2, b3, b4, b5, b6)
  
  betas <- read.csv(paste0("../results/", do_file)) %>%
    slice(2:n()) %>%  # ignore first row of NAs
    left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
    group_by(boot_series) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    dplyr::select(beta_mu)
  
  seasonal_betas <- tibble()
  for(i in 1:nrow(b_splines)){
    qis <- as.numeric(b_splines[i, ])
    beta_tmp <- as.numeric(betas[i, "beta_mu"])
    
    seasonal_tmp <- tibble(
      beta = as.numeric((1+exp(season %*% qis)) * beta_tmp),
      day = 1:365,
      boot = i
    )
    seasonal_betas <- bind_rows(seasonal_betas, seasonal_tmp)
  }  # end bootstrap loop
  
  seasonal_betas <- seasonal_betas %>%
    mutate(city = do_city)
  
  seasonal_functions <- bind_rows(seasonal_functions, seasonal_betas)
  
}  # end city loop


# Plot all seasonal transmission curves -----------------------------------

seasonal_functions <- seasonal_functions %>%
  mutate(
    date = as.Date(day, origin = "2016-12-31", tz = "UTC")
  )

ggplot(seasonal_functions, aes(x = date, y = beta, group = boot)) +
  geom_line(alpha = 0.3, color = "dodgerblue4") +
  labs(x = "Date", y = expression(paste("Tranmission rate, ",beta, " (", yr^-1, ")"))) +
  scale_x_date(date_labels = "%b", date_breaks = "2 months") +
  facet_wrap(~city) +
  theme_minimal()
