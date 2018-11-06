# plot-seasonality.R
#  Script to plot the MLE transmission seasonality curves for each city.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(pomp)


# Load MLEs ---------------------------------------------------------------

result_files <- list.files("../results/")
mle_files <- result_files[grep("initial-mif-lls", result_files)]

b_splines <- tibble()
betas <- tibble()
for(do_file in mle_files){
  tmp <- read.csv(paste0("../results/", do_file)) %>%
    drop_na() %>%
    filter(loglik == max(loglik)) %>%
    dplyr::select(b1, b2, b3, b4, b5, b6) %>%
    gather() %>%
    mutate(city = str_sub(do_file, 17, (nchar(do_file)-4)))
  
  b_splines <- bind_rows(b_splines, tmp)
  
  tmp2 <- read.csv(paste0("../results/", do_file)) %>%
    drop_na() %>%
    filter(loglik == max(loglik)) %>%
    dplyr::select(beta_mu) %>%
    gather() %>%
    mutate(city = str_sub(do_file, 17, (nchar(do_file)-4)))
  
  betas <- bind_rows(betas, tmp2)
}


# Load exemplar bases from pomp model -------------------------------------

measles_pomp <- readRDS("measles-pomp-object-Agadez.RDS")

bases <- as_tibble(measles_pomp@covar) %>%
  dplyr::select(starts_with("x")) %>%
  dplyr::slice(1:365) %>%
  mutate(
    day = 1:365
  ) %>%
  gather(key = base, value = value, -day)


# Calculate seasonality function ------------------------------------------

season <- bases %>%
  spread(key = base, value = value) %>%
  dplyr::select(-day) %>%
  as.matrix()

seasonal_beta <- tibble()
for(do_city in unique(b_splines$city)){
  qis <- b_splines %>%
    filter(city == do_city) %>%
    pull(value) 
  
  beta_tmp <- betas %>%
    filter(city == do_city) %>%
    pull(value)
  
  seasonal_tmp <- tibble(
    beta = as.numeric((1+exp(season %*% qis)) * beta_tmp),
    day = 1:365,
    city = do_city
  )
  
  seasonal_beta <- bind_rows(seasonal_beta, seasonal_tmp)
}



ggplot(seasonal_beta, aes(x = day, y = beta, color = city)) +
  geom_line(size = 1) +
  labs(x = "Day of year", y = expression(paste("Tranmission rate, ",beta, " (", yr^-1, ")"))) +
  scale_color_manual(values = ggthemes::ptol_pal()(4), name = NULL) +
  theme_minimal() +
  theme(legend.position = c(0.7, 0.8), legend.box.background = element_rect(color = "white"))

seasonal_beta %>%
  group_by(city) %>%
  summarise(minb = min(beta))
