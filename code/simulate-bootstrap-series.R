# simulate-bootstrap-series.R
#  Script to simulate replicate time series from the MLEs for each
#  city. These bootstrapped time series will be used to re-fit to models
#  to approximate confidence intervals for parameter estimates.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(pomp)
library(tidyverse)
library(lubridate)


# Simulate series from MLEs -----------------------------------------------

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  # Load pomp object 
  pomp_file <- paste0("measles-pomp-object-", do_city, ".RDS")
  measles_pomp <- readRDS(pomp_file)

  # Load MLEs
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>%
    as_tibble() %>%
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)

  # Simulate time series of incidence
  n_boots <- 100

  sims <- as_tibble(
    simulate(
      measles_pomp,
      params = mles,
      nsim = n_boots,
      as.data.frame = TRUE,
      include.data = FALSE,
      seed = 12345872
    )
  ) %>%
    filter(time >= 1995) %>%
    dplyr::select(sim, time, reports, S, E, I, RE_seas) %>%
    mutate(
      sim = as.numeric(sim)
    ) %>%
    nest(-sim)

  out_file <- paste0("../simulations/bootstrap-sims-", do_city, ".RDS")
  saveRDS(object = sims, file = out_file)
}
