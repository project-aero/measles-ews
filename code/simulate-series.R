# simulate-series.R
#  Script to simulate several incidence time series from a model with known
#  parameters. These simulation results are used to test the analysis
#  comparing EWS(t) and R_E(t) through time.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(pomp)
library(tidyverse)
library(lubridate)


for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  # Load pomp object --------------------------------------------------------
  
  pomp_file <- paste0("measles-pomp-object-", do_city, ".RDS")
  measles_pomp <- readRDS(pomp_file)
  
  
  # Load MLEs ---------------------------------------------------------------
  
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>%
    as_tibble() %>%
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  
  # Simulate time series of incidence ---------------------------------------
  n_sims <- 1000
  
  sims <- as_tibble(
    simulate(
      measles_pomp,
      params = mles,
      nsim = n_sims,
      as.data.frame = TRUE,
      include.data = FALSE,
      seed = 12345872
    )
  ) %>%
    filter(time >= 1995) %>%
    dplyr::select(sim, time, reports, RE_seas) %>%
    mutate(
      sim = as.numeric(sim)
    ) %>%
    nest(-sim)
  
  out_file <- paste0("../simulations/mle-sims-", do_city, ".RDS")
  saveRDS(object = sims, file = out_file)
}

