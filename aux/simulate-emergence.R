# simulate-emergence.R
#  Script to simulate the (re)emergence of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing increasing effective reproduction ratio, via transmission rate.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(spaero)
source("make-pomp-simulator-function.R")


all_sims <- tibble()

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model -----------------------------------
  
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se, -beta_mu)
  
  mle_beta <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    pull(beta_mu)
  
  pomp_file <- paste0("./measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  # Simulate from the new pomp object ---------------------------------------
  
  simulator_pomp <- make_pomp_simulator(
    do_city, 
    mles, 
    years_to_sim = 30, 
    initial_population_size = fitted_pomp@covar[1, "N"], 
    susc_discount = 0.0000001  # super big discount to S to simulate emergence
  )
  
  model_sims <- simulate(
    simulator_pomp,
    nsim = 500,
    as.data.frame = TRUE,
    include.data = FALSE,
    seed = 1843936) %>%
    as_tibble()
  
  tmp_re_sims <- model_sims %>%
    dplyr::select(sim, time, RE_seas, reports) %>%
    mutate(city = do_city) %>%
    nest(-city)
  
  all_sims <- bind_rows(all_sims, tmp_re_sims)
  
  outfile <- paste0("../simulations/emergence-simulations-", do_city, ".RDS")
  saveRDS(object = filter(all_sims, city == do_city) %>% unnest(), file = outfile)
}

