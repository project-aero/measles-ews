# simulate-emergence.R
#  Script to simulate the (re)emergence of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing increasing effective reproduction ratio, via transmission rate.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Set up grid of susceptible depletions -----------------------------------

discount_grid <- c(0.0001, seq(0.1, 1, length.out = 10))

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(spaero)
library(foreach)
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
  
  outsims <- foreach(i = discount_grid, .packages = c("pomp", "tidyverse", "dplyr"), .combine = "rbind"){
    simulator_pomp <- make_pomp_simulator(
      do_city, 
      mles, 
      years_to_sim = 30, 
      initial_population_size = fitted_pomp@covar[1, "N"], 
      susc_discount = i
    )
    
    model_sims <- simulate(
      simulator_pomp,
      nsim = 500,
      as.data.frame = TRUE,
      include.data = FALSE) %>%
      as_tibble()
    
    tmp_re_sims <- model_sims %>%
      dplyr::select(sim, time, RE_seas, reports) %>%
      mutate(
        city = do_city,
        susc_discount = i
      ) %>%
      nest(-city)
  }

  outfile <- paste0("../simulations/emergence-simulations-", do_city, ".RDS")
  saveRDS(object = filter(all_sims, city == do_city) %>% unnest(), file = outfile)
}

