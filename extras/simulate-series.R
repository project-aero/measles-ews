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
    years_to_sim = 1000, 
    initial_population_size = fitted_pomp@covar[1, "N"], 
    susc_discount = 1  # start at estimated susceptible level
  )
  
  model_sims <- simulate(
    simulator_pomp,
    nsim = 1,
    as.data.frame = TRUE,
    include.data = FALSE,
    seed = 18431236) %>%
    as_tibble()
  
  tmp_re_sims <- model_sims %>%
    dplyr::select(sim, time, RE_seas, reports, S) %>%
    mutate(city = do_city) %>%
    nest(-city)
  
  all_sims <- bind_rows(all_sims, tmp_re_sims)
  
  outfile <- paste0("../simulations/mle-sims-", do_city, ".RDS")
  saveRDS(object = filter(all_sims, city == do_city) %>% unnest(), file = outfile)
}







#### OLD, ONLY OVER DATA ###

# for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
#   # Load pomp object --------------------------------------------------------
#   
#   pomp_file <- paste0("measles-pomp-object-", do_city, ".RDS")
#   measles_pomp <- readRDS(pomp_file)
#   
#   
#   # Load MLEs ---------------------------------------------------------------
#   
#   mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
#   mles <- read.csv(mle_file) %>%
#     as_tibble() %>%
#     filter(loglik == max(loglik, na.rm = TRUE)) %>%
#     dplyr::select(-do_grid, -loglik, -loglik_se)
#   
#   
#   # Simulate time series of incidence ---------------------------------------
#   n_sims <- 1000
#   
#   sims <- as_tibble(
#     simulate(
#       measles_pomp,
#       params = mles,
#       nsim = n_sims,
#       as.data.frame = TRUE,
#       include.data = FALSE,
#       seed = 12345872
#     )
#   ) %>%
#     filter(time >= 1995) %>%
#     dplyr::select(sim, time, reports, RE_seas) %>%
#     mutate(
#       sim = as.numeric(sim)
#     ) %>%
#     nest(-sim)
#   
#   out_file <- paste0("../simulations/mle-sims-", do_city, ".RDS")
#   saveRDS(object = sims, file = out_file)
# }

