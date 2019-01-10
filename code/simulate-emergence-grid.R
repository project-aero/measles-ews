# simulate-emergence.R
#  Script to simulate the (re)emergence of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing increasing effective reproduction ratio, via transmission rate.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Get number of reps from command line ------------------------------------

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
nreps <- as.numeric(myargument)

# Set up grid of susceptible depletions -----------------------------------

discount_grid <- c(0.0001, seq(0.1, 1, length.out = 10))


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)  # functions for parallel computing

source("./code/make-pomp-simulator-function.R")

# Set number of cores based on machine ------------------------------------

if(parallel::detectCores() <= length(discount_grid)){
  registerDoParallel(cores = detectCores() - 1)
} else{
  registerDoParallel()
}


all_sims <- tibble()

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model -----------------------------------
  
  mle_file <- paste0("./results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  pomp_file <- paste0("./code/measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  
  # Simulate from the new pomp object ---------------------------------------
  
  outsims <- foreach(i = discount_grid,
                     .packages = c("pomp", "tidyverse", "dplyr"), 
                     .combine = "rbind") %dopar%
  {
    simulator_pomp <- make_pomp_simulator(
      do_city, 
      mles, 
      years_to_sim = 40, 
      initial_population_size = round(mean(fitted_pomp@covar[, "N"])), 
      susc_discount = i,
      vacc_coverage_ts = NULL
    )
    
    model_sims <- simulate(
      simulator_pomp,
      nsim = nreps,
      as.data.frame = TRUE,
      include.data = FALSE) %>%
      as_tibble()
    
    tmp_re_sims <- model_sims %>%
      dplyr::select(sim, time, RE_seas, reports) %>%
      mutate(
        city = do_city,
        susc_discount = i
      )
    
    outfile <- paste0("./simulations/emergence-simulations-grid-", do_city, "-", i, ".RDS")
    saveRDS(object = tmp_re_sims, file = outfile)
  }
  
}
