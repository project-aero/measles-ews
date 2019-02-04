# find-beta-random-walk-intensity.R
#  This script runs a particle filter at the MLEs with a Gamma random walk
#  added to the transmission rate at each time step. We look at the
#  log-likelihood across several RW intensities to find the optimal one.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(lubridate)
source("make-pomp-filtering-function.R")


# Loop over cities and perform particle filter ----------------------------

cities <- c("Agadez", "Maradi", "Niamey", "Zinder")

for(do_city in cities){
  # Load MLEs
  mles <- read.csv(paste0("../results/initial-mif-lls-", do_city, ".csv")) %>%
    filter(loglik == max(loglik, na.rm = T)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  # Make data table
  do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
  measles_data <- readRDS(do_file) %>%
    dplyr::filter(region == paste(do_city, "(City)")) 
  
  obs_data <- measles_data %>%
    dplyr::select(time, cases) %>%
    dplyr::rename(reports = cases)
  
  if(do_city == "Niamey"){
    # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
    obs_data$reports[275] <- round(mean(obs_data$reports[c(274,276)]))
  }
  
  # Make covariate table
  covar_file <- "../data/clean-data/annual-demographic-data-niger-cities-clean.RDS"
  covar_data <- readRDS(covar_file) %>%
    dplyr::filter(region == paste(do_city, "(City)")) %>%
    dplyr::select(time, population_size, birth_per_person_per_year) %>%
    dplyr::rename(
      N = population_size,
      mu = birth_per_person_per_year
    )
  
  # Generate basis functions for seasonality
  bspline_basis <- periodic.bspline.basis(
    covar_data$time,
    nbasis = 6,
    degree = 3,
    period = 1,
    names = "xi%d"
  ) %>%
    as_tibble()
  
  covar_data <- bind_cols(covar_data, bspline_basis)
  
  test_vals <- pretty(seq(0.01, 0.0001, length.out = 30), n = 20)[-1]
  llcity <- tibble()
  counter <- 1
  for(rwval in test_vals){
    measles_pomp <- make_pomp_filter(obs_data, covar_data, mles, rw_value = rwval)
    
    # Calculate log-likelihood of “new” model and save
    ll <- logmeanexp(
      replicate(n = 5,
                logLik(object = pfilter(object = measles_pomp, Np = 5000))), 
      se=FALSE
    )
    
    llcity <- bind_rows(llcity, tibble(loglik = ll, rw_value = rwval))
    print(counter)
    counter <- counter+1
  }  # end rw_value loop
  
  llcity <- llcity %>%
    mutate(city = do_city)
  
  ll_all <- bind_rows(ll_all, llcity)
}  # end city loop

