# estimate-transmission-state.R
#  This script runs a plain vanilla particle filter at the MLEs but allowing
#  mean transmission to talk a random walk in time. The filtering distribution
#  for the state of mean transmission at each time gives us the estimate and
#  distribution of transmission over the course of the time series.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(lubridate)
source("make-pomp-filtering-function.R")


# Functions to extract and summarize filtered states ----------------------

extract_from_list <- function(filter_list, state_number){
  out <- as_tibble(
    lapply(filter_list, `[`, state_number, )
  ) %>%
    mutate(
      particle = 1:n()
    ) %>%
    gather(key = time, value = value, -particle)
  
  return(out)
}

summarise_filtered_state <- function(df, state_name, observations = NA){
  out <- df %>%
    group_by(time) %>%
    summarise(
      med = median(value),
      upper_95 = quantile(value, 0.975),
      lower_95 = quantile(value, 0.025),
      upper_80 = quantile(value, 0.9),
      lower_80 = quantile(value, 0.1)
    ) %>%
    slice(2:n()) %>%
    mutate(
      week = 1:n(),
      date = measles_data$date[2:nrow(measles_data)],
      state = state_name,
      observation = observations
    )
  
  return(out)
}


# Loop over cities and perform particle filter ----------------------------

cities <- c("Agadez", "Maradi", "Niamey", "Zinder")

for(do_city in cities){
  # Load MLEs
  mles <- read.csv(paste0("../results/initial-mif-lls-", do_city, ".csv")) %>%
    filter(loglik == max(loglik, na.rm = T)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  # Make data table
  do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
  obs_data <- readRDS(do_file) %>%
    dplyr::filter(region == paste(do_city, "(City)")) %>%
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
  
  measles_pomp <- make_pomp_filter(obs_data, covar_data, mles)
  
  
  # Run particle filter
  measles_filter <- pfilter(object = measles_pomp, Np = 50000, 
                            save.states = TRUE, pred.mean = TRUE, 
                            pred.var = TRUE)
  states <- measles_filter@saved.states  # filtering distribution of states
  predicted_means <- measles_filter@pred.mean  # prediction distribution of states
  predicted_vars <- measles_filter@pred.var  # prediction variance of states
  
  # Extract prediction means and variances for cases
  pred_cases <- t(predicted_means)[,"cases_state"]
  pred_stdev <- sqrt(t(predicted_vars)[,"cases_state"])
  
  predictive_distribution <- tibble(
    mean_cases = pred_cases,
    sdev_cases = pred_stdev
  )
  
  pred_outfile <- paste0("../results/predictive-dist-states-", DO_CITY, ".RDS")
  saveRDS(
    object = predictive_distribution, 
    file = pred_outfile
  )
  
  
  # Extract and summarize states
  state_names <- names(states[[1]][,1])
  
  eff_rep_nonseasonal <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "RE")
  ) %>%
    summarise_filtered_state(
      state_name = "effective_r_nonseasonal",
      observations = NA
    )
  
  eff_rep_seasonal <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "RE_seas")
  ) %>%
    summarise_filtered_state(
      state_name = "effective_r_seasonal",
      observations = NA
    )
  
  transmission_rate <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "beta_t")
  ) %>%
    summarise_filtered_state(
      state_name = "transmission_rate",
      observations = NA
    )
  
  susceptibles <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "S")
  ) %>%
    summarise_filtered_state(
      state_name = "susceptibles",
      observations = NA
    )
  
  exposed <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "E")
  ) %>%
    summarise_filtered_state(
      state_name = "exposed",
      observations = NA
    )
  
  infected <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "I")
  ) %>%
    summarise_filtered_state(
      state_name = "infected",
      observations = NA
    )
  
  cases <- extract_from_list(
    filter_list =  states, 
    state_number = which(state_names == "cases_state")
  ) %>%
    summarise_filtered_state(
      state_name = "cases",
      observations = measles_data$cases[2:nrow(measles_data)]
    )
  
  
  # Combine state dfs and save ----------------------------------------------
  
  all_states <- eff_rep_nonseasonal %>%
    bind_rows(eff_rep_seasonal) %>%
    bind_rows(transmission_rate) %>%
    bind_rows(susceptibles) %>%
    bind_rows(exposed) %>%
    bind_rows(infected) %>%
    bind_rows(cases) %>%
    ungroup() %>%
    group_by(state) %>%
    nest()
  
  outfile <- paste0("../results/filtered-states-", do_city, ".RDS")
  saveRDS(
    object = all_states, 
    file = outfile
  )
  
  
  # Calculate log-likelihood of “new” model and save
  
  ll <- logmeanexp(
    replicate(
      n = 10,
      logLik(
        object = pfilter(
          object = measles_pomp, 
          Np = 10000
        )
      )
    ), 
    se=TRUE
  )
  
  ll_tbl <- ll %>%
    as_tibble() %>%
    mutate(
      var = c("loglik", "loglik_se")
    )
  
  ll_file <- paste0("../results/time-vary-beta-loglik-", do_city, ".csv")
  write.csv(x = ll_tbl, file = ll_file)
}

