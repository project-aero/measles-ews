# large-mif-search.R
#  R script to be run on a multi-core high performance machine. The script runs
#  through several algorithms for estimating parameters for a measles SIR
#  model fit to data from Niamey, Niger. Fitting is primarily done using
#  maximization of the likelihood via iterated filtering, implemented using the
#  `pomp` package.
#
#  The script runs in three main chunks:
#
#    1.  Initial parameter guesses from probe matching, i.e., fitting to a 
#        synthetic likelihood.
#    2a. A global parameter search from a Latin hypercube of parameter values
#        based on noised-up values from (1).
#    2b. A local parameter search around the MLE based on the highest likelihood
#        parameters from (2a).
#    3.  Particle MCMC for Bayesian inference of latent states, including
#        parameter uncertainty, with values starting at the MLE parameters
#        from (2b).
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)
#  Pejman Rohani


# Load packages -----------------------------------------------------------

library(tidyverse)   # data wrangling
library(lubridate)   # time and date functions
library(pomp)        # fitting state-space models with particle filtering
library(foreach)     # functions for parallel computing
library(doParallel)  # functions for parallel computing

registerDoParallel()
theme_set(theme_minimal())



# Load formatted data and pomp object -------------------------------------

biweek_data <- readRDS("niger_cleaned_biweekly_measles.RDS")
measles_pomp <- readRDS("measles-pomp.RDS")
param_names <- names(coef(measles_pomp))


# Probe model with synthetic likelihood -----------------------------------

probe.zeroes <- function(y){
  # number of zeros in the data
  
  xy <- y["cases", ]
  as.numeric(length(which(xy == 0)))
}
probe.max <- function(y){
  # max incidence in the data
  
  max(y["cases", ], na.rm = TRUE)
}
probe.cumsum <- function(y){
  # total number of cases in the data
  
  cases <- y["cases", ]
  max(cumsum(cases))
}

plist <- list(
  probe.zeroes,
  probe.max,
  probe.cumsum,
  probe.acf("cases", lags = 1,transform = sqrt, type = "correlation")
)

stew(
  file = "probe_match_results.rda", {
    pm <- probe.match(
      measles_pomp,
      probes = plist,
      est = param_names,
      nsim = 500,
      transform = TRUE,
      start = coef(measles_pomp),
      method = "Nelder-Mead",
      maxit = 10000
    )
  },
  seed = 290868325,
  kind = "L'Ecuyer"
)

