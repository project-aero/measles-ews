library(pomp)
library(tidyverse)
mf <- readRDS("~/Desktop/mf-test.RDS")

test_mif <- function(mf, iters_to_test = 75, alpha = 0.05){
  # Tests whether MIF algorithm has converged by fitting a linear trend
  # to the log likelihoods from the final n iterations and testing whether 
  # the trend is significant.
  
  # Args:
  #  mf: A MIF object
  #  iters_to_test: Number of final iterations to use in trend fitting
  #  alpha: Confidence level for significance test
  # Returns:
  #  A logical value indicating whether the trend is significant (TRUE) or
  #  non-significant (FALSE).
  
  lls <- as.numeric(mf@conv.rec[ , 1])  # extract log likelihoods
  lls <- lls[!is.na(lls)]  # ignore final MIF ll, which is always NA
  lls <- tail(lls, iters_to_test)  # grab final `iters_to_test` iterations
  mod <- summary(lm(lls ~ seq_along(lls)))  # fit linear model to trend
  pval <- mod$coefficients[2, 4]  # extract p-value for trend
  return(pval < alpha)  # return test for significant trend
}

run_mif <- function(mf, rw_sd_setup, iters = 100, particles = 30000, cooling_fraction = 1){
  # Runs n iterations of MIF2, starting from a previous MIF object (mf).
  #
  # Args:
  #  mf: A pomp::mif2 object
  #  rw_sd_setup: A function stored as an object that defines the random walk
  #    parameters for the pomp::mif2 algorithm. See ?pomp::mif2() for details.
  #  iters: Number of MIF iterations
  #  particles: Number of particles
  #  cooling_fraction: Rate at which to "cool" the random walk of parameters
  # Returns:
  #  A pomp::mif2 object
  
  out <- pomp::mif2(mf, Nmif = iters, Np = particles, transform = TRUE,
                    cooling.fraction.50 = cooling_fraction,
                    cooling.type = "geometric", rw.sd = rw_sd_setup
  )
  return(out)
}

test_mif(mf)
