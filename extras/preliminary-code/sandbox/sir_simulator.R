# stoch_sir.R:
#  Simulation version of a stochastic SIR model for measles dynamics in
#  Niamey, Niger. The goal here is just to hone in on model functional
#  form a what different parameter values do to dynamics.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Preliminaries -----------------------------------------------------------

rm(list = ls(all.names = ))
library(tidyverse)


# Define model function ---------------------------------------------------

sim_sir <- function(ntimes, mean_beta, psi, theta, epsilon, rho, I0, S0, births){
  I <- numeric(ntimes)
  S <- numeric(ntimes)
  gamma <- numeric(ntimes)
  I[1] <- I0
  S[1] <- S0
  gamma[1] <- mean_beta
  
  for(t in 2:ntimes){
    gamma[t] <- max(0.0000000001, psi*gamma[t-1] + (1 - psi) * rnorm(1, 0, 0.001))
    beta <- gamma[t] * (1 + epsilon * sin( (2*pi*t)/26 + theta ) )
    escape_prob = exp(-beta * (I[t-1] + 1))
    escapees = rbinom(1, S[t-1], escape_prob)
    I[t] = S[t-1] - escapees
    S[t] = births + escapees
  }
  return(list(I, S, gamma))
}
# 
# runout <- sim_sir(ntimes = 100, mean_beta = exp(-9), psi = 0.9, theta = 20, epsilon = 1, rho = 0.5, I0 = 40, S0 = 50000, births = 1000)
# 
# par(mfrow = c(3, 1))
# plot(runout[[1]], type = "l", ylab = "I")
# plot(runout[[2]], type = "l", ylab = "S")
# plot(runout[[3]], type = "l", ylab = expression(beta))




# Simulate with estimated parameters --------------------------------------

params <- read_csv("../mcmc_results.csv")

beta_vec <- params %>%
  filter(grepl("beta\\[", Parameter)) %>%
  pull(median_value)

susc_vec <- params %>%
  filter(grepl("S\\[", Parameter)) %>%
  pull(median_value)

infect_vec <- params %>%
  filter(grepl("I\\[", Parameter)) %>%
  pull(median_value)

sim_sir <- function(ntimes, beta, I0, S0, births){
  I <- numeric(ntimes)
  S <- numeric(ntimes)
  I[1] <- I0
  S[1] <- S0
  
  for(t in 2:ntimes){
    escape_prob = exp(-beta[t] * (I[t-1] + 1))
    escapees = rbinom(1, S[t-1], escape_prob)
    I[t] = S[t-1] - escapees
    S[t] = births + escapees
  }
  return(list(I, S))
}

testout <- sim_sir(ntimes = length(beta_vec), beta = beta_vec, I0 = infect_vec[1], S0 = susc_vec[1], births = 1000)

par(mfrow = c(2,1))
plot(testout[[1]], type = "l")
plot(testout[[2]], type = "l")
