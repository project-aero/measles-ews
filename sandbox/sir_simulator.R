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
    gamma[t] <- max(0.0000000001, psi*gamma[t-1] + (1 - psi) * rnorm(1, 0, 0.0002))
    beta <- gamma[t] * (1 + epsilon * sin( (2*pi*t)/26 + theta ) )
    escape_prob = exp(-beta * (I[t-1] + 1))
    escapees = rbinom(1, S[t-1], escape_prob)
    I[t] = S[t-1] - escapees
    S[t] = births + escapees
  }
  return(list(I, S, gamma))
}

runout <- sim_sir(ntimes = 100, mean_beta = 0.0000005, psi = 0.9, theta = 20, epsilon = 1, rho = 0.5, I0 = 40, S0 = 20000, births = 1000)

par(mfrow = c(3, 1))
plot(runout[[1]], type = "l", ylab = "I")
plot(runout[[2]], type = "l", ylab = "S")
plot(runout[[3]], type = "l", ylab = expression(beta))



beta <- 0.00001
I <- 50
exp(-beta*I)
exp(log(beta)+log(I))


