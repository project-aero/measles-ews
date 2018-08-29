# Example of genetic algorithm for nonlinear population growth

library(tidyverse)
library(GA)


# Define population growth function ---------------------------------------

sim_growth <- function(theta){
  r <- theta[1]
  N0 <- theta[2]
  K <- theta[3]
  sigma <- theta[4]
  ngens <- theta[5]
  
  N <- numeric(ngens)
  N[1] <- N0
  
  for(t in 2:ngens){
    rnow <- r + rnorm(1,0,sigma)
    delta <- rnow*N[t-1]*(1 - (N[t-1]/K)) 
    N[t] <- N[t-1] + delta
  }
  return(N)
}

theta <- c(r = 0.1, N0 = 100, K = 1000, sigma = 0.25, ngens = 50)
outsim <- sim_growth(theta)
plot(outsim)

# Define “fitness” function -----------------------------------------------

fitnessLL <- function(y, theta){
  -sum((y - sim_growth(theta))^2)
}


# Run genetic learning algorithm ------------------------------------------

ga_out <- ga(
  type = "real-valued",
  fitness = fitnessLL,
  y = outsim,
  lower = c(0.01, 10, 100, 0, length(outsim)),
  upper = c(2, 1000, 100000, 1, length(outsim)),
  popSize = 500,
  crossover  = gareal_blxCrossover,
  maxiter = 5000,
  run = 200,
  names = c("r", "N0", "K", "sigma", "ngens")
)

summary(ga_out)
plot(outsim, cex = 0)
for(i in 1:1000){
  lines(sim_growth(theta = ga_out@solution), col = "red", lwd = 0.2)
}
points(outsim, pch = 19)
lines(outsim)

