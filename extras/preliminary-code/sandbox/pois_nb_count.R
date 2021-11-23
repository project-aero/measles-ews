# pois_nb_count.R
#  Script to fit a model to a simulated Poisson count process of Malthusian growth.
#  Models are fit assuming Poisson and Negative Binomial data models.
#
# Author: Andrew Tredennick


# Libraries ---------------------------------------------------------------

library(rjags)
library(coda)


# Simulate Malthusian growth ----------------------------------------------

get_malth <- function(r = 0.1, noise = 1, years = 40){
  N <- numeric(years)
  N[1] <- rpois(1, 10)
  for(t in 2:years){
    rnow <- rnorm(1, r, noise)
    # N[t] <- rpois(1,(1+rnow)*N[t-1])
    mu <- (1+rnow)*N[t-1]
    # if(t > 40 & t < 45) mu <- rpois(1, 10)
    N[t] <- rnegbin(1, mu = mu, theta = 10)
  }
  return(N)
}

pop <- get_malth(noise = 0.1)
plot(pop)


# Define JAGS Poisson model ------------------------------------------------

pois_model <- "
model{
  for(i in 1:n){
    y[i] ~ dpois(z[i])
  }

  for(t in 2:n){
    z[t] ~ dnorm(mu[t], tau)
    mu[t] <- (1+r)*z[t-1]
  }

  z[1] ~ dunif(0, 20)
  tau ~ dgamma(0.001, 0.001)
  r ~ dnorm(0, 0.1)
}
"


# Fit JAGS Poisson model --------------------------------------------------

jags_data <- list(y = pop, n = length(pop))

# Initialize model in JAGS
model <- jags.model(
  textConnection(pois_model),
  data = jags_data,
  n.chains = 1,
  n.adapt = 2000
)

update(model, n.iter = 10000)

mcmc_results <- coda.samples(
  model, 
  variable.names = c("r", "z"), 
  n.iter = 10000, 
  n.thin = 1
)

zm_quants <- summary(mcmc_results)$quantile
states <- zm_quants[grep("z",rownames(zm_quants)),]

plot(pop)
lines(states[,3])
lines(states[,1], lty=2)
lines(states[,5], lty=2)


# Define JAGS NB model -----------------------------------------------------

nb_model <- "
model{
  for(i in 1:n){
    p[i] <- eta/(eta+z[i])
    y[i] ~ dnegbin(p[i], eta)
  }

  for(t in 2:n){
    z[t] ~ dnorm(mu[t], tau)
    mu[t] <- (1+r)*z[t-1]
  }
  
  z[1] ~ dunif(0, 20)
  tau ~ dgamma(0.001, 0.001)
  r ~ dnorm(0, 0.1)
  eta ~ dunif(0, 50)
}
"


# Fit JAGS Poisson model --------------------------------------------------

jags_data <- list(y = pop, n = length(pop))

# Initialize model in JAGS
model <- jags.model(
  textConnection(nb_model),
  data = jags_data,
  n.chains = 1,
  n.adapt = 2000
)

update(model, n.iter = 10000)

mcmc_results <- coda.samples(
  model, 
  variable.names = c("r", "z"), 
  n.iter = 10000, 
  n.thin = 1
)

zm_quants <- summary(mcmc_results)$quantile
states <- zm_quants[grep("z",rownames(zm_quants)),]

lines(states[,3], col = "red")
lines(states[,1], lty=2, col = "red")
lines(states[,5], lty=2, col = "red")

