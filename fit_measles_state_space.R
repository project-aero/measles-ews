# fit_measles_state_space.R:
#  R script to fit a Bayesian state space model to the measles case data from
#  Niamey, Niger from 1995-2005. This script fits the specified model using
#  MCMC via JAGS. The posterior distributions are summarized and saved. 
#  NOTE: Full MCMC output is NOT saved currently, as this is a development
#  script.
#
# Authors:
#  Andrew Tredennick (atredenn@gmail.com)
#  Pejman Rohani
#  Andrew Park


# Load libraries ----------------------------------------------------------

library(tidyverse)  # data wrangling
library(lubridate)  # time and date functions
library(rjags)      # MCMC sampling algorithms for Bayesian models
library(ggmcmc)     # quick conversions of MCMC output to data frames
library(parallel)   # functions for parallel computing on multiple cores


# Read in data ------------------------------------------------------------

file_name <- "../niger_measles/niger_regional_1995_2005.csv"
niger_measles_raw <- read_csv(file_name, col_types = cols())

num_regions <- nrow(niger_measles_raw)
num_weeks <- ncol(niger_measles_raw) - 1  # subtract 1 from ncol() because first column are regions
weeks_per_year <- num_weeks/11
weeks <- rep(1:52, times = 11)

# Create a vector of years for all num_weeks
years <- rep(1995:2005, each = weeks_per_year)

# Function for calculating start of week based on week number and year
calculate_start_of_week = function(week, year) {
  date <- ymd(paste(year, 1, 1, sep="-"))
  week(date) = week
  return(date)
}

# Clean up the data frame
measles_data <- niger_measles_raw %>%
  gather(key = week, value = cases, -X1) %>%
  mutate(
    week_num = rep(1:num_weeks, each = num_regions),
    year = rep(years, each  = num_regions),
    week = rep(weeks, each = num_regions),
    date = calculate_start_of_week(week, year)
  ) %>%
  dplyr::rename(region = X1) %>%
  filter(region == "Niamey (City)")


# Specify the JAGS model --------------------------------------------------

measles_model <- "
model{
  # Likelihood/Data Model
  for(i in 1:nobs){
    p[i] <- eta/(eta + rho*I[i])
    y[i] ~ dnegbin(p[i], eta)
  }
  eta ~ dunif(0, 50)
  
  # Process model
  for(t in 1:(ntimes-1)){
    # lambda[t] = exp(-beta[t] * (I[t] + psi))
    lambda[t] = (1 + (beta[t] * (I[t] + psi)) / kappa)^(-kappa)
    Delta[t] ~ dbin(lambda[t], S[t])
    I[t+1] = max(0.000001, S[t] - Delta[t])
    S[t+1] = b[t] + Delta[t]
  }
  psi ~ dunif(0, 100)
  kappa ~ dunif(0, 100)
  
  # Parameter model
  for(t in 1:(ntimes-1)){
    beta[t] = exp(gamma[t]) * (1 + upsilon * sin( (2*pi*t)/26 + phi ) )
  }

  for(t in 2:ntimes){
    epsilon[t] ~ dnorm(0, tau_gamma)
    gamma[t] = gamma0 + t * log(1+r) + epsilon[t]
  }
  
  epsilon[1] ~ dnorm(0, tau_gamma)
  gamma[1] = gamma0 + 1 * log(1+r) + epsilon[1]
  gamma0 ~ dunif(-12, -9)  # semi-informed prior based on Ferrari et al. 2008 model results (sd*2)
  tau_gamma = pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 5)
  r ~ dunif(-0.2, 0.2)
  
  upsilon ~ dunif(0, 1)
  phi ~ dunif(0, 26)
  rho ~ dnorm(0.48, 1000) T(0, 1) # informed prior based on Ferrari et al. 2008 model results
  
  # Initial conditions
  S0 ~ dunif(5000, 70000)  # informed prior based on Ferrari et al. 2008 model results
  S[1] <- trunc(S0)
  I0 ~ dpois(initI/rho)
  I[1] = I0
  
  # Derived quantities
  for(i in 1:nobs){
    Iobs[i] = I[i]*rho 
    Rnaught[i] = exp(gamma[i])*N[i]
  }
}
"


# Format data for model fitting -------------------------------------------

population <- read_csv("../niger_measles/district_pops.csv") %>%
  gather(key = year, value = population, -X1) %>%
  rename(district = X1) %>%
  filter(district == "Niamey I")

birth_rates <- read_csv("../niger_measles/niger_crude_birth_rates.csv") %>%
  mutate(
    date = mdy(date),  # lubridate prefixes any 2digit year 00-68 with 20, not a problem for us though
    year = as.character(year(date)),
    rate_per_person = births_per_thousand/1000
  ) %>%
  select(year, rate_per_person)

births <- population %>%
  left_join(birth_rates, by = "year") %>%
  mutate(
    births_per_year = population * rate_per_person,
    births_per_week = births_per_year / 52
  )

# Aggregate weekly data to biweeks
biweek_data <- measles_data %>%
  mutate(biweek = rep(1:(n()/2), each = 2)) %>%
  group_by(biweek) %>%
  summarise(
    cases = sum(cases),
    date = min(date)
  ) %>%
  mutate(
    births = round(rep(births$births_per_week*2, each = 26)),
    population = round(rep(population$population, each = 26))
  ) %>%
  filter(year(date) < 2001)


# Gather data
y <- biweek_data$cases
initI <- biweek_data$cases[1]
b <- biweek_data$births
ntimes <- nrow(biweek_data)
nobs <- nrow(biweek_data)


# Set up JAGS variables ---------------------------------------------------

# List up data for JAGS
jags_data <- list(
  y = y,
  initI = initI,
  b = b,
  ntimes = ntimes,
  nobs = nobs,
  pi = pi,
  N = biweek_data$population
)

# Create initial value function for parallel MCMC
generate_initial_values <- function(){
  init_list <- list(
    Iobs = rpois(nobs, 100),  # observed cases state vector
    I = rpois(nobs, 100/0.5),  # latent infected class state vector
    S = rpois(nobs, 30000),  # latent susceptible class state vector
    gamma = runif(nobs, log(1e-07), log(1e-04)),  # time-varying transmission rate vector
    beta = runif(nobs, 1e-07, 1e-04),  # time-varying seasonal transmission rate vector
    rho = rnorm(1, 0.16, 0.00001),  # reporting fraction, centered on ~0.166
    escape_prob = runif((ntimes-1), 0.9, 0.99),  # time-varying escape-from-infection probability vector
    theta = runif(1, 10, 14),  # phase of sin wave seasonality
    r = runif(1, 0.001, 5),  # dispersion of negative binomial observation process
    sigma_gamma = runif(1, 0.1, 0.6),  # std. dev. of noise on transmission rate
    epsilon = runif(1, 0.4, 0.8),  # amplitude of sin wave seasonality
    S0 = rpois(1, 30000),  # initial condition for susceptible class
    I0 = rpois(1, initI),  # initial condition for infected class
    gamma0 = runif(1, -13, -9),  # initial condition for transmission rate
    m = rpois(1, 10),  # susceptible immigration
    rg = runif(1, -0.006, -0.0009)  # exponential growth/decay rate of transmission rate through time
  )
}


# Run the MCMC ------------------------------------------------------------

# Set MCMC parameters
n_adapt <- 50000  # iterations for adaptation phase
n_update <- 200000  # iterations for burn-in to stationary distributions
n_sample <- 50000  # iterations for sampling from posterior

# Set up cluster
if(detectCores() < 4){  # make sure there are enough cores
  stop("Too few cores on this machine for parallel MCMC.")
}

cl <- makeCluster(3)

# Run JAGS in parallel
clusterExport(
  cl,
  c("jags_data", "generate_initial_values", "n_adapt", "n_update", "n_sample", "measles_model")
)

mcmc_results <- clusterEvalQ(
  cl, {
    library(rjags)  # reload rjags for each core
    set.seed(1)
    
    # Initialize model in JAGS (adaptation phase)
    model <- jags.model(
      textConnection(measles_model),
      data = jags_data,
      n.chains = 1,
      n.adapt = n_adapt
    )
    
    # Update chain with long burn-in
    update(model, n.iter = n_update)
    
    # Collect samples from posterior
    vars_to_collect <- c(
      "Iobs", "I", "S", "Rnaught", "gamma", "beta", "rho", "lambda", 
      "psi", "upsilon", "theta", "eta", "r", "sigma_gamma", "S0", "I0", 
      "gamma0", "kappa"
    )
    
    mcmc_core <- coda.samples(
      model, 
      variable.names = vars_to_collect, 
      n.iter = n_sample, 
      n.thin = 10
    )
    
    return(as.mcmc(mcmc_core))
  }  # end cluster-specific calls
)  # end cluster definition

mcmc_all <- mcmc.list(mcmc_results)  # grab all chains
stopCluster(cl)  # stop the cluster


# Summarize and save output -----------------------------------------------

ggs(mcmc_all) %>%
  group_by(Parameter) %>%
  summarise(
    median_value = median(value),
    upper_95 = quantile(value, 0.975),
    lower_95 = quantile(value, 0.025),
    upper_50 = quantile(value, 0.75),
    lower_50 = quantile(value, 0.25)
  ) %>%
  write_csv(path = "./mcmc_results.csv")
