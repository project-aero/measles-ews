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


# Load libraries ----------------------------------------------------------

library(tidyverse)  # data wrangling
library(lubridate)  # time and date functions
library(rjags)      # MCMC sampling algorithms for Bayesian models
library(ggmcmc)     # quick conversions of MCMC output to data frames


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
    p[i] <- r/(r + rho*I[i])
    y[i] ~ dnegbin(p[i], r)
  }
  r ~ dunif(0, 50)
  
  # Process model
  for(t in 1:(ntimes-1)){
    escape_prob[t] = exp(-beta[t] * (I[t] + m))
    escapees[t] ~ dbin(escape_prob[t], S[t])
    I[t+1] = max(0.000001, S[t] - escapees[t])
    S[t+1] = b[t] + escapees[t]
  }
  m ~ dunif(0, 100)
  
  # Parameter model
  for(t in 1:(ntimes-1)){
    beta[t] = exp(gamma[t]) * (1 + epsilon * sin( (2*pi*t)/26 + theta ) )
  }

  for(t in 2:ntimes){
    noise[t] ~ dnorm(0, tau_gamma)
    gamma[t] = psi*gamma[t-1] + psi2*noise[t]
  }
  
  noise[1] ~ dnorm(0, tau_gamma)
  gamma[1] = psi*gamma0 + psi2*noise[1]
  gamma0 ~ dunif(-12, -9)  # prior ranging from Re = 1 to Re = 25
  tau_gamma = pow(sigma_gamma, -2)
  sigma_gamma ~ dunif(0, 2)
  psi ~ dunif(0.8, 1.3)  # slightly weird range b/c AR(1) is on log scale now
  psi2 ~ dunif(0,1)
  
  epsilon ~ dunif(0, 1)
  theta ~ dunif(0, 26)
  rho ~ dbeta(3.6, 6.3)
  
  # Initial conditions
  S0 ~ dunif(50000, 400000)
  S02 ~ dpois(S0)
  S[1] <- S02
  I0 ~ dpois(initI/rho)
  I[1] = I0
  
  # Derived quantities
  for(i in 1:nobs){
    Iobs[i] = I[i]*rho 
    Rnaught[i] = exp(gamma[i])*N[i]
  }
}
"


# Fit the model with JAGS -------------------------------------------------

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
  filter(year(date) < 2000)


# Gather data
y <- biweek_data$cases
initI <- biweek_data$cases[1]
b <- biweek_data$births
ntimes <- nrow(biweek_data)
nobs <- nrow(biweek_data)

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

# Initialize model in JAGS
model <- jags.model(
  textConnection(measles_model),
  data = jags_data,
  n.chains = 1,
  n.adapt = 75000
)

vars_to_collect <- c(
  "Iobs", "I", "S", "Rnaught", "gamma", "beta", "rho", "escape_prob", 
  "psi", "psi2", "theta", "r", "sigma_gamma","epsilon", "S02", "I0", 
  "gamma0", "m"
)

mcmc_results <- coda.samples(
  model, 
  variable.names = vars_to_collect, 
  n.iter = 200000, 
  n.thin = 10
)


# Summarize and save output -----------------------------------------------

ggs(mcmc_results) %>%
  group_by(Parameter) %>%
  summarise(
    median_value = median(value),
    upper_95 = quantile(value, 0.975),
    lower_95 = quantile(value, 0.025),
    upper_50 = quantile(value, 0.75),
    lower_50 = quantile(value, 0.25)
  ) %>%
  write_csv(path = "./mcmc_results.csv")
