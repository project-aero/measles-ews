# define-measles-pomp.R
#  R script to define the measles pomp object.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(pomp)

# Data --------------------------------------------------------------------

biweek_data <- readRDS("niger_cleaned_biweekly_measles.RDS")


# Process model -----------------------------------------------------------

measles_process <- Csnippet(
  "
  double dN[2];
  
  double beta_r = beta_mu * (1 + dot_product(K, &xi1, &b1));
  double var_epsilon = sigma_env;
  beta_r *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon, var_epsilon);
  double inv_lambda = exp(-beta_r * (I/N + psi));
  
  if(inv_lambda > 1){
    inv_lambda = 1;
  }
  double delta = rbinom(nearbyint(S), inv_lambda);
  dN[0] = delta + births;
  dN[1] = S - delta;
  S = nearbyint(dN[0]);
  I = nearbyint(dN[1]);
  beta_t = beta_r;
  "
)


# Measurement model -------------------------------------------------------

rmeas <- Csnippet(
  "
  cases = rnbinom_mu(1/tau, rho*I);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
  "
)


# Likelihood model --------------------------------------------------------

dmeas <- Csnippet(
  "
  double mean_cases = nearbyint(rho*I);
  lik = dnbinom_mu(cases, 1/tau, mean_cases, give_log);
  "
)


# Parameter mapping -------------------------------------------------------

init <- Csnippet(
  "
  I = nearbyint(I_0);
  S = nearbyint(S_0);
  beta_t = 0;
  "
)

to_est <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Trho = logit(rho);
  Ttau = log(tau);
  Tsigma_env = log(sigma_env);
  TS_0 = log(S_0);
  TI_0 = log(I_0);
  Tpsi = log(psi);
  "
)

from_est <- Csnippet(
  "
  Tbeta_mu = exp(beta_mu);
  Trho = expit(rho);
  Ttau = exp(tau);
  Tsigma_env = exp(sigma_env);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  Tpsi = exp(psi);
  "
)


# Define data and covariate tables ----------------------------------------

# Fetch observations
obs_data <- biweek_data %>%
  dplyr::select(biweek, cases) %>%
  mutate(cases = round(cases, 0))

# Generate basis functions for seasonality
bspline_basis <- periodic.bspline.basis(
  obs_data$biweek,
  nbasis = 6,
  degree = 3,
  period = 26,
  names = "xi%d"
) %>%
  as_tibble()

# Births smooth
births_oneweek <- biweek_data %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    births = unique(births),
    biweek = min(biweek)
  )
  
births_interp <- predict(
  smooth.spline(
    x = births_oneweek$biweek, 
    y = births_oneweek$births
  ),
  x = biweek_data$biweek)$y %>%
  as_tibble() %>%
  rename(births = value) %>%
  mutate(
    births = round(births),
    biweek = biweek_data$biweek
  )

# Population smooth
pop_oneweek <- biweek_data %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    population = unique(population),
    biweek = min(biweek)
  )

pop_interp <- predict(
  smooth.spline(
    x = pop_oneweek$biweek, 
    y = pop_oneweek$population
  ),
  x = biweek_data$biweek)$y %>%
  as_tibble() %>%
  rename(population = value) %>%
  mutate(
    population = round(population),
    biweek = biweek_data$biweek
  )
  

# Combine with births data for covariate table
covar_data <- births_interp %>%
  arrange(biweek) %>%
  left_join(pop_interp, by = "biweek") %>%
  rename(N = population) %>%
  bind_cols(bspline_basis)


# Set up pomp container ---------------------------------------------------

# Set some realistic parameter values for testing via simulation
params <- c(
  beta_mu = 0.001,
  b1 = 3,
  b2 = 0,
  b3 = 1.5,
  b4 = 6,
  b5 = 5,
  b6 = 3,
  psi = 1,
  rho = 0.5,
  tau = 1,
  sigma_env = 0.01,
  S_0 = 30000, 
  I_0 = 100
)

# Generate pomp object
measles_pomp <- pomp(
  data = obs_data,
  time = "biweek",
  covar = covar_data,
  tcovar = "biweek",
  t0 = 1,
  rprocess = discrete.time.sim(measles_process, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  initializer = init,
  statenames = c("S", "I", "beta_t"),
  toEstimationScale = to_est,
  fromEstimationScale = from_est,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;"
)

saveRDS(measles_pomp, file = "measles-pomp.RDS")
