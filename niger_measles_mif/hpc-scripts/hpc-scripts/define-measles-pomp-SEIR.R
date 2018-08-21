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

model_data <- readRDS("weekly_data_cleaned.RDS")


# Process model -----------------------------------------------------------

measles_process <- Csnippet(
  "
  double dN[3];
  
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
  dN[2] = E;
  S = nearbyint(dN[0]);
  E = nearbyint(dN[1]);
  I = nearbyint(dN[2]);
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
  E = nearbyint(E_0);
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
  TE_0 = log(S_0);
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
  TE_0 = exp(E_0);
  TI_0 = exp(I_0);
  Tpsi = exp(psi);
  "
)


# Define data and covariate tables ----------------------------------------

# Fetch observations
obs_data <- model_data %>%
  dplyr::select(week_num, cases) %>%
  mutate(cases = round(cases, 0))

# Generate basis functions for seasonality
bspline_basis <- periodic.bspline.basis(
  obs_data$week_num,
  nbasis = 6,
  degree = 3,
  period = 52,
  names = "xi%d"
) %>%
  as_tibble()

# Births smooth
births_oneweek <- model_data %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    births = unique(births),
    week = min(week_num)
  )
  
births_interp <- predict(
  smooth.spline(
    x = births_oneweek$week, 
    y = births_oneweek$births
  ),
  x = model_data$week_num)$y %>%
  as_tibble() %>%
  rename(births = value) %>%
  mutate(
    births = round(births),
    week = model_data$week_num
  )

# Population smooth
pop_oneweek <- model_data %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    population = unique(population),
    week = min(week_num)
  )

pop_interp <- predict(
  smooth.spline(
    x = pop_oneweek$week, 
    y = pop_oneweek$population
  ),
  x = model_data$week_num)$y %>%
  as_tibble() %>%
  rename(population = value) %>%
  mutate(
    population = round(population),
    week = model_data$week_num
  )
  

# Combine with births data for covariate table
covar_data <- births_interp %>%
  arrange(week) %>%
  left_join(pop_interp, by = "week") %>%
  rename(N = population) %>%
  bind_cols(bspline_basis)


# Set up pomp container ---------------------------------------------------

# Set some realistic parameter values for testing via simulation
params <- c(
  beta_mu = 0.01,
  b1 = 3,
  b2 = 0,
  b3 = 1.5,
  b4 = 6,
  b5 = 5,
  b6 = 3,
  psi = 1,
  rho = 0.5,
  tau = 0.001,
  sigma_env = 0.1,
  S_0 = 30000, 
  E_0 = 100,
  I_0 = 100
)

# Generate pomp object
measles_pomp <- pomp(
  data = obs_data,
  time = "week_num",
  covar = covar_data,
  tcovar = "week",
  t0 = 1,
  rprocess = discrete.time.sim(measles_process, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  initializer = init,
  statenames = c("S", "I", "E", "beta_t"),
  toEstimationScale = to_est,
  fromEstimationScale = from_est,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;"
)

# simulate(
#   measles_pomp, 
#   nsim = 9,
#   as.data.frame = TRUE,
#   include.data = TRUE) %>%
#   ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
#   geom_line() +
#   scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
#   guides(color = FALSE) +
#   facet_wrap(~sim, ncol = 2) +
#   theme(strip.text=element_blank())

saveRDS(measles_pomp, file = "measles-pomp-seir.RDS")
