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
  double inv_lambda = exp(-beta_r * (I + psi));
  
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

# Combine with births data for covariate table
covar_data <- biweek_data %>%
  dplyr::select(biweek, births) %>%
  arrange(biweek) %>%
  bind_cols(bspline_basis)


# Set up pomp container ---------------------------------------------------

# Set some realistic parameter values for testing via simulation
params <- c(
  beta_mu = 0.0000006,
  b1 = 3,
  b2 = 0,
  b3 = 1.5,
  b4 = 6,
  b5 = 5,
  b6 = 3,
  psi = 0.1,
  rho = 0.5,
  tau = 0.0000001,
  sigma_env = 0.01,
  S_0 = 50000, 
  I_0 = 200
)

# Generate pomp object
measles_pomp <- pomp(
  data = obs_data,
  time = "biweek",
  covar = covar_data,
  tcovar = "biweek",
  t0 = 1,
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/26),
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

simulate(
  measles_pomp, 
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())


saveRDS(measles_pomp, file = "measles-pomp.RDS")
