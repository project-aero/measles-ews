# estimate-transmission-state.R
#  This script runs replicate particle filters at the MLEs but allowing
#  mean transmission to talk a random walk in time. The filtering distribution
#  for the state of mean transmission at each time gives us the estimate and
#  distribution of transmission over the course of the time series.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(lubridate)
library(ggridges)


# Load MLEs ---------------------------------------------------------------

mles <- read_csv("../results/initial-mif-lls.csv") %>%
  filter(loglik == max(loglik, na.rm = T)) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se)


# Update procress C snippet with beta random walk -------------------------

measles_process <- Csnippet(
  "
  // Define the variables
  int nrate = 2;        // number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];	// transition numbers
  double lambda;        // force of infection
  double beta;          // transmission rate
  double beta_t;
  double gamma = 365/14;   // recovery rate (14 days)
  double dW;            // white noise
  double seas;          // seasonality term
  double dN0S, dN0I, dNSI, dNIR;  // transitions
  
  // Beta random walk
  beta_t = exp(rnorm(0, 1))*beta_mu;

  // Calculate force of infection
  seas = (1 + exp(dot_product(K, &xi1, &b1)));
  beta = beta_t*seas;
  lambda = beta*I/N;
  
  // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
  dW = rgammawn(beta_sd, dt);
  
  // Compute the transition rates
  rate[0] = lambda*dW/dt;                 // force of infection
  rate[1] = gamma;	                      // recovery
  
  // Compute the state transitions
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, I, &rate[1], dt, &trans[1]);
  
  // Transitions
  dN0S = rpois(0.3 * mu * N * dt);
  dN0I = rpois(iota * dt);
  dNSI = trans[0];
  dNIR = trans[1];
  
  // Balance the equations
  S += dN0S - dNSI;
  I += dN0I + dNSI - dNIR;
  
  cases += dNSI;  // cases are cumulative infections (S->I)
  if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;
  RE = (beta * dW/dt) / gamma;
  beta_state = beta_t;
  "
)

# Define likelihood function ----------------------------------------------

measles_dmeasure <- Csnippet(
  "
  double mean;
  double f;
  mean = cases*rho;
  // f = dnbinom_mu(reports, 1/tau, mean, give_log);  // negative binomial likelihood
  f = dpois(reports, mean, give_log);  // poisson likelihood
  
  lik = (give_log) ? log(f) : f;
  "
)


# Define process simulator for observations -------------------------------

measles_rmeasure <- Csnippet(
  "
  // reports = rnbinom_mu(1/tau, rho*cases);  // negative binomial measurement process
  reports = rpois(rho*cases);  // poisson measurement process
  if (reports > 0.0) {
  reports = nearbyint(reports);
  } else {
  reports = 0.0;
  }
  "
)


# Define parameter transformation scales ----------------------------------

from_estimation <- Csnippet(
  "
  Tbeta_mu = exp(beta_mu);
  Tiota = exp(iota);
  Trho = expit(rho);
  Tbeta_sd = exp(beta_sd);
  TS_0 = expit(S_0);
  TI_0 = expit(I_0);
  "
)

to_estimation <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Tiota = log(iota);
  Trho = logit(rho);
  Tbeta_sd = log(beta_sd);
  TS_0 = logit(S_0);
  TI_0 = logit(I_0);
  "
)

initial_values <- Csnippet(
  "
  S = nearbyint(N*S_0);
  I = nearbyint(N*I_0);
  cases = 0.5*N*I_0;
  W = 0;
  RE = 0; 
  beta_state = 0;
  "
)


# Make data tables --------------------------------------------------------

do_city <- "Niamey (City)"
do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
measles_data <- readRDS(do_file) %>%
  dplyr::filter(region == do_city)

obs_data <- measles_data %>%
  dplyr::select(time, cases) %>%
  dplyr::rename(
    reports = cases
  )

if(do_city == "Niamey (City)"){
  # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
  obs_data$reports[275] <- round(mean(obs_data$reports[c(274,276)]))
}

covar_file <- "../data/clean-data/annual-demographic-data-niger-cities-clean.RDS"
covar_data <- readRDS(covar_file) %>%
  dplyr::filter(region == do_city) %>%
  dplyr::select(time, population_size, birth_per_person_per_year) %>%
  dplyr::rename(
    N = population_size,
    mu = birth_per_person_per_year
  )

# Generate basis functions for seasonality
bspline_basis <- periodic.bspline.basis(
  covar_data$time,
  nbasis = 6,
  degree = 3,
  period = 1,
  names = "xi%d"
) %>%
  as_tibble()

covar_data <- bind_cols(covar_data, bspline_basis)

params <- unlist(mles)

measles_pomp <- pomp(
  data = obs_data,
  times = "time",
  covar = covar_data,
  tcovar = "time",
  t0 = 1995.000,
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
  rmeasure = measles_rmeasure,
  dmeasure = measles_dmeasure,
  initializer = initial_values,
  statenames = c("S", "I", "cases", "W", "RE", "beta_state"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;",
  zeronames = c("cases", "W")
)


test <- pfilter(object = measles_pomp, Np=20000, save.states = TRUE)
states <- test@saved.states

out <- as_tibble(lapply(states, `[`,6,)) %>%
  mutate(
    particle = 1:n()
  ) %>%
  gather(key = time, value = beta, -particle)

transmission_ts <- out %>%
  group_by(time) %>%
  summarise(
    med = median(beta),
    upper = quantile(beta, 0.80),
    lower = quantile(beta, 0.20)
  ) %>%
  slice(2:n()) %>%
  mutate(
    week = 1:n()
  )

ggplot(transmission_ts, aes(x = week, y = med)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  geom_line() +
  theme_minimal() +
  scale_y_log10() +
  labs(x = "Week", y = expression(paste("Trasnmission rate, log(",beta,")")))

cor.test(transmission_ts$week, transmission_ts$med, method = "kendall")
