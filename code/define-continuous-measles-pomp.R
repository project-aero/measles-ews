# define-continuous-measles-pomp.R:
#  Script to define a 'pomp' model object for the dynamics of measles in
#  Niger. The model is a continuous-time SIR model with a time-varying
#  transmission rate. The model is specified as a set of nonlinear stochastic
#  differential equations (SDEs).
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


##  https://github.com/theresasophia/pomp-astic/blob/master/mle_stst%2B_study.R


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(pomp)


# Define stochastic process (SDEs) ----------------------------------------

measles_process <- Csnippet(
  "
  // Define the variables
  int nrate = 2;        // number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];	// transition numbers
  double lambda;        // force of infection
  double beta;          // transmission rate
  double dW;            // white noise
  double seas;          // seasonality term
  double dN0S, dN0I, dNSI, dNIR;  // transitions

  // Calculate force of infection
  seas = (1 + exp(dot_product(K, &xi1, &b1)));
  beta = beta_mu*seas;
  lambda = beta*I/pop;

  // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
  dW = rgammawn(beta_sd, dt);

  // Compute the transition rates
  rate[0] = lambda*dW/dt;                 // force of infection
  rate[1] = gamma;	                      // recovery

  // Compute the state transitions
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, I, &rate[1], dt, &trans[1]);

  // Transitions
  dN0S = rpois(births * dt);
  dN0I = rpois(iota * dt);
  dNSI = trans[0];
  dNIR = trans[1];

  // Balance the equations
  S += dN0S - dNSI;
  I += dN0I + dNSI - dNIR;
  R +=               dNIR;

  cases += dNSI;  // cases are cumulative infections (I->R)
  if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;
  RE = (beta * dW/dt) / gamma;
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
  Tgamma = exp(gamma);
  Tiota = exp(iota);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  "
)

to_estimation <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Tgamma = log(gamma);
  Tiota = log(iota);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  TS_0 = log(S_0);
  TI_0 = log(I_0);
  "
)

initial_values <- Csnippet(
  "
  S = nearbyint(S_0);
  I = nearbyint(I_0);
  R = nearbyint(R_0);
  cases = 0;
  W = 0;
  RE = 0;
  "
)


# Make data tables --------------------------------------------------------

do_city <- "Niamey (City)"
do_file <- "../data/clean-data/weekly-measles-and-demog-niger-cities-clean.RDS"
measles_data <- readRDS(do_file) %>%
  dplyr::filter(region == do_city)

obs_data <- measles_data %>%
  dplyr::select(obs_week, cases) %>%
  dplyr::rename(
    time = obs_week,
    reports = cases
  )

if(do_city == "Niamey (City)"){
  # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
  obs_data$reports[275] <- round(mean(obs_data$reports[c(274,276)]))
}

covar_data <- measles_data %>%
  dplyr::select(obs_week, population_smooth, births_per_week_smooth) %>%
  dplyr::rename(
    time = obs_week,
    pop = population_smooth,
    births = births_per_week_smooth
  )

# Generate basis functions for seasonality
bspline_basis <- periodic.bspline.basis(
  obs_data$time,
  nbasis = 6,
  degree = 3,
  period = 52,
  names = "xi%d"
) %>%
  as_tibble()

covar_data <- bind_cols(covar_data, bspline_basis)


# Combine everything into a pomp object -----------------------------------

params <- c(
  beta_mu = 25,
  gamma = 16,
  beta_sd = 0.001,
  b1 = 3,
  b2 = 3,
  b3 = 6,
  b4 = 5,
  b5 = 2,
  b6 = 0,
  iota = 2,
  rho = 0.5,
  S_0 = 20000, 
  I_0 = 100
)

measles_pomp <- pomp(
  data = obs_data,
  times = "time",
  covar = covar_data,
  tcovar = "time",
  t0 = 1,
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
  rmeasure = measles_rmeasure,
  dmeasure = measles_dmeasure,
  initializer = initial_values,
  statenames = c("S", "I", "R", "cases", "W", "RE"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;
             double R_0 = 700000;",
  zeronames = c("cases", "W")
)

saveRDS(object = measles_pomp, file = "measles-pomp-object.RDS")


# Extra plotting code for testing -----------------------------------------

# simulate(
#   measles_pomp,
#   nsim = 9,
#   as.data.frame = TRUE,
#   include.data = TRUE) %>%
#   ggplot(aes(x = time, y = reports, group = sim, color = (sim == "data"))) +
#   geom_line() +
#   scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
#   guides(color = FALSE) +
#   facet_wrap(~sim, ncol = 2) +
#   scale_y_sqrt() +
#   theme(strip.text=element_blank()) +
#   geom_hline(aes(yintercept = 0))
# 
# logLik(pfilter(measles_pomp, Np = 1000))
