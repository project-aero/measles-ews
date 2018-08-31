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


# Define deterministic skeleton -------------------------------------------

measles_skeleton <- Csnippet(
  "
  // Define the variables
  int nrate = 6;        // number of rates
  double rate[nrate];		// transition rates
  double term[nrate];		// terms in the equations
  double lambda;        // force of infection
  double beta;          // static transmission rate, since deterministic
  
  // Calculate force of infection
  beta = beta_mu * (1 + dot_product(K, &xi1, &b1));
  lambda = (iota + beta*I) / pop;

  // Compute the transition rates
  rate[0] = births;	 // birth into susceptible class
  rate[1] = lambda;  // force of infection
  rate[2] = 0;       // death from susceptible class
  rate[3] = gamma;   // recovery
  rate[4] = 0;       // death from infectious class
  rate[5] = 0;       // death from recovered class

  // Compute the state transitions
  term[0] = rate[0];
  term[1] = rate[1]*S;
  term[2] = rate[2]*S;
  term[3] = rate[3]*I;
  term[4] = rate[4]*I;
  term[5] = rate[5]*R;

  // Balance the equations
  DS = term[0]-term[1]-term[2];
  DI = term[1]-term[3]-term[4];
  DR = term[3]-term[5];
  Dcases = term[1];  // accumulate the new S->I transitions
  DW = 0;            // no noise, so no noise accumulation
  "
)


# Define stochastic process (SDEs) ----------------------------------------

measles_process <- Csnippet(
  "
  // Define the variables
  int nrate = 6;        // number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];	// transition numbers
  double lambda;        // force of infection
  double beta;          // transmission rate
  double dW;            // white noise

  // Calculate force of infection
  beta = beta_mu * (1 + dot_product(K, &xi1, &b1));
  lambda = (iota + beta*I);

  // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
  dW = rgammawn(beta_sd, dt);

  // Compute the transition rates
  rate[0] = births;		                    // birth into susceptible class
  rate[1] = (iota + beta*I*dW/dt) / pop;  // force of infection
  rate[2] = 0;			                      // death from susceptible class
  rate[3] = gamma;	                      // recovery
  rate[4] = 0;			                      // death from infectious class
  rate[5] = 0; 		                        // death from recovered class

  // Compute the state transitions
  trans[0] = rpois(rate[0]*dt);	 // births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // Balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  cases += trans[1];  // cases are cumulative infections
  if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;
  Re += (beta*dW/dt)/gamma/pop;
  "
)


# Define likelihood function ----------------------------------------------

measles_dmeasure <- Csnippet(
  "
  double mean;
  double f;
  mean = cases*rho;
  f = dnbinom_mu(reports, 1/tau, mean, give_log);

  lik = (give_log) ? log(f) : f;
  "
)


# Define process simulator for observations -------------------------------

measles_rmeasure <- Csnippet(
  "
  reports = rnbinom_mu(1/tau, rho*cases);
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
  Ttau = exp(tau);
  Tiota = exp(iota);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  TR_0 = exp(R_0);
  "
)

to_estimation <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Tgamma = log(gamma);
  Ttau = log(tau);
  Tiota = log(iota);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  TS_0 = log(S_0);
  TI_0 = log(I_0);
  TR_0 = log(R_0);
  "
)

initial_values <- Csnippet(
  "
  S = nearbyint(S_0);
  I = nearbyint(I_0);
  R = nearbyint(R_0);
  cases = 0;
  W = 0;
  Re = 0;
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

#### REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS ####
obs_data$reports[275] <- round(mean(obs_data$reports[c(274,276)]))

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
  beta_mu = 1.1,
  gamma = 0.6,
  beta_sd = 3e-1,
  tau = 1e-1,
  b1 = 5,
  b2 = 3,
  b3 = 6,
  b4 = 5,
  b5 = 2,
  b6 = 0,
  iota = 2,
  rho = 0.3,
  S_0 = 80000, 
  I_0 = 100,
  R_0 = covar_data$pop[1] - 80000 - 100
)

measles_pomp <- pomp(
  data = obs_data,
  times = "time",
  covar = covar_data,
  tcovar = "time",
  t0 = 1,
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/10),
  skeleton = vectorfield(measles_skeleton),
  rmeasure = measles_rmeasure,
  dmeasure = measles_dmeasure,
  initializer = initial_values,
  statenames = c("S", "I", "R", "cases", "W", "Re"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;",
  zeronames = c("cases", "W", "Re")
)

saveRDS(object = measles_pomp, file = "measles-pomp-object.RDS")


# Extra plotting code for testing -----------------------------------------

# traj <- trajectory(measles_pomp, times=seq(1,52*11,by=1))
# plot(traj[2,1,], type = "l")
# 
simulate(
  measles_pomp,
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = reports, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank()) +
  geom_hline(aes(yintercept = 0))
# 
# logLik(pfilter(measles_pomp, Np = 100))