# define-continuous-measles-pomp.R:
#  Script to define a 'pomp' model object for the dynamics of measles in
#  Niger. The model is a continuous-time SIR model with a time-varying
#  transmission rate. The model is specified as a set of nonlinear stochastic
#  difference equations (SDEs).
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(pomp)


# Define stochastic process (SDEs) ----------------------------------------

measles_process <- Csnippet(
  "
  // Define the variables
  int nrate = 3;          // number of rates
  double rate[nrate];	  	// transition rates
  double trans[nrate];   	// transition numbers
  double lambda;          // force of infection
  double beta;            // transmission rate
  double eta = 365/8;     // infectious rate (8 days latent)
  double gamma = 365/5;   // recovery rate (5 days infectious)
  double dW;              // white noise
  double seas;            // seasonality term
  double dN0S, dN0I, dNSE, dNEI, dNIR;  // transitions

  // Calculate force of infection
  seas = (1 + exp(dot_product(K, &xi1, &b1)));
  beta = beta_mu*seas;
  lambda = beta*I/N;

  // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
  dW = rgammawn(beta_sd, dt);

  // Compute the transition rates
  rate[0] = lambda*dW/dt; // force of infection
  rate[1] = eta;          // infectious rate from latent
  rate[2] = gamma;	      // recovery from infectious

  // Compute the state transitions
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, E, &rate[1], dt, &trans[1]);
  reulermultinom(1, I, &rate[2], dt, &trans[2]);

  // Transitions
  dN0S = rpois(0.3 * mu * N * dt);
  dN0I = rpois(iota * dt);
  dNSE = trans[0];
  dNEI = trans[1];
  dNIR = trans[2];

  // Balance the equations
  S += dN0S - dNSE;
  E +=        dNSE - dNEI;
  I += dN0I        + dNEI - dNIR;

  cases += dNIR;  // cases are cumulative reports at end of infectious period (I->R)
  if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;
  RE = (beta_mu / gamma) * (S / N);
  RE_seas = (beta / gamma) * (S / N);
  "
)


# Define likelihood function ----------------------------------------------

measles_dmeasure <- Csnippet(
  "
  double mean;
  double f;
  mean = cases*rho;
  if (ISNA(reports)) {  // for catching missing observations
    lik = (give_log) ? 0 : 1;
  } else {
    f = dnbinom_mu(reports, 1/tau, mean, give_log);  // negative binomial likelihood
    // f = dpois(reports, mean, give_log);
  }

  lik = (give_log) ? log(f) : f;
  "
)


# Define process simulator for observations -------------------------------

measles_rmeasure <- Csnippet(
  "
  reports = rnbinom_mu(1/tau, rho*cases);  // negative binomial measurement process
  // reports = rpois(rho*cases);

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
  TE_0 = expit(E_0);
  TI_0 = expit(I_0);
  Ttau = exp(tau);
  "
)

to_estimation <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Tiota = log(iota);
  Trho = logit(rho);
  Tbeta_sd = log(beta_sd);
  TS_0 = logit(S_0);
  TE_0 = logit(E_0);
  TI_0 = logit(I_0);
  Ttau = log(tau);
  "
)

initial_values <- Csnippet(
  "
  S = nearbyint(N*S_0);
  E = nearbyint(N*E_0);
  I = nearbyint(N*I_0);
  cases = 0.5*N*I_0;
  W = 0;
  RE = 0;
  RE_seas = 0;
  "
)


# Make data tables --------------------------------------------------------

do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"

all_cities <- readRDS(do_file) %>%
  pull(region) %>%
  unique()

for(do_city in all_cities){
  measles_data <- readRDS(do_file) %>%
    dplyr::filter(region == do_city)
  
  obs_data <- measles_data %>%
    dplyr::select(time, cases) %>%
    dplyr::rename(
      reports = cases
    )
  
  if(do_city == "Niamey (City)"){
    # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
    obs_data$reports[276] <- round(mean(obs_data$reports[c(275,277)]))
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
  
  params <- c(
    beta_mu = 500,
    beta_sd = 0.001,
    b1 = 3,
    b2 = 3,
    b3 = 6,
    b4 = 5,
    b5 = 2,
    b6 = 0,
    iota = 2,
    rho = 0.5,
    S_0 = 0.03, 
    E_0 = 0.00032*0.5,
    I_0 = 0.00032*0.5,
    tau = 0.001
  )
  
  measles_pomp <- pomp(
    data = obs_data,
    times = "time",
    covar = covar_data,
    tcovar = "time",
    t0 = min(obs_data$time),
    rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
    rmeasure = measles_rmeasure,
    dmeasure = measles_dmeasure,
    initializer = initial_values,
    statenames = c("S", "E", "I", "cases", "W", "RE", "RE_seas"),
    toEstimationScale = to_estimation,
    fromEstimationScale = from_estimation,
    paramnames = names(params),
    params = params,
    globals = "int K = 6;",
    zeronames = c("cases", "W")
  )

  city_abb <- substr(do_city, start = 1, stop = 6)
  outfile <- paste0("measles-pomp-object-", city_abb, ".RDS")
  saveRDS(object = measles_pomp, file = outfile)
}



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
# logLik(pfilter(measles_pomp, Np = 2000))
