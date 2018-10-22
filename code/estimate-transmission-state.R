# estimate-transmission-state.R
#  This script runs a plain vanilla particle filter at the MLEs but allowing
#  mean transmission to talk a random walk in time. The filtering distribution
#  for the state of mean transmission at each time gives us the estimate and
#  distribution of transmission over the course of the time series.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


DO_CITY <- "Niamey"
pomp_city <- "Niamey (City)"

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(lubridate)
library(ggthemes)


# Load MLEs ---------------------------------------------------------------

mles <- read.csv(paste0("../results/initial-mif-lls-", DO_CITY, ".csv")) %>%
  filter(loglik == max(loglik, na.rm = T)) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se)


# Update procress C snippet with beta random walk -------------------------

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

  // Beta random walk
  beta_t *= rgammawn(0.001, dt)/dt;
  
  // Calculate force of infection
  seas = (1 + exp(dot_product(K, &xi1, &b1)));
  beta = beta_t*seas;
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
  
  cases += dNIR;  // caseas are cumulative reports at end of infectious period (I->R)
  if (beta_sd > 0.0)  W += (dW-dt)/beta_sd;
  RE = (beta_t / gamma) * (S / N);
  RE_seas = (beta / gamma) * (S / N);
  cases_state = rnbinom_mu(1/tau, rho*cases);
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
  cases = nearbyint(N*I_0);
  cases_state = nearbyint(N*I_0);
  W = 0;
  RE = 0;
  RE_seas = 0;
  beta_t = beta_mu;
  "
)


# Make data tables --------------------------------------------------------

do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
measles_data <- readRDS(do_file) %>%
  dplyr::filter(region == pomp_city)

obs_data <- measles_data %>%
  dplyr::select(time, cases) %>%
  dplyr::rename(
    reports = cases
  )

if(pomp_city == "Niamey (City)"){
  # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
  obs_data$reports[275] <- round(mean(obs_data$reports[c(274,276)]))
}

covar_file <- "../data/clean-data/annual-demographic-data-niger-cities-clean.RDS"
covar_data <- readRDS(covar_file) %>%
  dplyr::filter(region == pomp_city) %>%
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
  t0 = min(obs_data$time),
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
  rmeasure = measles_rmeasure,
  dmeasure = measles_dmeasure,
  initializer = initial_values,
  statenames = c("S", "E", "I", "cases", "cases_state", "W", "RE", "RE_seas", "beta_t"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;",
  zeronames = c("cases", "W", "cases_state")
)


# Run particle filter -----------------------------------------------------

measles_filter <- pfilter(object = measles_pomp, Np=50000, save.states = TRUE, pred.mean = TRUE, pred.var = TRUE)
states <- measles_filter@saved.states  # filtering distribution of states
predicted_means <- measles_filter@pred.mean  # prediction distribution of states
predicted_vars <- measles_filter@pred.var  # prediction variance of states


# Functions to extract and summarize filtered states ----------------------

extract_from_list <- function(filter_list, state_number){
  out <- as_tibble(
    lapply(filter_list, `[`, state_number, )
  ) %>%
    mutate(
      particle = 1:n()
    ) %>%
    gather(key = time, value = value, -particle)
  
  return(out)
}

summarise_filtered_state <- function(df, state_name, observations = NA){
  out <- df %>%
    group_by(time) %>%
    summarise(
      med = median(value),
      upper_95 = quantile(value, 0.975),
      lower_95 = quantile(value, 0.025),
      upper_80 = quantile(value, 0.9),
      lower_80 = quantile(value, 0.1)
    ) %>%
    slice(2:n()) %>%
    mutate(
      week = 1:n(),
      date = measles_data$date[2:nrow(measles_data)],
      state = state_name,
      observation = observations
    )
  
  return(out)
}


# Extract prediction means and variances for cases ------------------------

pred_cases <- t(predicted_means)[,"cases_state"]
pred_stdev <- sqrt(t(predicted_vars)[,"cases_state"])

predictive_distribution <- tibble(
  mean_cases = pred_cases,
  sdev_cases = pred_stdev
)

pred_outfile <- paste0("../results/predictive-dist-states-", DO_CITY, ".RDS")
saveRDS(
  object = predictive_distribution, 
  file = pred_outfile
)



# Extract and summarize states --------------------------------------------

state_names <- names(states[[1]][,1])

eff_rep_nonseasonal <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "RE")
) %>%
  summarise_filtered_state(
    state_name = "effective_r_nonseasonal",
    observations = NA
  )

eff_rep_seasonal <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "RE_seas")
) %>%
  summarise_filtered_state(
    state_name = "effective_r_seasonal",
    observations = NA
  )

transmission_rate <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "beta_t")
) %>%
  summarise_filtered_state(
    state_name = "transmission_rate",
    observations = NA
  )

susceptibles <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "S")
) %>%
  summarise_filtered_state(
    state_name = "susceptibles",
    observations = NA
  )

exposed <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "E")
) %>%
  summarise_filtered_state(
    state_name = "exposed",
    observations = NA
  )

infected <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "I")
) %>%
  summarise_filtered_state(
    state_name = "infected",
    observations = NA
  )

cases <- extract_from_list(
  filter_list =  states, 
  state_number = which(state_names == "cases_state")
) %>%
  summarise_filtered_state(
    state_name = "cases",
    observations = measles_data$cases[2:nrow(measles_data)]
  )


# Combine state dfs and save ----------------------------------------------

all_states <- eff_rep_nonseasonal %>%
  bind_rows(eff_rep_seasonal) %>%
  bind_rows(transmission_rate) %>%
  bind_rows(susceptibles) %>%
  bind_rows(exposed) %>%
  bind_rows(infected) %>%
  bind_rows(cases) %>%
  ungroup() %>%
  group_by(state) %>%
  nest()

outfile <- paste0("../results/filtered-states-", DO_CITY, ".RDS")
saveRDS(
  object = all_states, 
  file = outfile
)

# ggplot(all_states %>% unnest(), aes(x = date)) +
#   geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = ptol_pal()(2)[1], alpha = 0.3) +
#   geom_line(aes(y = med), color = ptol_pal()(2)[1]) +
#   geom_point(aes(y = observation), color = ptol_pal()(2)[2], size = 0.3) +
#   facet_wrap(~state, scales = "free_y") +
#   labs(y = "Filtered median (+/- 95% CI)", x = "Date")


# Calculate log-likelihood of “new” model and save ------------------------
ll <- logmeanexp(
  replicate(
    n = 10,
    logLik(
      object = pfilter(
        object = measles_pomp, 
        Np = 10000
      )
    )
  ), 
  se=TRUE
)

ll_tbl <- ll %>%
  as_tibble() %>%
  mutate(
    var = c("loglik", "loglik_se")
  )

ll_file <- paste0("../results/time-vary-beta-loglik-", DO_CITY, ".csv")
write.csv(x = ll_tbl, file = ll_file)
