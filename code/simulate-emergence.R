# simulate-emergence.R
#  Script to simulate the (re)emergence of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing increasing effective reproduction ratio, via transmission rate.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


DO_CITY <- "Niamey"  # the city to model


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(pomp)
library(spaero)


# Load fitted parameters and pomp model -----------------------------------

mle_file <- paste0("../results/initial-mif-lls-", DO_CITY, ".csv")
mles <- read.csv(mle_file) %>% 
  slice(2:n()) %>%  # ignore first row of storage NAs
  filter(loglik == max(loglik, na.rm = TRUE)) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se, -beta_mu)

pomp_file <- paste0("./measles-pomp-object-", DO_CITY, ".RDS")
fitted_pomp <- readRDS(pomp_file)


# Update pomp object parameters and covars --------------------------------

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
  S = nearbyint(N*0.00001);
  E = nearbyint(N*E_0);
  I = nearbyint(N*I_0);
  cases = rho*N*I_0;
  W = 0;
  RE_seas = 0;
  "
)


# Make data tables --------------------------------------------------------
do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"

all_cities <- readRDS(do_file) %>%
  pull(region) %>%
  unique()

measles_data <- readRDS(do_file) %>%
  dplyr::filter(region == "Niamey (City)")

obs_data <- measles_data %>%
  dplyr::select(time, cases) %>%
  dplyr::rename(
    reports = cases
  )

if(DO_CITY == "Niamey"){
  # REMOVE SUSPICIOUS DATA POINT AND REPLACE WITH MEAN OF NEIGHBORS #
  obs_data$reports[276] <- round(mean(obs_data$reports[c(275,277)]))
}
  
covar_table <- fitted_pomp@covar %>%
  as_tibble() %>%
  mutate(
    time = fitted_pomp@tcovar,
    beta_mu = 370.629
    # beta_mu = seq(1, 200, length.out = n())
  )

simulator_pomp <- pomp(
  data = obs_data,
  times = "time",
  covar = covar_table,
  tcovar = "time",
  t0 = fitted_pomp@t0,
  rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
  rmeasure = measles_rmeasure,
  dmeasure = measles_dmeasure,
  initializer = initial_values,
  statenames = c("S", "E", "I", "cases", "W", "RE_seas"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(mles),
  params = unlist(mles),
  globals = "int K = 6;",
  zeronames = c("cases", "W")
)




# Simulate from the new pomp object ---------------------------------------

model_sims <- simulate(
  simulator_pomp,
  nsim = 500,
  as.data.frame = TRUE,
  include.data = FALSE) %>%
  as_tibble()

re_time_avg <- model_sims %>%
  dplyr::select(time, RE_seas, sim) %>%
  group_by(time) %>%
  summarise(mean_re = mean(RE_seas)) %>%
  mutate(
    year = round(time)
  ) %>%
  group_by(year) %>%
  summarise(time_mean_re = mean(mean_re))

re_one_year <- pull(filter(re_time_avg, time_mean_re > 1), year)

ggplot(re_time_avg, aes(x = year, y = time_mean_re)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_segment(aes(x = re_one_year, xend = re_one_year, y = 0, yend = 1), color = "red") +
  geom_line(data = model_sims, aes(x = time, y = RE_seas, group = sim), alpha = 0.05, color = "grey") +
  geom_line(color = "blue") +
  geom_point(size = 3, color = "blue") +
  scale_x_continuous(breaks = 1995:2006) +
  labs(x = "Date", y = expression(R[E](t))) +
  theme_minimal()


# Calculate EWS over different sections -----------------------------------

ews_times <- model_sims %>%
  filter(round(time) < re_one_year) %>%
  pull(time) %>%
  unique()

ews_time_ids <- tibble(
  time = ews_times,
  half = c(rep("first", length(ews_times)/2), rep("second", length(ews_times)/2))
)

data_for_ews <- model_sims %>%
  dplyr::select(time, reports, sim) %>%
  left_join(ews_time_ids, by = "time")

window_bandwidth <- length(ews_times)/2
ews_out <- {}
for(i in unique(data_for_ews$sim)){
  tmp_data <- data_for_ews %>%
    filter(sim == i)
  
  tmp_first <- spaero::get_stats(
    x = pull(filter(tmp_data, half == "first"), reports), 
    center_trend = "local_constant", 
    center_kernel = "uniform", 
    center_bandwidth = window_bandwidth, 
    stat_trend = "local_constant", 
    stat_kernel = "uniform", 
    stat_bandwidth = window_bandwidth, 
    lag = 1, 
    backward_only = FALSE)$stats
  
  tmp_second <- spaero::get_stats(
    x = pull(filter(tmp_data, half == "second"), reports), 
    center_trend = "local_constant", 
    center_kernel = "uniform", 
    center_bandwidth = window_bandwidth, 
    stat_trend = "local_constant", 
    stat_kernel = "uniform", 
    stat_bandwidth = window_bandwidth, 
    lag = 1, 
    backward_only = FALSE)$stats
  
  tmp_out1 <- as_tibble(tmp_first) %>%
    summarise_all(funs(mean, .args = list(na.rm = TRUE))) %>%
    mutate(
      half = "first",
      sim = i
    )
  
  tmp_out2<- as_tibble(tmp_second) %>%
    summarise_all(funs(mean, .args = list(na.rm = TRUE))) %>%
    mutate(
      half = "second",
      sim = i
    )
  
  tmp_out <- bind_rows(tmp_out1, tmp_out2)
  ews_out <- bind_rows(ews_out, tmp_out)
}


# Plot the results --------------------------------------------------------

ews_long <- ews_out %>%
  gather(key = metric, value = value, -half, -sim)

ggplot(ews_long, aes(x = half, y = value)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free_y")
