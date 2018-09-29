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
library(ggthemes)


# Load MLEs ---------------------------------------------------------------

mles <- read_csv("../results/initial-mif-lls.csv") %>%
  filter(loglik == max(loglik, na.rm = T)) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se)

mles <- read_csv("~/Desktop/initial-mif-lls.csv") %>%
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
  
  cases += dNIR;  // cases are cumulative reports at end of infectious period (I->R)
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
  f = dnbinom_mu(reports, 1/tau, mean, give_log);  // negative binomial likelihood
  // f = dpois(reports, mean, give_log);
  
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
  beta_t = 475;
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
  statenames = c("S", "E", "I", "cases", "W", "RE", "beta_t"),
  toEstimationScale = to_estimation,
  fromEstimationScale = from_estimation,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;",
  zeronames = c("cases", "W")
)



# Run particle filter -----------------------------------------------------

test <- pfilter(object = measles_pomp, Np=5000, save.states = TRUE)
states <- test@saved.states  # save the states separately


# Extract transmission rate -----------------------------------------------

out <- as_tibble(lapply(states, `[`,7,)) %>%
  mutate(
    particle = 1:n()
  ) %>%
  gather(key = time, value = beta, -particle)

transmission_ts <- out %>%
  group_by(time) %>%
  summarise(
    med = median(beta),
    upper = quantile(beta, 0.975),
    lower = quantile(beta, 0.025),
    upper2 = quantile(beta, 0.9),
    lower2 = quantile(beta, 0.1)
  ) %>%
  slice(2:n()) %>%
  mutate(
    week = 1:n(),
    date = measles_data$date[2:nrow(measles_data)]
  )


# Plot the transmission time series ---------------------------------------

trans_plot <- ggplot(transmission_ts, aes(x = date, y = med)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = ptol_pal()(1)) +
  # geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.4, fill = "coral") +
  geom_hline(aes(yintercept = mles$beta_mu), color = ptol_pal()(2)[2], linetype = 2) +
  geom_line(color = ptol_pal()(1)) +
  theme_minimal() +
  labs(x = "Date", y = expression(paste("Mean trasnmission rate, ",beta," (",yr^-1,")"))) 

ggsave(
  plot = trans_plot,
  filename = "../figures/transmission-ts-posts.pdf", width = 5, height = 4
)


# Calculate correlation with time -----------------------------------------

cor.test(transmission_ts$week, transmission_ts$med, method = "kendall")


# Save the transmission time series ---------------------------------------

write_csv(transmission_ts, path = "../results/transmission-posteriors.csv")


# Extract S, I states -----------------------------------------------------

out_S <- as_tibble(lapply(states, `[`,1,)) %>%
  mutate(
    particle = 1:n()
  ) %>%
  gather(key = time, value = susceptible, -particle)

out_I <- as_tibble(lapply(states, `[`,4,)) %>%
  mutate(
    particle = 1:n()
  ) %>%
  gather(key = time, value = infected, -particle) %>%
  mutate(
    `reported cases` = rpois(n(), infected*as.numeric(mles["rho"]))
    # infected = rnbinom(n(), size = 20, mu = infected*as.numeric(mles["rho"]))
  ) %>%
  dplyr::select(-infected)


out_si <- out_S %>%
  left_join(out_I, by = c("particle", "time")) %>%
  gather(key = state, value = abundance, -particle, -time)

states_filtered <- out_si %>%
  group_by(state, time) %>%
  summarise(
    med_abundance = median(abundance),
    upper_abundance = quantile(abundance, 0.975),
    lower_abundance = quantile(abundance, 0.025)
  ) %>%
  ungroup() %>%
  mutate(
    date = rep(measles_data$date, times = 2),
    observations = c(measles_data$cases, rep(NA, nrow(measles_data)))
  )

ggplot(data = states_filtered, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_abundance, ymax = upper_abundance), alpha = 0.5) +
  geom_point(aes(y = observations), color = "red", size = 0.3) +
  geom_line(aes(y = med_abundance)) +
  labs(x = "Date", y = "Number of persons") +
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  scale_y_sqrt() +
  theme_minimal()


