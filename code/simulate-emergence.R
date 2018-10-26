# simulate-emergence.R
#  Script to simulate the (re)emergence of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing increasing effective reproduction ratio, via transmission rate.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(pomp)
library(pROC)
library(spaero)


all_sims <- tibble()

for(DO_CITY in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model -----------------------------------
  
  mle_file <- paste0("../results/initial-mif-lls-", DO_CITY, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se, -beta_mu)
  
  mle_beta <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    pull(beta_mu)
  
  pomp_file <- paste0("./measles-pomp-object-", DO_CITY, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  
  # Update pomp object parameters and covars --------------------------------
  
  # Define stochastic process (SDEs) 
  
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
  
  
  # Define likelihood function 
  
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
  
  
  # Define process simulator for observations 
  
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
  
  
  # Define parameter transformation scales 
  
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
    S = nearbyint(N*sfact);
    E = nearbyint(N*E_0);
    I = nearbyint(N*I_0);
    cases = rho*N*I_0;
    W = 0;
    RE_seas = 0;
    "
  )
  
  
  # Make data tables
  do_file <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
  
  all_cities <- readRDS(do_file) %>%
    pull(region) %>%
    unique()
  
  measles_data <- readRDS(do_file) %>%
    dplyr::filter(region == paste0(DO_CITY, " (City)"))
  
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
      beta_mu = mle_beta
      # beta_mu = seq(1, 200, length.out = n())
    )
  
  sfact = 0.00001
  global_str <- paste0("int K = 6; ", "double sfact = ", sfact, ";")
  
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
    globals = global_str,
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
    summarise(time_mean_re = mean(mean_re)) %>%
    mutate(
      diff_one = 1 - time_mean_re,
      city = DO_CITY
    ) %>%
    nest(-city)
  
  tmp_re_sims <- model_sims %>%
    dplyr::select(sim, time, RE_seas, reports) %>%
    mutate(city = DO_CITY) %>%
    nest(-city)
  
  all_sims <- bind_rows(all_sims, tmp_re_sims)
}


# Save the simulations ----------------------------------------------------

saveRDS(object = all_sims, file = "../simulations/emergence-simulations.RDS")





#
# 
# 
# # Calculate AUC -----------------------------------------------------------
# 
# cats <- tibble(
#   half = c("first", "second"),
#   cat = c(0, 1)
# )
# 
# ews_long <- ews_long %>%
#   left_join(cats, by = "half")
# 
# auc_tbl <- {}
# for(do_metric in unique(ews_long$metric)){
#   tmp <- filter(ews_long, metric == do_metric)
#   roc_obj <- roc(tmp$cat, tmp$value)
#   tmp_auc <- auc(roc_obj)
#   plot(roc_obj, main = do_metric)
#   tmp_tbl <- tibble(
#     metric = do_metric,
#     auc = as.numeric(tmp_auc)
#   )
#   
#   auc_tbl <- bind_rows(auc_tbl, tmp_tbl)
# }
# 
# plt_tbl <- auc_tbl %>%
#   filter(metric %in% c("variance", "autocovariance", "autocorrelation", 
#                        "decay_time", "mean"))
#   
# ggplot(plt_tbl, aes(x = metric, y = auc-0.5)) +
#   geom_col() +
#   coord_flip() +
#   scale_y_continuous(limits = c(0,0.5)) +
#   theme_minimal()
