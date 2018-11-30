# simulate-elimination-grid.R
#  Script to simulate the elimination of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing decreasing susceptible pool due to vaccination.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Define function for vaccination campaign --------------------------------

rho_curve_ramp <- function(t, start = 52*4, speed = 0.0015){
  ifelse(
    t <= start,
    rho <- 0.7,
    # rho <- 0.3 * (1 - exp((t - start) * -0.015 speed)) + 0.7  # exponential
    rho <- min(0.7 + (t - start)*speed, 1)  # linear
  )
  return(rho)
}


# Define function to calculate R0 from seasonal params --------------------

calc_R0 <- function(beta, qis, season, eta = (365/8), 
                    mu = 0.05, nu = 0.05, gamma = (365/5)){
  B <- as.numeric((1 + exp(season %*% qis)) * beta)
  R0 <- (eta*B*mu) / (nu*(eta + nu)*(gamma + nu))
  return(R0)
}


# Make example plot of vaccination curve ----------------------------------

test <- sapply(0:520, FUN = rho_curve_ramp)
# pdf(file = "../figures/vaccination-coverage-example.pdf", width = 6, height = 4.5)
plot(test, type = "l", xlab = "Time (weeks)",
     ylab = "Vaccination coverage", col = "dodgerblue4", lwd = 2)
abline(h = 0.7, lty = 2)
abline(h = 0.95, lty = 3)
text(420, 0.71, labels = "Current vaccination coverage in Niger", cex = 0.7)
text(125, 0.94, labels = "Vaccination coverage for herd immunity", cex = 0.7)
# dev.off()


# Set up grid of vaccination roll out speeds ------------------------------

# speed_grid <- seq(-5e-04, -1e-04, length.out = 6)  # exponential
speed_grid <- seq(0.000015, 0.0001, length.out = 6)  # linear


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)  # functions for parallel computing

if(parallel::detectCores() <= 4){
  registerDoParallel(cores = 4)
} else{
  registerDoParallel()
}

source("make-pomp-simulator-function.R")

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model -----------------------------------
  
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  pomp_file <- paste0("./measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  # Calculte R0 -------------------------------------------------------------
  qis <- mles %>%
    dplyr::select(b1, b2, b3, b4, b5, b6) %>%
    as.numeric()
  
  beta <- mles %>%
    pull(beta_mu)
  
  bases <- as_tibble(fitted_pomp@covar) %>%
    dplyr::select(starts_with("x")) %>%
    dplyr::slice(1:365) %>%
    mutate(
      day = 1:365
    ) %>%
    gather(key = base, value = value, -day)
  
  season <- bases %>%
    spread(key = base, value = value) %>%
    dplyr::select(-day) %>%
    as.matrix()
  
  N <- round(mean(fitted_pomp@covar[, "N"]))
  
  R0 <- calc_R0(beta = beta, qis = qis, season = season)
  crit_vacc_cover <- 1 - (1/max(R0))
  
  # Simulate from the new pomp object ---------------------------------------
  
  outsims <- foreach(i = speed_grid,
                     .packages = c("pomp", "tidyverse", "dplyr"), 
                     .combine = "rbind") %dopar%
  {
    years <- 100
    weeks <- years*52
    days <- years*365
    vacc_coverage_ts <- sapply(0:days, FUN = rho_curve_ramp, 
                               start = 50*365, speed = i)
    
    simulator_pomp <- make_pomp_simulator(
      do_city, 
      mles, 
      years_to_sim = years, 
      initial_population_size = round(mean(fitted_pomp@covar[, "N"])), 
      susc_discount = 1,
      vacc_coverage_ts = vacc_coverage_ts
    )
    
    model_sims <- simulate(
      simulator_pomp,
      nsim = 1,
      as.data.frame = TRUE,
      include.data = FALSE) %>%
      as_tibble() 
    
    # summ <- model_sims %>%
    #   filter(time > 0) %>%
    #   mutate(year = trunc(time)) %>%
    #   group_by(year) %>%
    #   summarise(avg_re = mean(RE_seas))
    par(mfrow = c(1, 2))
    plot(model_sims$S+model_sims$I+model_sims$E+model_sims$R, type = "l")
    vacc_start <- 50*365/7
    tcrit <- which(vacc_coverage_ts >= crit_vacc_cover)[1]/7
    window_start <- vacc_start - (tcrit - vacc_start)
    plot(model_sims$reports, type = "l", xlab = "week", ylab = "reports", col = "grey45")
    abline(v = vacc_start, col = "red", lwd =2, lty = 2)
    abline(v = tcrit, col = "dodgerblue4", lty = 2, lwd = 2)
    abline(v = window_start, col = "dodgerblue4", lty = 2, lwd = 2)
    # par(mfrow = c(1,2))
    # acf(model_sims$reports[window_start:vacc_start])
    # acf(model_sims$reports[vacc_start:tcrit])
    # plot(summ$avg_re, type = "l", col = "grey35", xlab = "year", ylab = expression(R[E]))
    # abline(h = 1, col = "red", lty = 2, lwd = 2)
    # abline(v = 50, col = "dodgerblue4", lty = 2, lwd = 2)
    # 
    # 
    # plot_sim <- model_sims %>%
    #   filter(sim == 1 & time > 0)
    # weekly_vacc <- vacc_coverage_ts[seq(1, length(vacc_coverage_ts), 7)]
    # 
    # plot(plot_sim$time, plot_sim$RE_seas, type = "l", col = "grey")
    # lines(plot_sim$time, predict(loess(RE_seas~time, data = plot_sim, span = 0.75)), col = "red")
    # abline(h = 1)
    # abline(v = 5)
    # 
    # par(mar = c(5, 4, 4, 4) + 0.1)
    # plot(plot_sim$time, plot_sim$reports, type = "l", xlab = "Time (years)", ylab = "Reported infections")
    # par(new = TRUE)
    # plot(plot_sim$time, weekly_vacc[2:length(weekly_vacc)], xlab = "", ylab = "", xlim = par("usr")[1:2], xaxs = "i",
    #      ylim = c(0.7, 1), type = "l", lwd = 1, axes = FALSE, col = "red")
    # axis(4)
    # mtext("Vaccination coverage", side = 4, line = 2.8)
    
    tmp_re_sims <- model_sims %>%
      dplyr::select(sim, time, RE_seas, reports, vacc_discount) %>%
      mutate(
        vacc_coverage = 1 - vacc_discount,
        vacc_speed = i,
        city = do_city
      ) %>%
      dplyr::select(-vacc_discount)

    outfile <- paste0("../simulations/elimination-simulations-grid-",
                      do_city, "-", i, ".RDS")
    saveRDS(object = tmp_re_sims, file = outfile)
  }
  
}
