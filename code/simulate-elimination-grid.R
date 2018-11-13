# simulate-emergence.R
#  Script to simulate the elimination of measles in four Nigerien cities
#  based on empirically-fitted parameter values. Simulations are driven by
#  a slowing decreasing susceptible pool due to vaccination.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)



# Define function for vaccination campaign --------------------------------

rho_curve_ramp <- function(t, start = 52*4, speed = -0.015){
  ifelse(
    t <= start,
    rho <- 0,
    rho <- 1 * (1 - exp((t - start) * speed))
  )
  return(rho)
}


# Make example plot of vaccination curve ----------------------------------

test <- sapply(0:520, FUN = rho_curve_ramp)
pdf(file = "../figures/vaccination-coverage-example.pdf", width = 6, height = 4.5)
plot(test, type = "l", xlab = "Time (weeks)", ylab = "Vaccination coverage", col = "dodgerblue4", lwd = 2)
abline(h = 0.7, lty = 2)
abline(h = 0.95, lty = 3)
text(120, 0.65, labels = "Current vaccination coverage in Niger", cex = 0.7)
text(125, 0.9, labels = "Vaccination coverage for herd immunity", cex = 0.7)
dev.off()


# Set up grid of vaccination roll out speeds ------------------------------

speed_grid <- seq(-0.01, -0.001, length.out = 6)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)  # functions for parallel computing

registerDoParallel()
source("make-pomp-simulator-function.R")


all_sims <- tibble()

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model -----------------------------------
  
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se, -beta_mu)
  
  mle_beta <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    pull(beta_mu)
  
  pomp_file <- paste0("./measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  # Simulate from the new pomp object ---------------------------------------
  
  outsims <- foreach(i = speed_grid, .packages = c("pomp", "tidyverse", "dplyr"), .combine = "rbind") %dopar%
  {
    years <- 100
    weeks <- years*52
    days <- years*365
    vacc_coverage_ts <- sapply(0:days, FUN = rho_curve_ramp, start = 5*365, speed = i)
    vacc_coverage_ts[vacc_coverage_ts < 0.7] <- 0.7
    
    simulator_pomp <- make_pomp_simulator(
      do_city, 
      mles, 
      years_to_sim = years, 
      initial_population_size = fitted_pomp@covar[1, "N"], 
      susc_discount = 1,
      vacc_coverage_ts = vacc_coverage_ts
    )
    
    model_sims <- simulate(
      simulator_pomp,
      nsim = 10,
      as.data.frame = TRUE,
      include.data = FALSE) %>%
      as_tibble()
    

    plot_sim <- model_sims %>%
      filter(sim == 1 & time > 0)
    weekly_vacc <- vacc_coverage_ts[seq(1, length(vacc_coverage_ts), 7)]
    
    plot(plot_sim$time, plot_sim$RE_seas, type = "l", col = "grey")
    lines(plot_sim$time, predict(loess(RE_seas~time, data = plot_sim, span = 0.05)), col = "red")
    abline(h = 1)
    abline(v = 5)
    
    par(mar = c(5, 4, 4, 4) + 0.1)
    plot(plot_sim$time, plot_sim$reports, type = "l", xlab = "Time (years)", ylab = "Reported infections")
    par(new = TRUE)
    plot(plot_sim$time, weekly_vacc, xlab = "", ylab = "", xlim = par("usr")[1:2], xaxs = "i",
         ylim = c(0.7, 1), type = "l", lwd = 1, axes = FALSE, col = "red")
    axis(4)
    mtext("Vaccination coverage", side = 4, line = 2.8)
    
    summ <- model_sims %>%
      mutate(year = trunc(time)) %>%
      group_by(year) %>%
      summarise(avg_re = max(RE_seas))
    
    plot(summ$avg_re, type = "l")
    
    tmp_re_sims <- model_sims %>%
      dplyr::select(sim, time, RE_seas, reports) %>%
      mutate(
        city = do_city,
        susc_discount = i
      )
    
    outfile <- paste0("../simulations/emergence-simulations-grid-", do_city, "-", i, ".RDS")
    saveRDS(object = tmp_re_sims, file = outfile)
  }
  
}
