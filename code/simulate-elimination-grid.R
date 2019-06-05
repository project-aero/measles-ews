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


# Set up grid of vaccination roll out speeds ------------------------------

# speed_grid <- seq(-5e-04, -1e-04, length.out = 6)  # exponential
speed_grid <- seq(0.000015, 0.0001, length.out = 6)  # linear


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)  # functions for parallel computing

if(parallel::detectCores() <= length(speed_grid)){
  registerDoParallel(cores = detectCores() - 1)
} else{
  registerDoParallel(cores = length(speed_grid))
}

source("make-pomp-simulator-function.R")

for(do_city in c("Niamey")){
  
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
  
  
  # Simulate from the new pomp object ---------------------------------------
  
  outsims <- foreach(i = speed_grid,
                     .packages = c("pomp", "tidyverse", "dplyr"), 
                     .combine = "rbind") %do%
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
      exposed_discount = 1,
      infected_discount = 1,
      vacc_coverage_ts = vacc_coverage_ts
    )
    
    model_sims <- simulate(
      simulator_pomp,
      nsim = 500,
      as.data.frame = TRUE,
      include.data = FALSE,
      seed = 172856) %>%
      as_tibble() 
    
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
    #saveRDS(object = tmp_re_sims, file = outfile)
  }
  
}


# Make maps of skips to S

prevacc <- model_sims %>% filter(time >= 20 & time < 50)
postvac <- model_sims %>% filter(time >= 50 & time < 80)

presplit_sim <- split(prevacc, prevacc$sim)
postsplit_sim <- split(postvac, postvac$sim)

simsum <- function(sm){
  Slo <- loess(S~time, data = sm, span = 52 / (nrow(sm) * 2))
  sm$Ssmooth <- predict(Slo)
  sdrop <- diff(sm$Ssmooth) < 0
  is_epi <- c(FALSE, sdrop)
  epirle <- rle(is_epi)
  is_scrap <- epirle$lengths[!epirle$values] < 5
  epirle$values[!epirle$values][is_scrap] <- TRUE 
  epirle$values[!epirle$values] <- 1 + seq(1, sum(!epirle$values))
  iepids <- inverse.rle(epirle)
  
  miep <- max(iepids)
  splt <- split(sm$S, iepids)
  splt_time <- split(sm$time, iepids)
  trimsplt <- splt[-miep][-1]
  trimsplt_time <- splt_time[-miep][-1]
  
  iep_lens <- sapply(trimsplt, length)
  iep_s0 <- sapply(trimsplt, "[", 1)
  iep_sfin <- mapply(function(x, y) x[y], x = trimsplt, y = iep_lens)
  iep_tfin <- mapply(function(x, y) x[y], x = trimsplt_time, y = iep_lens)
  epi_sizes <- c(iep_sfin[-length(iep_lens)] - iep_s0[-1], NA)
  data.frame(s0 = iep_s0, sfin = iep_sfin, tfin = iep_tfin, 
             iep_weeks = iep_lens, epi_sizes = epi_sizes)
}


presimsums <- lapply(presplit_sim, simsum)  
postsimsums <- lapply(postsplit_sim, simsum)

premapdata <- do.call(rbind, presimsums)
postmapdata <- do.call(rbind, postsimsums)

plot(iep_weeks ~ s0, data = premapdata)
points(iep_weeks ~s0, data = postmapdata, col = 2)

premapdata %>% filter(s0 < 25000 & s0 > 20000) %>% pull("iep_weeks") %>% density %>% plot
postmapdata %>% filter(s0 < 25000 & s0 > 20000) %>% pull("iep_weeks") %>% density %>% plot

# Make map of epi size to S and Time of start

plot(epi_sizes ~ sfin, data = premapdata, xlim = c(2e4, 8e4), ylim = c(0, 7e4))
plot(epi_sizes ~ sfin, data = postmapdata, xlim = c(2e4, 8e4), ylim = c(0, 7e4))

premapdata$season <- premapdata$tfin + 0.5 - trunc(premapdata$tfin + 0.5)
postmapdata$season <- postmapdata$tfin + 0.5 - trunc(postmapdata$tfin + 0.5)

plot(epi_sizes ~ season, data = postmapdata)
plot(epi_sizes ~ season, data = premapdata)

# Compare relationship of epi size and time of year by vaccination program

premapdata$vacc <- 0
postmapdata$vacc <- 1

mapdata <- rbind(premapdata, postmapdata)

mapgam <- mgcv::gam(epi_sizes ~ s(sfin) + s(season) + vacc, data = mapdata)
anova(mapgam)

### The increasing vaccination rates reduce epidemic sizes but this effect is 
### much smaller than S0 and seasonality. Leaving it out reduces the deviance negligably (<1% change).

library(ggplot2)

g <- ggplot(data = mapdata, aes(x = season, y = sfin)) + geom_point() + 
  geom_density2d() + facet_wrap(~vacc)
plot(g)

## The distribution of Sfin is shifted up a few thousand in the initial vaccinatino program.

g <- ggplot(data = mapdata, aes(x = iep_weeks)) + geom_histogram() + facet_grid(vacc~.)
g

mapdata %>% group_by(vacc) %>% summarise(meaniep = mean(iep_weeks))

## There are more 1-year ieps in the initial vaccination regime and the mean iep is acctualy slightly larger.

simsplt <- split(model_sims, model_sims$sim)

filtcalc <- function(df, filter_len = 10){
  df$year <- round(df$time)
  splt <- split(df, df$year)
  yrd <- sapply(splt, function(x) sum(x$cases))
  yf <- stats::filter(yrd, filter = filter_len:1, sides = 1,
                      method = "convolution")
  yr <- as.integer(names(splt))
  data.frame(year = yr[-1], tot_cases = yrd[-1], max_cases = yrd[1],
             filt = as.numeric(yf[-length(yf)]))
}

simproc <- lapply(simsplt, filtcalc)
simyr <- bind_rows(simproc, .id = "sim")

simyr$camp <- 1
simyr$camp[simyr$year < 50] <- 0
simyrwin <- simyr %>% filter(year > 20 & year < 80)

filtstat <- function(df){
  m <- lm(I(tot_cases - max_cases) ~ 0 + filt + filt:camp, data = df)
  m2 <- lm(tot_cases ~ camp, data = df)
  c(coef(m)["filt:camp"], coef(m2)["camp"])
}

splityr <- split(simyrwin, simyrwin$sim)
inhibs <- sapply(splityr, filtstat)

rowMeans(inhibs < 0)

# Extra code --------

# Define function to calculate R0 from seasonal params

# calc_R0 <- function(beta, qis, season, eta = (365/8), 
#                     mu = 0.05, nu = 0.05, gamma = (365/5), p = 0.7){
#   B <- as.numeric((1 + exp(season %*% qis)) * beta)
#   R0 <- (eta*B*mu*) / (nu*(eta + nu)*(gamma + nu))
#   return(R0)
# }


# Make example plot of vaccination curve 

# R0 <- calc_R0(beta = beta, qis = qis, season = season)
# crit_vacc_cover <- 1 - (1/max(R0))

# vacc_ex_tbl <- tibble(
#   week = 0:520,
#   coverage = sapply(0:520, FUN = rho_curve_ramp)
# )
# 
# ggplot(vacc_ex_tbl, aes(week, coverage)) +
#   geom_line(color = "dodgerblue4", size = 1) +
#   geom_hline(aes(yintercept = 0.7), linetype = 2, 
#              size = 0.25, color = "grey45") +
#   geom_hline(aes(yintercept = 0.95), linetype = 2, 
#              size = 0.25, color = "grey45") +
#   labs(x = "Week", y = "Vaccination coverage") +
#   annotate(geom = "text", x = 420, y = 0.71, color = "grey45",
#            label = "Current vaccination coverage in Niger", size = 3) +
#   annotate(geom = "text", x = 125, y = 0.94, color = "grey45",
#            label = "Vaccination coverage for herd immunity", size = 3) +
#   theme_classic(base_size = 14)
# ggsave(filename =  "../figures/vaccination-coverage-example.pdf", 
#        width = 6, height = 3.5, units = "in")


# summ <- model_sims %>%
#   filter(time > 0) %>%
#   mutate(year = trunc(time)) %>%
#   group_by(year) %>%
#   summarise(avg_re = mean(RE_seas))
# par(mfrow = c(1, 2))
# plot(model_sims$S+model_sims$I+model_sims$E+model_sims$R, type = "l")
# vacc_start <- 50*365/7
# tcrit <- which(vacc_coverage_ts >= crit_vacc_cover)[1]/7
# window_start <- vacc_start - (tcrit - vacc_start)
# plot(model_sims$reports, type = "l", xlab = "week", ylab = "reports", col = "grey45")
# abline(v = vacc_start, col = "red", lwd =2, lty = 2)
# abline(v = tcrit, col = "dodgerblue4", lty = 2, lwd = 2)
# abline(v = window_start, col = "dodgerblue4", lty = 2, lwd = 2)
# par(mfrow = c(1,2))
# stats::var(model_sims$reports[window_start:vacc_start])
# stats::var(model_sims$reports[vacc_start:tcrit])
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
