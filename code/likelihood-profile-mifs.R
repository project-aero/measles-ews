# likelihood-profile-mifs.R
#  Script to run profiles over the likelihood surface for selected parameters
#  to calculate confidence intervals and to investigate the robustness
#  of estimates maximizimum likelihood estimates from initial MIF runs.
#  
# NOTE:
#  This script is designed to be run on a high performance cluster.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Get command line argument for do_grid -----------------------------------

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)


# Set city to model -------------------------------------------------------

DO_CITY <- "Niamey"  # which city to model
do_param <- "rho"


# Load libraries ----------------------------------------------------------

# On PC
# library(tibble)
# library(magrittr)
# library(pomp)
# library(foreach)
# library(doParallel)
# library(dplyr)
# library(lhs)

# On HPC
library("tibble", lib.loc="~/myRlib/")
library("magrittr", lib.loc="~/myRlib/")
library("pomp", lib.loc="~/myRlib/")
library("foreach", lib.loc="~/myRlib/")
library("doParallel", lib.loc="~/myRlib/")
library("dplyr", lib.loc="~/myRlib/")
library("lhs", lib.loc="~/myRlib/")


# Read in pomp object -----------------------------------------------------

pomp_file <- paste0("measles-pomp-object-", DO_CITY, ".RDS")
measles_pomp <- readRDS(pomp_file)


# Load MLEs ---------------------------------------------------------------

mle_file <- paste0("initial-mif-lls-", DO_CITY, ".csv")
mles <- read.csv(mle_file) %>%
  slice(2:n()) %>%
  na.omit() %>%
  arrange(-loglik) %>%
  slice(1:4000)  # leave out really out there values


# Make grid for profile ---------------------------------------------------

grid_search_size <- 100

highest_mles <- mles %>%
  filter(loglik == max(loglik)) %>%
  bind_rows(replicate(grid_search_size-1, ., simplify = FALSE))

params <- colnames(mles)[4:ncol(mles)]
params <- params[which(!params %in% c("b1","b2","b3","b4","b5","b6","tau","E_0","I_0","beta_sd"))]

if(do_param == "beta_mu"){
  tmp_values <- pull(mles, var = do_param)
  sd_values <- sd(tmp_values)*1.5
  mu_values <- mean(tmp_values)
  alpha <- mu_values^2 / sd_values^2
  beta <- mu_values / sd_values^2
  set.seed(1234572)  # make sure each worker simulates the same distribution
  tmp_profile <- rgamma(grid_search_size, alpha, beta)
  
  # min_value <- as.numeric(quantile(tmp_values, 0.01))
  # max_value <- as.numeric(quantile(tmp_values, 0.99))
  # tmp_profile <- exp(seq(from = log(min_value), to = log(max_value), length.out = grid_search_size))
  
  tmp_grid <- highest_mles %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  tmp_grid[ , do_param] <- tmp_profile
  tmp_grid <- tmp_grid %>%
    dplyr::select(
      beta_mu,beta_sd,b1,b2,b3,b4,b5,b6,iota,rho,S_0,E_0,I_0,tau
    ) %>%
    mutate(
      profiled_param = do_param
    ) %>%
    bind_rows(replicate(9, ., simplify = FALSE))
  
  large_profile_grid <- tmp_grid
}

if(do_param == "rho"){
  tmp_values <- pull(mles, var = do_param)
  sd_values <- sd(tmp_values)*2
  mu_values <- mean(tmp_values)
  alpha <- (((1-mu_values)/(sd_values^2)) - (1/mu_values)) * mu_values^2
  beta <- alpha*((1/mu_values)-1)
  
  set.seed(1234572)  # make sure each worker simulates the same distribution
  tmp_profile <- rbeta(grid_search_size, shape1 = alpha, shape2 = beta)
  tmp_profile <- seq(min(tmp_profile), max(tmp_profile), length.out = grid_search_size)
  tmp_profile <- round(tmp_profile, 3)
  
  tmp_grid <- highest_mles %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  tmp_grid[ , do_param] <- tmp_profile
  tmp_grid <- tmp_grid %>%
    dplyr::select(
      beta_mu,beta_sd,b1,b2,b3,b4,b5,b6,iota,rho,S_0,E_0,I_0,tau
    ) %>%
    mutate(
      profiled_param = do_param
    ) %>%
    bind_rows(replicate(9, ., simplify = FALSE))
  
  large_profile_grid <- tmp_grid
}

if(do_param == "iota"){
  tmp_values <- pull(mles, var = do_param)
  tmp_values <- tmp_values[which(tmp_values < max(tmp_values))]
  sd_values <- sd(tmp_values)*4
  mu_values <- mean(tmp_values)
  alpha <- mu_values^2 / sd_values^2
  beta <- mu_values / sd_values^2
  set.seed(1234572)  # make sure each worker simulates the same distribution
  tmp_profile <- rgamma(grid_search_size, alpha, beta)
  tmp_profile <- seq(min(tmp_profile), max(tmp_profile), length.out = grid_search_size)
  # tmp_profile <- seq(1e-5, 200, length.out = grid_search_size)
  
  tmp_grid <- highest_mles %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  tmp_grid[ , do_param] <- tmp_profile
  tmp_grid <- tmp_grid %>%
    dplyr::select(
      beta_mu,beta_sd,b1,b2,b3,b4,b5,b6,iota,rho,S_0,E_0,I_0,tau
    ) %>%
    mutate(
      profiled_param = do_param
    ) %>%
    bind_rows(replicate(9, ., simplify = FALSE))
  
  large_profile_grid <- tmp_grid
}

if(do_param == "S_0"){
  tmp_values <- pull(mles, var = do_param)
  sd_values <- sd(tmp_values)*1.5
  mu_values <- mean(tmp_values)
  alpha <- (((1-mu_values)/(sd_values^2)) - (1/mu_values)) * mu_values^2
  beta <- alpha*((1/mu_values)-1)
  
  set.seed(1234572)  # make sure each worker simulates the same distribution
  tmp_profile <- rbeta(grid_search_size, shape1 = alpha, shape2 = beta)
  
  tmp_grid <- highest_mles %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  tmp_grid[ , do_param] <- tmp_profile
  tmp_grid <- tmp_grid %>%
    dplyr::select(
      beta_mu,beta_sd,b1,b2,b3,b4,b5,b6,iota,rho,S_0,E_0,I_0,tau
    ) %>%
    mutate(
      profiled_param = do_param
    ) %>%
    bind_rows(replicate(4, ., simplify = FALSE))
  
  large_profile_grid <- tmp_grid
}


# Perform MIF -------------------------------------------------------------

set.seed(NULL)
particles <- 30000
mif_iters <- 200

profile_params <- large_profile_grid[do_grid, ]
profile_over <- profile_params[ , "profiled_param"]
profile_params <- profile_params[ , which(colnames(profile_params) != "profiled_param")]

if(profile_over == "beta_mu"){
  rw_sd_setup <- rw.sd(beta_mu = 0, beta_sd = 0.02, iota = 0.02, rho = 0.02,
                       b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02,
                       b6 = 0.02, I_0 = ivp(0.1), E_0 = ivp(0.1), 
                       S_0 = ivp(0.1), tau = 0.02)
}

if(profile_over == "iota"){
  rw_sd_setup <- rw.sd(beta_mu = 0.2, beta_sd = 0.02, iota = 0, rho = 0.02,
                       b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02,
                       b6 = 0.02, I_0 = ivp(0.1), E_0 = ivp(0.1), 
                       S_0 = ivp(0.1), tau = 0.02)
}

if(profile_over == "rho"){
  rw_sd_setup <- rw.sd(beta_mu = 0.2, beta_sd = 0.02, iota = 0.02, rho = 0,
                       b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02,
                       b6 = 0.02, I_0 = ivp(0.1), E_0 = ivp(0.1), 
                       S_0 = ivp(0.1), tau = 0.02)
}

if(profile_over == "S_0"){
  rw_sd_setup <- rw.sd(beta_mu = 0.2, beta_sd = 0.02, iota = 0.02, rho = 0.02,
                       b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02,
                       b6 = 0.02, I_0 = ivp(0.1), E_0 = ivp(0.1), 
                       S_0 = ivp(0), tau = 0.02)
}

mf <- measles_pomp %>% 
  mif2(
    start = unlist(profile_params),
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 1,
    cooling.type = "geometric",
    rw.sd = rw_sd_setup
  ) %>%
  mif2(
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 0.9,
    cooling.type = "geometric",
    rw.sd = rw_sd_setup
  )

ll <- logmeanexp(replicate(10, logLik(pfilter(mf, Np = particles))), se=TRUE)

outdf <- data.frame(
  do_grid = do_grid,
  loglik = as.numeric(ll[1]),
  loglik_se = as.numeric(ll[2]),
  value = profile_params[profile_over],
  parameter = profile_over
) 

out_file <- paste0("loglik-profile.csv")
write.table(outdf, out_file, sep = ",", col.names = F, append = T, row.names = FALSE)

