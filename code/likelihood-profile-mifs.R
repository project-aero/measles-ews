# likelihood-profile-mifs.R
#  Script to run profiles over the likelihood surface for selected parameters
#  to calculate confidence intervals and to investigate the robustness
#  of estimates maximizimum likelihood estimates from initial MIF runs.
#  
# NOTE:
#  This script is primarily designed to be run on a high performance cluster.
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
  slice(1:2000)


# Make grid for profile ---------------------------------------------------

min_value <- range(mles$rho)[1] - (range(mles$rho)[1]*0.5)
max_value <- range(mles$rho)[2] + (range(mles$rho)[2]*0.5)
profile_grid <- seq(from = min_value, to = max_value, length.out = 1000)

# S0_grid <- seq(from = 0.1, to = 0.4, length.out = 50)
# profile_grid <- expand.grid(rho_grid, S0_grid)
# colnames(profile_grid) <- c("rho", "S_0") 


# Combine profile grid with MLEs ------------------------------------------

highest_mles <- mles %>%
  filter(loglik == max(loglik)) %>%
  bind_rows(replicate(length(profile_grid)-1, ., simplify = FALSE))

profile_params <- highest_mles %>%
  dplyr::select(-do_grid, -loglik, -loglik_se) %>%
  mutate(
    rho = profile_grid
  ) %>%
  dplyr::select(
    beta_mu,beta_sd,b1,b2,b3,b4,b5,b6,iota,rho,S_0,E_0,I_0,tau
  ) %>%
  as.matrix()

colnames(profile_params) <- names(coef(measles_pomp))


# Perform MIF -------------------------------------------------------------

particles <- 10000
mif_iters <- 100

mf <- measles_pomp %>% 
  mif2(
    start = unlist(profile_params[do_grid,]),
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 1,
    cooling.type = "geometric",
    rw.sd = rw.sd(
      beta_mu = 0.02,
      beta_sd = 0.02,
      iota = 0.02,
      rho = 0.0,  # rho does not vary
      b1 = 0.02,
      b2 = 0.02,
      b3 = 0.02,
      b4 = 0.02,
      b5 = 0.02,
      b6 = 0.02,
      I_0 = ivp(0.1),
      E_0 = ivp(0.1),
      S_0 = ivp(0.1),
      tau = 0.02
    )
  ) %>%
  mif2(
    start = unlist(profile_params[do_grid,]),
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 0.9,
    cooling.type = "geometric",
    rw.sd = rw.sd(
      beta_mu = 0.02,
      beta_sd = 0.02,
      iota = 0.02,
      rho = 0.0,  # rho does not vary
      b1 = 0.02,
      b2 = 0.02,
      b3 = 0.02,
      b4 = 0.02,
      b5 = 0.02,
      b6 = 0.02,
      I_0 = ivp(0.1),
      E_0 = ivp(0.1),
      S_0 = ivp(0.1),
      tau = 0.02
    )
  )

ll <- logmeanexp(replicate(50, logLik(pfilter(mf, Np = particles))), se=TRUE)

outdf <- data.frame(
  do_grid = do_grid,
  loglik = as.numeric(ll[1]),
  loglik_se = as.numeric(ll[2]),
  rho_value = profile_grid[do_grid]
) 

out_file <- "rho-profile.csv"
write.table(outdf, out_file, sep = ",", col.names = F, append = T, row.names = FALSE)

