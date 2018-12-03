# bootstrap-fit-mif.R
#  Script to fit the SEIR model to replicate (bootstrapped) time series
#  for each city. Models were simulated at the MLE to produce 100 replicates.
#  For each replicate, this script fits the model using MIF from 50
#  unique initial conditions for the parameters.
#
#  NOTE: This script is written for running on a High Performance Cluster,
#        and certain variables are defined from the command line.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)



# Load packages -----------------------------------------------------------

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


# Define computation grid -------------------------------------------------

nboots <- 100
nmifs <- 50
comp_grid <- expand.grid(1:nboots, 1:50)
colnames(comp_grid) <- c("boot_series", "param_set")


# Load pomp object --------------------------------------------------------

do_city <- "Agadez"
measles_pomp <- readRDS(paste0("measles-pomp-object-",do_city, ".RDS"))
start_population <- as.numeric(measles_pomp@covar[1,1])


# Load bootstrapped series ------------------------------------------------

bootstraps <- readRDS(paste0("../simulations/bootstrap-sims-", do_city, ".RDS"))
boot_data <- bootstraps %>%
  filter(sim == comp_grid[do_grid, "boot_series"])
boot_data <- boot_data$data[[1]]$reports  # extract from nested df
boot_reports <- t(as.matrix(c(NA, boot_data)))  # make a row vector for pomp
row.names(boot_reports) <- "reports"


# Redefine pomp data object -----------------------------------------------

measles_pomp@data <- boot_reports


# Define Latin hypercube parameter space ----------------------------------

# Number of random parameters to generate
grid_size <- nmifs

# Upper bounds for random parameters
param_uppers <- tibble(
  beta_mu = 1000,
  beta_sd = 5,
  b1 = 10,
  b2 = 10,
  b3 = 10,
  b4 = 10,
  b5 = 10,
  b6 = 10,
  iota = 500,
  rho = 0.999,
  S_0 = 0.2,
  E_0 = 0.008,
  I_0 = 0.008,
  tau = 50
) %>%
  as.numeric()
names(param_uppers) <- names(coef(measles_pomp))

# Lower bounds for random parameters
param_lowers <- tibble(
  beta_mu = 5,
  beta_sd = 0.00001,
  b1 = -10,
  b2 = -10,
  b3 = -10,
  b4 = -10,
  b5 = -10,
  b6 = -10,
  iota = 0.001,
  rho = 0.001,
  S_0 = 0.00016,
  E_0 = 0.000016,
  I_0 = 0.000016,
  tau = 0.00001
) %>%
  as.numeric()
names(param_lowers) <- names(coef(measles_pomp))

# Construct random latin hypercube sample
set.seed(123471246)  # get same LHS every time
lhs_grid <- randomLHS(n = grid_size, k = length(coef(measles_pomp)))

for(i in 1:ncol(lhs_grid)){
  lhs_grid[ , i] <- qunif(lhs_grid[ , i], param_lowers[i], param_uppers[i])
}

colnames(lhs_grid) <- names(coef(measles_pomp))


# Perform initial MIF -----------------------------------------------------

particles <- 10000
mif_iters <- 100

mf <- measles_pomp %>% 
  mif2(
    start = unlist(lhs_grid[comp_grid[do_grid, "param_set"],]),
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 1,
    cooling.type = "geometric",
    rw.sd = rw.sd(
      beta_mu = 0.02,
      beta_sd = 0.02,
      iota = 0.02,
      rho = 0.02,
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
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 0.9,
    cooling.type = "geometric",
    rw.sd = rw.sd(
      beta_mu = 0.02,
      beta_sd = 0.02,
      iota = 0.02,
      rho = 0.02,
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
coef_ests <- data.frame(t(coef(mf)))

outdf <- data.frame(
  do_grid = do_grid,
  loglik = as.numeric(ll[1]),
  loglik_se = as.numeric(ll[2])
) %>%
  bind_cols(coef_ests)


# Write results to file ---------------------------------------------------

ll_file <- "bootstrap-mif-lls.csv"
write.table(outdf, ll_file, sep = ",", col.names = F, append = T, row.names = FALSE)

outmif <- data.frame(
  do_grid = rep(do_grid, mif_iters+1),
  iteration = seq(0:mif_iters)
) %>%
  bind_cols(as.data.frame(conv.rec(mf)))

trace_file <- "bootstrap-mif-traces.csv"
write.table(outmif, trace_file, sep = ",", col.names = F, append = T, row.names = FALSE)



