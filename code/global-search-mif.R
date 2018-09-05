# global-search-mif.R
#  This script conducts a global parameter search in a very large parameter
#  parameter space for the MLE of the SIR model, given the measles data.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

# On PC
# library(tibble)
# library(magrittr)
# library(pomp)
# library(foreach)
# library(doParallel)
# library(dplyr)

# On HPC
library("tibble", lib.loc="~/myRlib/")
library("magrittr", lib.loc="~/myRlib/")
library("pomp", lib.loc="~/myRlib/")
library("foreach", lib.loc="~/myRlib/")
library("doParallel", lib.loc="~/myRlib/")
library("dplyr", lib.loc="~/myRlib/")





# Load pomp object --------------------------------------------------------

measles_pomp <- readRDS("measles-pomp-object.RDS")
start_population <- as.numeric(measles_pomp@covar[1,1])


# Define Latin hypercube parameter space ----------------------------------

grid_size <- 5000

param_uppers <- tibble(
  beta_mu = 50,
  gamma = 50,
  beta_sd = 2,
  b1 = 10,
  b2 = 10,
  b3 = 10,
  b4 = 1,
  b5 = 10,
  b6 = 1,
  iota = 10,
  rho = 0.8,
  S_0 = 100000,
  I_0 = 500,
  R_0 = start_population - 100000 - 500
) %>%
  as.numeric()
names(param_uppers) <- names(coef(measles_pomp))

param_lowers <- tibble(
  beta_mu = 0.001,
  gamma = 0.001,
  beta_sd = 0.00001,
  b1 = -10,
  b2 = -10,
  b3 = -10,
  b4 = -10,
  b5 = -10,
  b6 = -10,
  iota = 0.001,
  rho = 0.1,
  S_0 = 1000,
  I_0 = 10,
  R_0 = start_population - 1000 - 10
) %>%
  as.numeric()
names(param_lowers) <- names(coef(measles_pomp))

guesses <- sobolDesign(
  lower = param_lowers,
  upper = param_uppers,
  nseq = grid_size
) 


# Plot and save the parameter space ---------------------------------------

# pdf("../results/lhs_global.pdf", width = 10, height = 10)
# plot(guesses, pch = ".", las = 1)
# dev.off()


# Perform initial parameter search with pfilter ---------------------------

particles <- 20
cores <- 4
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterCall(cl, function(x) .libPaths(x), .libPaths("~/myRlib"))

if(file.exists("initial-search-lls.RDS") == FALSE){

  foreach(
    guess = iter(guesses[1:4,],"row"),
    .combine = rbind,
    .packages = c("pomp","magrittr","dplyr"),
    .errorhandling = "remove",
    .export = "measles_pomp",
    .inorder = FALSE
  ) %dopar% {
    ll <- logmeanexp(replicate(10, logLik(pfilter(measles_pomp, Np = particles, params = guess))), se=TRUE)
    data.frame(
      loglik = ll[1],
      loglik_se = ll[2]
    ) %>%
      bind_cols(guess)
     
  } -> initial_lls

  saveRDS(object = initial_lls, file = "initial-search-lls.RDS")

} else{
  warning("Output file already exists.")
}




# Perform initial MIF search ----------------------------------------------

# tic <- Sys.time()
# ll <- logmeanexp(replicate(10, logLik(pfilter(measles_pomp, params = guesses[1,], Np = 2000))), se=TRUE)
# Sys.time()-tic # 40.5 seconds
# 
# 
# tic <- Sys.time()
# measles_pomp %>% 
#   mif2(
#     start = unlist(guesses[1,]),
#     Nmif = 25,
#     Np = 200,
#     transform = TRUE,
#     cooling.fraction.50 = 1,
#     cooling.type = "geometric",
#     rw.sd = rw.sd(
#       beta_mu = 0.02, 
#       gamma = 0.02,
#       rho = 0.02, 
#       tau = 0.02, 
#       beta_sd = 0.02,
#       b1 = 0.02, 
#       b2 = 0.02, 
#       b3 = 0.02, 
#       b4 = 0.02, 
#       b5 = 0.02, 
#       b6 = 0.02,
#       I_0 = ivp(0.1, 52), 
#       S_0 = ivp(0.1, 52),
#       R_0 = ivp(0.1, 52),
#       iota = 0.02
#     )
#   ) -> mf
# 
# toc <- Sys.time()
# timeit <- toc-tic
# 
# saveRDS(mf, "testmf.RDS")
# saveRDS(timeit, "proctime.RDS")

# 
# if(file.exists("global-search.RDS") == FALSE){
#   
#   foreach(
#     guess = iter(guesses[1:10,],"row"),
#     .combine = rbind,
#     .packages = c("pomp","magrittr"),
#     .errorhandling = "remove",
#     .export = "measles_pomp",
#     .inorder = FALSE
#   ) 
#   %dopar% {
#     measles_pomp %>% 
#       mif2(
#         start = unlist(guess),
#         Nmif = 25,
#         Np = 2000,
#         transform = TRUE,
#         cooling.fraction.50 = 1,
#         cooling.type = "geometric",
#         rw.sd = rw.sd(
#           beta_mu = 0.02, 
#           gamma = 0.02,
#           rho = 0.02, 
#           tau = 0.02, 
#           beta_sd = 0.02,
#           b1 = 0.02, 
#           b2 = 0.02, 
#           b3 = 0.02, 
#           b4 = 0.02, 
#           b5 = 0.02, 
#           b6 = 0.02,
#           I_0 = ivp(0.1, 52), 
#           S_0 = ivp(0.1, 52),
#           R_0 = ivp(0.1, 52),
#           iota = 0.02
#         )
#       ) -> mf
#     
#     ll <- logmeanexp(replicate(10,logLik(pfilter(mf))), se=TRUE)
#     tibble(
#       loglik = ll[1],
#       loglik_se = ll[2], 
#       as.list(coef(mf))
#     )
#   } -> global_mles
#   
#   saveRDS(object = global_mles, file = "global-search.RDS")
#   
# } else{
#   global_mles <- readRDS("global-search.RDS")
# }
