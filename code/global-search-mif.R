# global-search-mif.R
#  This script conducts a global parameter search in a very large parameter
#  parameter space for the MLE of the SIR model, given the measles data.
#
# Author:
#  Andrew Tredennick


args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)


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


# # Set up output files -----------------------------------------------------
# 
# mles <- data.frame(
#   do_grid = NA,
#   loglik = NA,
#   loglik_se = NA,
#   beta_mu = NA,
#   beta_sd = NA,
#   b1 = NA,
#   b2 = NA,
#   b3 = NA,
#   b4 = NA,
#   b5 = NA,
#   b6 = NA,
#   iota = NA,
#   S_0 = NA,
#   I_0 = NA,
#   R_0 = NA
# )
# 
# ll_file <- "initial-mif-lls.csv"
# write.table(mles, ll_file, sep = ",", col.names = T, append = T, row.names = FALSE)
# 
# mf_traces <- data.frame(
#   do_grid = NA,
#   iteration = NA,
#   loglik = NA,
#   nfail = NA,
#   beta_mu = NA,
#   beta_sd = NA,
#   b1 = NA,
#   b2 = NA,
#   b3 = NA,
#   b4 = NA,
#   b5 = NA,
#   b6 = NA,
#   iota = NA,
#   S_0 = NA,
#   I_0 = NA,
#   R_0 = NA
# )
# 
# trace_file <- "initial-mif-traces.csv"
# write.table(mf_traces, trace_file, sep = ",", col.names = T, append = T, row.names = FALSE)


# Load pomp object --------------------------------------------------------

measles_pomp <- readRDS("measles-pomp-object.RDS")
start_population <- as.numeric(measles_pomp@covar[1,1])


# Define Latin hypercube parameter space ----------------------------------

# Number of random parameters to generate
grid_size <- 1000

# Upper bounds for random parameters
param_uppers <- tibble(
  beta_mu = 1000,
  beta_sd = 5,
  b1 = 3,
  b2 = 3,
  b3 = 3,
  b4 = 3,
  b5 = 3,
  b6 = 3,
  iota = 50,
  S_0 = 0.16,
  I_0 = 0.0008
) %>%
  as.numeric()
names(param_uppers) <- names(coef(measles_pomp))[1:11]

# Lower bounds for random parameters
param_lowers <- tibble(
  beta_mu = 26,
  beta_sd = 0.00001,
  b1 = -3,
  b2 = -3,
  b3 = -3,
  b4 = -3,
  b5 = -3,
  b6 = -3,
  iota = 0.001,
  S_0 = 0.00016,
  I_0 = 0.000016
) %>%
  as.numeric()
names(param_lowers) <- names(coef(measles_pomp))[1:11]

# Construct random latin hypercube sample
set.seed(123471246)  # get same LHS every time
lhs_grid <- randomLHS(n = grid_size, k = length(coef(measles_pomp)))

for(i in 1:ncol(lhs_grid)){
  lhs_grid[ , i] <- qunif(lhs_grid[ , i], param_lowers[i], param_uppers[i])
}

colnames(lhs_grid) <- names(coef(measles_pomp))
lhs_grid[,12] <- 1 - (lhs_grid[,10] + lhs_grid[,11])

# Perform initial MIF -----------------------------------------------------

particles <- 2000
mif_iters <- 50

mf <- measles_pomp %>% 
  mif2(
    start = unlist(lhs_grid[do_grid,]),
    Nmif = mif_iters,
    Np = particles,
    transform = TRUE,
    cooling.fraction.50 = 1,
    cooling.type = "geometric",
    rw.sd = rw.sd(
      beta_mu = 0.02,
      beta_sd = 0.02,
      iota = 0.02,
      b1 = 0.02,
      b2 = 0.02,
      b3 = 0.02,
      b4 = 0.02,
      b5 = 0.02,
      b6 = 0.02,
      I_0 = ivp(0.1),
      S_0 = ivp(0.1),
      R_0 = ivp(0.1)
    )
  ) 

ll <- logmeanexp(replicate(10,logLik(pfilter(mf, Np = particles))), se=TRUE)
coef_ests <- data.frame(t(coef(mf)))

outdf <- data.frame(
  do_grid = do_grid,
  loglik = as.numeric(ll[1]),
  loglik_se = as.numeric(ll[2])
) %>%
  bind_cols(coef_ests)


# Write results to file ---------------------------------------------------

ll_file <- "initial-mif-lls.csv"
write.table(outdf, ll_file, sep = ",", col.names = F, append = T, row.names = FALSE)

outmif <- data.frame(
  do_grid = rep(do_grid, mif_iters+1),
  iteration = seq(0:mif_iters)
) %>%
  bind_cols(as.data.frame(conv.rec(mf)))

trace_file <- "initial-mif-traces.csv"
write.table(outmif, trace_file, sep = ",", col.names = F, append = T, row.names = FALSE)

saveRDS(object = mf, file = paste0("./mif-objects/mifobject-", do_grid, ".RDS"))





# Perform initial parameter search with pfilter ---------------------------

# if(file.exists("initial-search-lls.RDS") == FALSE){
#   
#   particles <- 2000
#   cores <- 20
#   cl <- makeCluster(cores)
#   registerDoParallel(cl)
#   
#   clusterCall(cl, function(x) .libPaths(x), .libPaths("~/myRlib"))
# 
#   foreach(
#     guess = iter(guesses,"row"),
#     .combine = rbind,
#     .packages = c("pomp","magrittr","dplyr"),
#     .errorhandling = "remove",
#     .export = "measles_pomp",
#     .inorder = FALSE
#   ) %dopar% {
#     ll <- logmeanexp(replicate(10, logLik(pfilter(measles_pomp, Np = particles, params = guess))), se=TRUE)
#     data.frame(
#       loglik = ll[1],
#       loglik_se = ll[2]
#     ) %>%
#       bind_cols(guess)
#      
#   } -> initial_lls
# 
#   saveRDS(object = initial_lls, file = "initial-search-lls.RDS")
# 
# } else{
#   warning("Output file already exists.")
# }




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
