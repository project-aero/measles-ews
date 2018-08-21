# large-mif-search.R
#  R script to be run on a multi-core high performance machine. The script runs
#  through several algorithms for estimating parameters for a measles SIR
#  model fit to data from Niamey, Niger. Fitting is primarily done using
#  maximization of the likelihood via iterated filtering, implemented using the
#  `pomp` package.
#
#  The script runs in three main chunks:
#
#    1.  Initial parameter guesses from probe matching, i.e., fitting to a 
#        synthetic likelihood.
#    2a. A global parameter search from a Latin hypercube of parameter values
#        based on noised-up values from (1).
#    2b. A local parameter search around the MLE based on the highest likelihood
#        parameters from (2a).
#    3.  Particle MCMC for Bayesian inference of latent states, including
#        parameter uncertainty, with values starting at the MLE parameters
#        from (2b).
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)
#  Pejman Rohani


# Load packages -----------------------------------------------------------

library(tidyverse)   # data wrangling
library(lubridate)   # time and date functions
library(pomp)        # fitting state-space models with particle filtering
library(foreach)     # functions for parallel computing
library(doParallel)  # functions for parallel computing

registerDoParallel(cores = 25)
theme_set(theme_minimal())



# Load formatted data and pomp object -------------------------------------

biweek_data <- readRDS("niger_cleaned_biweekly_measles.RDS")

if(file.exists("measles-pomp.RDS")){
  measles_pomp <- readRDS("measles-pomp.RDS")
} else{
  source("define-measles-pomp.R")
  measles_pomp <- readRDS("measles-pomp.RDS")
}

param_names <- names(coef(measles_pomp))


# Probe model with synthetic likelihood -----------------------------------

probe.zeroes <- function(y){
  # number of zeros in the data
  
  xy <- y["cases", ]
  as.numeric(length(which(xy == 0)))
}
probe.max <- function(y){
  # max incidence in the data
  
  max(y["cases", ], na.rm = TRUE)
}
probe.cumsum <- function(y){
  # total number of cases in the data
  
  cases <- y["cases", ]
  max(cumsum(cases))
}

plist <- list(
  probe.zeroes,
  probe.max,
  probe.cumsum,
  probe.acf("cases", lags = 1,transform = sqrt, type = "correlation")
)

stew(
  file = "probe-match-results.rda", {
    pm <- probe.match(
      measles_pomp,
      probes = plist,
      est = param_names,
      nsim = 500,
      transform = TRUE,
      start = coef(measles_pomp),
      method = "Nelder-Mead",
      maxit = 10000
    )
  },
  seed = 290868325,
  kind = "L'Ecuyer"
)


# Set up Latin hypercube for sampling space -------------------------------
search_factor <- 10
uppers <- coef(pm)*search_factor
lowers <- coef(pm)*(1/search_factor)
uppers["rho"] <- 1
uppers["S_0"] <- 100000
lowers["S_0"] <- 2000

guesses <- sobolDesign(
  lower = lowers,
  upper = uppers,
  nseq = 200
) %>%
  as_tibble()


# Run global MLE search in parallel ---------------------------------------

if(file.exists("global-search.RDS") == FALSE){
  
  foreach(
    guess = iter(guesses,"row"),
    .combine = rbind,
    .packages = c("pomp","magrittr"),
    .errorhandling = "remove",
    .export = "measles_pomp",
    .inorder = FALSE
  ) %dopar% {
    measles_pomp %>% 
      mif2(
        start = unlist(guess),
        Nmif = 25,
        Np = 2000,
        transform = TRUE,
        cooling.fraction.50 = 1,
        cooling.type = "geometric",
        rw.sd = rw.sd(
          beta_mu = 0.02, 
          rho = 0.02, 
          tau = 0.02, 
          sigma_env = 0.02,
          b1 = 0.02, 
          b2 = 0.02, 
          b3 = 0.02, 
          b4 = 0.02, 
          b5 = 0.02, 
          b6 = 0.02,
          I_0 = ivp(0.1, 27), 
          S_0 = ivp(0.1, 27),
          psi = 0.02
        )
      ) %>%
      mif2(cooling.fraction.50 = 0.95) -> mf
    
    ll <- logmeanexp(replicate(10,logLik(pfilter(mf))), se=TRUE)
    data.frame(
      loglik = ll[1],
      loglik_se = ll[2], 
      as.list(coef(mf))
    )
  } -> global_mles
  
  saveRDS(object = global_mles, file = "global-search.RDS")
  
} else{
  global_mles <- readRDS("global-search.RDS")
}
 


simulate(
  measles_pomp, 
  nsim = 9,
  params = filter(global_mles, loglik == max(loglik)) %>% dplyr::select(-loglik, -loglik_se) %>% unlist(),
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  theme(strip.text=element_blank())
