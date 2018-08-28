# probe-matching.R:
#  Script to conduct initial parameter estimation based on probe matching,
#  where a synthetic likelihood is maximized. The synthetic likelihood is 
#  made up of statistical summaries of the data set rather than the data
#  points themselves. Probe matching is especially useful for delineating
#  the MLE parameter space when working with nonlinear (or chaotic) dynamics.
#  See Wood 2010, Nature: https://www.nature.com/articles/nature09319.
#
#  This script runs a probe matching algorithm using pomp, and then saves
#  the estimated parameter values as starting points for the global search
#  that is implemented next (on a high performance cluster).
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)


# Load pomp object --------------------------------------------------------

measles_pomp <- readRDS("measles-pomp-object.RDS")



# Define probes for synthetic likelihood ----------------------------------

probe_zeroes <- function(y){
  # number of zeros in the data
  xy <- y["reports", ]
  as.numeric(length(which(xy == 0)))
}

probe_max <- function(y){
  # max incidence in the data
  max(y["reports", ], na.rm = TRUE)
}

probe_cusum <- function(y){
  # total number of cases in the data
  cases <- y["reports", ]
  max(cumsum(cases))
}

plist <- list(
  probe_zeroes,
  probe_cusum,
  probe_max,
  probe.acf("reports", lags = 1,transform = sqrt, type = "correlation")
)


# Run probe matching algorithm in pomp ------------------------------------

param_names <- names(coef(measles_pomp))

tic <- Sys.time()

pm <- probe.match(
  measles_pomp,
  probes = plist,
  est = param_names,
  nsim = 10,
  transform = TRUE,
  start = coef(measles_pomp),
  method = "SANN",
  maxit = 100
)

toc <- Sys.time()
toc-tic


simulate(
  measles_pomp,
  params = coef(pm),
  nsim = 19,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = reports, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())


