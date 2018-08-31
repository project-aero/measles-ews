# simulate-series.R
#  Script to simulate several incidence time series from a model with known
#  parameters. These simulation results are used to test the analysis
#  comparing EWS(t) and R_E(t) through time.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(pomp)
library(tidyverse)
library(lubridate)


# Load pomp object --------------------------------------------------------

pomp_file <- "./measles-pomp-object.RDS"
measles_pomp <- readRDS(pomp_file)


# Simulate time series of incidence ---------------------------------------

sims <- as_tibble(
  simulate(
    measles_pomp,
    nsim = 1,
    as.data.frame = TRUE,
    include.data = FALSE
  )
)