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

DO_CITY <- "Niamey"  # which city to model

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)


# Load MLEs ---------------------------------------------------------------

mle_file <- paste0("../results/initial-mif-lls-", DO_CITY, ".csv")
mles <- read_csv(mle_file) %>%
  slice(2:n()) %>%
  arrange(-loglik)