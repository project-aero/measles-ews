# calculate-boot-cis.R
#  Script to calculate approximate 95% confidence intervals for all parameters
#  from bootstrapped fits of model simulations.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(dplyr)



# Define computation grid -------------------------------------------------

nboots <- 100
nmifs <- 50
comp_grid <- expand.grid(1:nboots, 1:50)
colnames(comp_grid) <- c("boot_series", "param_set")
comp_grid$do_grid <- 1:nrow(comp_grid)


# Load parameter estimate grid --------------------------------------------

param_cis <- read.csv("../results/bootstrap-mif-lls-Agadez.csv") %>%
  slice(2:n()) %>%  # ignore first row of NAs
  left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
  group_by(boot_series) %>%
  filter(loglik == max(loglik)) %>%
  ungroup() %>%
  dplyr::select(-do_grid, -loglik, -loglik_se, -boot_series, -param_set,
                -b1, -b2, -b3, -b4, -b5, -b6) %>%
  gather(key = parameter, value = mle_value) %>%
  group_by(parameter) %>%
  summarise(
    upper_95_ci = quantile(mle_value, 0.975),
    lower_95_ci = quantile(mle_value, 0.025)
  )

