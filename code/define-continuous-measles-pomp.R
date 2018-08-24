# define-continuous-measles-pomp.R:
#  Script to define a 'pomp' model object for the dynamics of measles in
#  Niger. The model is a continuous-time SIR model with a time-varying
#  transmission rate. The model is specified as a set of nonlinear stochastic
#  differential equations (SDEs).
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(pomp)
