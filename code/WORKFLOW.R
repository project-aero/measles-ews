# WORKFLOW.R
#  This script documents the workflow for the `measles-ews` project. The
#  workflow is implemented by sourcing scripts in the `./code/` subdirectory,
#  but several of the scripts are specifically designed for running
#  on a high performance computing cluster; these scripts are denoted
#  with `**HPC**`.
#
#  This workflow should be considered as a guide through the analysis rather
#  than a script that reproduces the analysis on the fly. The analysis, and
#  results can be reproduced, of course, but it will take a while :).
#
# Primary author:
#  Andrew Tredennick (atredenn@gmail.com)
#
# Contributors:
#  Eamon O'Dea
#  Pejman Rohani
# =============================================================================


# Load required packages --------------------------------------------------

library(tidyverse)
library(pomp)
library(ggthemes)
library(spaero)
library(pROC)


# Model fitting scripts ---------------------------------------------------

