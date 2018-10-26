# plot-ews-re-simulated-corrs.R
#  Script to plot the correlations between EWS and effective R over the
#  replicate simulations from the model.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(ggthemes)


# Load simulation results -------------------------------------------------

infile <- paste0("../results/sim-corrs-ews-re.RDS")
all_corrs <- readRDS(file = infile) %>%
  mutate(
    ews = ifelse(ews == "variance", "Variance", ews),
    ews = ifelse(ews == "variance_first_diff", "First differenced variance", ews),
    ews = ifelse(ews == "autocovariance", "Autocovariance", ews),
    ews = ifelse(ews == "autocorrelation", "Autocorrelation", ews),
    ews = ifelse(ews == "decay_time", "Decay time", ews),
    ews = ifelse(ews == "mean", "Mean", ews),
    ews = ifelse(ews == "index_of_dispersion", "Index of dispersion", ews),
    ews = ifelse(ews == "coefficient_of_variation", "Coefficient of variation", ews),
    ews = ifelse(ews == "skewness", "Skewness", ews),
    ews = ifelse(ews == "kurtosis", "Kurtosis", ews)
  )


# Make the plot -----------------------------------------------------------

ggplot(all_corrs, aes(x = ews, y = spearman_value, fill = ews)) +
  geom_hline(aes(yintercept = 0), color = "darkred", linetype = 1) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~city, ncol = 1) +
  labs(y = expression(paste("Spearman's ",rho)), x = "Early warning signal") +
  scale_fill_manual(values = ptol_pal()(10), guide = FALSE) +
  coord_flip() +
  theme_minimal()
