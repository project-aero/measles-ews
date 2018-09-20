# calculate-beta-ews-correlations.R
#  Script to calculate the correlations of transmission rate with several
#  EWS over time. 
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(lubridate)
library(spaero)


# Load previous results ---------------------------------------------------

ews_results <- readRDS("../results/ews-niger-cities.RDS") %>%
  filter(city == "Niamey (City)")

beta_results <- read_csv("../results/transmission-posteriors.csv", col_types = cols()) %>%
  dplyr::select(time, med) %>%
  mutate(
    date = unique(ews$date)[2:length(unique(ews$date))]
  ) %>%
  dplyr::select(-time) %>%
  rename(beta_value = med) %>%
  filter(beta_value != max(beta_value, na.rm = TRUE))


scale_it <- function(x){
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

ews_beta <- ews_results %>%
  left_join(beta_results, by = "date") %>%
  group_by(ews) %>%
  mutate(
    scaled_ews = scale_it(value),
    scaled_beta = scale_it(beta_value)
  )

ews_beta_corrs <- ews_beta %>%
  group_by(ews) %>%
  summarise(
    spearman_value = SpearmanRho(value, beta_value, use = c("pairwise.complete.obs"),  conf.level = 0.95)[1],
    spearman_lwr = SpearmanRho(value, beta_value, use = c("pairwise.complete.obs"),  conf.level = 0.95)[2],
    spearman_upr = SpearmanRho(value, beta_value, use = c("pairwise.complete.obs"),  conf.level = 0.95)[3], 
    spearman_pvalue = cor.test(value, beta_value, use = "pairwise.complete.obs", method = "spearman")[["p.value"]]
  ) %>%
  mutate(
    pos = spearman_value > 0,
    sig = spearman_pvalue < 0.05,
    color_id_final = "cnull",
    color_id_final = ifelse(pos == TRUE & sig == TRUE, "apos", color_id_final),
    color_id_final = ifelse(pos == FALSE & sig == TRUE, "bneg", color_id_final)
  )

ggplot(ews_beta_corrs, aes(x = ews, y = spearman_value, color = color_id_final)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_errorbar(aes(ymin = spearman_lwr, ymax = spearman_upr), width = 0.1) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "grey35"), guide = FALSE) +
  labs(y = expression(paste("Spearman's ", rho)), x = "Early warning signal") +
  theme_minimal() +
  coord_flip()

ggplot(ews_beta, aes(x = date)) +
  geom_line(aes(y = scaled_ews, color = "EWS")) +
  geom_line(aes(y = scaled_beta, color = "beta")) +
  facet_wrap(~ews) +
  labs(x = "Date", y = "Early warning signal\nor\nTransmission rate") +
  scale_color_manual(values = ptol_pal()(2), labels = c(expression(paste(beta,"     ")), "EWS"), name = NULL) +
  theme_minimal() +
  theme(legend.position = "top")

# ggplot(ews_beta, aes(x = scaled_ews, y = scaled_beta)) +
#   stat_smooth(method = "lm", se = FALSE, color = "red") +
#   facet_wrap(~ews) +
#   theme_minimal()
