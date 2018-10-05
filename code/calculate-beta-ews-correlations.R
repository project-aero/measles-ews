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
library(DescTools)

DO_CITY = "Niamey"

# Load previous results ---------------------------------------------------

# Early warning signals
ews_results <- readRDS("../results/ews-niger-cities.RDS") %>%
  filter(str_detect(city, paste0("^", DO_CITY))) %>%  # city starts (^) with DO_CITY
  rename(ews_value = value)

# Modeled states of interest
focal_states <- c(
  "effective_r_nonseasonal",
  "effective_r_seasonal",
  "transmission_rate"
)

# Results from particle filtering
filtered_states <- readRDS(paste0("../results/filtered-states-", DO_CITY, ".RDS")) %>%
  filter(state %in% focal_states) %>%
  unnest() %>%
  dplyr::select(date, med, state) %>%
  rename(state_value = med)


# Define scaling function for plotting EWS and states ---------------------

scale_it <- function(x){
  # Scales values to (0,1) range
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


# Combine EWS and stats to calculate correlations -------------------------

ews_states <- ews_results %>%
  left_join(filtered_states, by = "date") %>%
  group_by(ews, state) %>%
  mutate(
    scaled_ews_value = scale_it(ews_value),
    scaled_state_value = scale_it(state_value)
  )

ews_state_corrs <- ews_states %>%
  group_by(ews, state) %>%
  summarise(
    spearman_value = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[1],
    spearman_lwr = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[2],
    spearman_upr = SpearmanRho(ews_value, state_value, use = "pairwise.complete.obs",  conf.level = 0.95)[3], 
    spearman_pvalue = cor.test(ews_value, state_value, use = "pairwise.complete.obs", method = "spearman")[["p.value"]]
  ) %>%
  mutate(
    pos = spearman_value > 0,
    sig = spearman_pvalue < 0.05,
    color_id_final = "cnull",
    color_id_final = ifelse(pos == TRUE & sig == TRUE, "apos", color_id_final),
    color_id_final = ifelse(pos == FALSE & sig == TRUE, "bneg", color_id_final)
  )

my_labs <- c(
  "effective_r_nonseasonal" = "R (nonseasonal)",
  "effective_r_seasonal" = "R (seasonal)",
  "transmission_rate" = "Transmission rate"
)

ggplot(ews_state_corrs, aes(x = ews, y = spearman_value, color = color_id_final)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_errorbar(aes(ymin = spearman_lwr, ymax = spearman_upr), width = 0.2) +
  geom_point() +
  facet_wrap(~state, labeller = as_labeller(my_labs)) +
  scale_color_manual(values = c(ptol_pal()(2)[1], ptol_pal()(2)[2], "grey35"), guide = FALSE) +
  labs(y = expression(paste("Spearman's ", rho)), x = "Early warning signal") +
  theme_bw() +
  coord_flip() +
  theme(strip.background = element_blank()) 

# ggplot(ews_beta, aes(x = date)) +
#   geom_line(aes(y = scaled_ews, color = "EWS")) +
#   geom_line(aes(y = scaled_beta, color = "beta")) +
#   facet_wrap(~ews) +
#   labs(x = "Date", y = "Early warning signal\nor\nTransmission rate") +
#   scale_color_manual(values = ptol_pal()(2), labels = c(expression(paste(beta,"     ")), "EWS"), name = NULL) +
#   theme_minimal() +
#   theme(legend.position = "top")

# ggplot(ews_beta, aes(x = scaled_ews, y = scaled_beta)) +
#   stat_smooth(method = "lm", se = FALSE, color = "red") +
#   facet_wrap(~ews) +
#   theme_minimal()
