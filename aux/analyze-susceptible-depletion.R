# analyze-susceptible-depletion.R
#  Script to calculate the fraction of susceptible depletion following
#  outbreaks. The "data" are replicate simulations from the MLEs, which
#  are also used for the bootsrapped CI analysis.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(pomp)
library(ggthemes)



# Define threshold for outbreak year --------------------------------------

# Outbreaks defined as: an annual case count exceeding x% of the maximum
# annual case count for a region, where x is the threshold set below.
outbreak_threshold <- 0.8  # 80% of max = outbreak year


# Load simulated data -----------------------------------------------------

cities <- c("Agadez", "Maradi", "Niamey", "Zinder")
sim_data <- tibble()
for(do_city in cities){
  tmp_file <- paste0("../simulations/bootstrap-sims-", do_city, ".RDS")
  tmp_data <- readRDS(tmp_file) %>%
    unnest() %>%
    mutate(city = do_city)
  sim_data <- bind_rows(sim_data, tmp_data)
}


# Identify outbreak years and calculate S depletion -----------------------

outbreak_years <- sim_data %>%
  mutate(year = trunc(time)) %>%
  group_by(city, sim, year) %>%
  summarise(total_cases = sum(reports),
            max_S = max(S),
            min_S = min(S)) %>%
  mutate(outbreak = ifelse(total_cases > outbreak_threshold*max(total_cases), 
                           TRUE, FALSE),
         max_S_t_minus_1 = lag(max_S),
         susc_depletion = min_S / max_S_t_minus_1) %>%
  filter(outbreak == TRUE) %>%
  drop_na()


# Calculate fraction of depletions less than 0.5 --------------------------

fracs_less_than_half <- outbreak_years %>%
  ungroup() %>%
  group_by(city) %>%
  mutate(
    deplete_half = ifelse(susc_depletion < 0.5, TRUE, FALSE)
  ) %>%
  group_by(city, deplete_half) %>%
  count() %>%
  group_by(city) %>%
  mutate(
    frac_less = n/sum(n)
  ) %>%
  filter(deplete_half == TRUE) %>%
  dplyr::select(city, frac_less)

ggplot(fracs_less_than_half, aes(x = city, y = frac_less*100, color = city)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = city), yend=0) +
  expand_limits(y=0) +
  labs(x = "City", y = "Percentage of outbreaks with\nsusceptible depletion less than 0.5") +
  theme_classic(base_size = 12) +
  scale_color_colorblind() +
  guides(color = FALSE) +
  coord_flip()

# Plot boxplots of susceptible depletion ----------------------------------

ggplot(data = outbreak_years, aes(x = city, y = susc_depletion)) +
  geom_boxplot(aes(fill = city, color = city), alpha = 0.4, width = 0.5) +
  geom_hline(aes(yintercept = 0.5), color = "grey") +
  labs(x = "City", y = "Level of susceptible depletion") +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_classic(base_size = 12) +
  coord_flip() +
  guides(color = FALSE, fill = FALSE)

ggplot() +
  geom_histogram(data = filter(outbreak_years, city == "Agadez"), 
                 aes(x = susc_depletion, fill = city), 
                 color = "white", binwidth = 0.05) +
  geom_histogram(data = filter(outbreak_years, city == "Maradi"), 
                 aes(x = susc_depletion, fill = city), 
                 color = "white", binwidth = 0.05) +
  geom_histogram(data = filter(outbreak_years, city == "Niamey"), 
                 aes(x = susc_depletion, fill = city), 
                 color = "white", binwidth = 0.05) +
  geom_histogram(data = filter(outbreak_years, city == "Zinder"), 
                 aes(x = susc_depletion, fill = city), 
                 color = "white", binwidth = 0.05) +
  geom_vline(aes(xintercept = 0.5), color = "grey25", linetype = 2) +
  labs(y = "Number of outbreaks", x = "Level of susceptible depletion") +
  scale_fill_colorblind(name = NULL) +
  theme_classic(base_size = 12)

ggplot(data = outbreak_years, aes(x = susc_depletion)) +
  geom_freqpoly(aes(color = city), binwidth = 0.05) +
  labs(y = "Number of outbreaks", x = "Level of susceptible depletion") +
  scale_color_colorblind() +
  theme_classic(base_size = 12)
  