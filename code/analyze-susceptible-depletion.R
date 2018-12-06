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
  mutate(outbreak = ifelse(total_cases > outbreak_threshold*max(total_cases), TRUE, FALSE),
         max_S_t_minus_1 = lag(max_S),
         susc_depletion = min_S / max_S_t_minus_1) %>%
  filter(outbreak == TRUE) %>%
  drop_na()


# Plot boxplots of susceptible depletion ----------------------------------

ggplot(data = outbreak_years, aes(x = city, y = susc_depletion)) +
  geom_boxplot(fill = "tan", color = "tan", alpha = 0.4, width = 0.5) +
  labs(x = "City", y = "Level of susceptible depletion") +
  theme_classic(base_size = 12) +
  coord_flip()
