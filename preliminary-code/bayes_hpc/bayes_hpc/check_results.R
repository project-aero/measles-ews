# Check MCMC results

library(tidyverse)
library(lubridate)

# Read in data ------------------------------------------------------------

file_name <- "niger_regional_1995_2005.csv"
niger_measles_raw <- read_csv(file_name, col_types = cols())

num_regions <- nrow(niger_measles_raw)
num_weeks <- ncol(niger_measles_raw) - 1  # subtract 1 from ncol() because first column are regions
weeks_per_year <- num_weeks/11
weeks <- rep(1:52, times = 11)

# Create a vector of years for all num_weeks
years <- rep(1995:2005, each = weeks_per_year)

# Function for calculating start of week based on week number and year
calculate_start_of_week = function(week, year) {
  date <- ymd(paste(year, 1, 1, sep="-"))
  week(date) = week
  return(date)
}

# Clean up the data frame
measles_data <- niger_measles_raw %>%
  gather(key = week, value = cases, -X1) %>%
  mutate(
    week_num = rep(1:num_weeks, each = num_regions),
    year = rep(years, each  = num_regions),
    week = rep(weeks, each = num_regions),
    date = calculate_start_of_week(week, year)
  ) %>%
  dplyr::rename(region = X1) %>%
  filter(region == "Niamey (City)")

biweek_data <- measles_data %>%
  mutate(biweek = rep(1:(n()/2), each = 2)) %>%
  group_by(biweek) %>%
  summarise(
    cases = sum(cases),
    date = min(date)
  )


summaries <- read_csv("mcmc_results.csv")

cases <- summaries %>%
  filter(grepl("Iobs\\[", Parameter)) %>%
  mutate(
    biweek = 1:n(),
    observations = pull(biweek_data, cases)
  )

ggplot(cases, aes(x = biweek, y = median_value)) +
  geom_point(aes(y = observations)) +
  geom_ribbon(aes(ymax = upper_95, ymin = lower_95), alpha = 0.2) +
  geom_line() +
  theme_minimal()

latents <- summaries %>%
  filter(grepl("S\\[", Parameter)) %>%
  mutate(
    biweek = 1:n()
  )

plot(latents$median_value, type = "l")


chains <- readRDS("mcmc_chains.RDS")
