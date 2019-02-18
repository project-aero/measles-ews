# plot-data.R
#  R script to plot the weekly case reports from four cities in Niger.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggthemes)


# Load data ---------------------------------------------------------------

file_name <- "../data/raw-data/niger_regional_1995_2005.csv"
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
  filter(grepl("City", region)) %>%
  separate(region, " ", into = c("region", "toss")) %>%
  dplyr::select(-toss)


# Make the plot and save --------------------------------------------------

ggplot(measles_data, aes(x = date, y = cases, color = region)) +
  geom_path() +
  facet_wrap(~region, scales = "free_y", nrow = 1) +
  scale_color_manual(values = ptol_pal()(4)) +
  guides(color = FALSE) +
  labs(x = "Date", y = "Reported measles cases") +
  theme_minimal()

ggsave(
  filename = "../figures/data-timeseries.pdf", 
  height = 3, 
  width = 8.5, 
  units = "in"
)
