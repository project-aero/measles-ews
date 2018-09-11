# fetch-clean-data.R:
#  Script to load and clean up raw data, which are weekly incidence of measles 
#  in 40 regions of Niger. This script adds some date information and then
#  subsets out just the major cities (four of the 40 regions). The data are
#  saved as a tibble in .RDS format.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)


# Read in data and format -------------------------------------------------

file_name <- "../data/raw-data/niger_regional_1995_2005.csv"  # fromt the AERO data repo
niger_measles_raw <- suppressWarnings(
  read_csv(file_name, col_types = cols())
) 


# Generate some information based on the data
num_regions <- nrow(niger_measles_raw)  # one region per row (columns are weeks)
num_weeks <- ncol(niger_measles_raw) - 1  # subtract 1 from ncol() because first column are regions
weeks_per_year <- num_weeks/11  # we know where are 11 years
weeks <- rep(1:52, times = 11)  # data preformatted to all years having 52 weeks

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
    obs_week = rep(1:num_weeks, each = num_regions),
    year = rep(years, each  = num_regions),
    week_of_year = rep(weeks, each = num_regions),
    date = calculate_start_of_week(week_of_year, year)
  ) %>%
  dplyr::rename(region = X1) %>%
  filter(grepl("City", region)) %>%
  dplyr::select(region, date, year, week_of_year, obs_week, cases) %>%
  arrange(region, date) %>%
  mutate(
    time = decimal_date(date)
  )


# Read in demographic data and format ------------------------------------

birth_file <- "../data/raw-data/niger_crude_birth_rates.csv"
birth_data <- read_csv(birth_file, col_types = cols()) %>%
  mutate(
    date = mdy(date),  # lubridate prefixes any 2digit year 00-68 with 20, not a problem for us though
    year = as.character(year(date)),
    rate_per_person = births_per_thousand/1000
  ) %>%
  dplyr::select(year, rate_per_person) %>%
  mutate(year = as.numeric(year)) %>%
  filter(year > 1990 & year < 2010)

pop_file <- "../data/raw-data/district_pops.csv"
city_strings <- str_sub(unique(measles_data$region), start = 1, end = 6)
pop_data <- suppressWarnings(
  read_csv(pop_file, col_types = cols())
  )  %>%
  gather(key = year, value = population, -X1) %>%
  rename(region = X1) %>%
  mutate(
    region = ifelse(region == "Niamey I", "Niamey", region),
    year = as.numeric(year)
  ) %>%
  filter(region %in% city_strings)

demog_data <- pop_data %>%
  left_join(birth_data, by = "year") %>%
  mutate(
    per_capita_birth_rate = rate_per_person,
    births_per_year = population*per_capita_birth_rate,
    births_per_week = births_per_year/52
  ) %>%
  dplyr::select(-rate_per_person) %>%
  arrange(region, year) %>%
  group_by(region) %>%
  mutate(
    obs_week = rep(seq(from = 1, to = max(measles_data$obs_week), by = 52))
  )


# Fit splines to smooth annual data to weekly -----------------------------

smoothed_demog <- tibble()
all_weeks <- unique(measles_data$obs_week)

for(do_city in unique(demog_data$region)){
  tmp <- demog_data %>%
    filter(region == do_city)
  
  b <- predict(
    smooth.spline(x = tmp$obs_week, y = tmp$births_per_week), 
    x = all_weeks)$y
  
  N <- predict(
    smooth.spline(x = tmp$obs_week, y = tmp$population), 
    x = all_weeks)$y
  
  out <- tibble(
    region = paste(do_city, "(City)"),  # to match up with case dataframe
    obs_week = all_weeks,
    population_smooth = N,
    births_per_week_smooth = b
  )
  
  smoothed_demog <- bind_rows(
    smoothed_demog,
    out
  )
}
  


# Save the cleaned data ---------------------------------------------------

saveRDS(
  object = demog_data,
  file = "../data/clean-data/annual-demographic-data-niger-cities-clean.RDS"
)

measles_full <- measles_data %>%
  left_join(smoothed_demog,  by = c("region", "obs_week"))

saveRDS(
  object = measles_full,
  file = "../data/clean-data/weekly-measles-and-demog-niger-cities-clean.RDS"
)

