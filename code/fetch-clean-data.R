# fetch-clean-data.R:
#  Script to load and clean up raw data, which are weekly incidence of measles 
#  in 40 regions of Niger. This script reads in the raw case data and reformats
#  it to include date information. Annual birth rate (national) and total
#  population size (regional) are also read in, interpolated to be match
#  the weekly time interval of case reports and added to the data frame.
#  Note that we only save output for our four focal cities: Agadez, Maradi,
#  Niamey, and Zinder.
#
#  Two outputs are produced:
#    demog_data: Birth rates and population size for each city at each week.
#      file: ../data/clean-data/annual-demographic-data-niger-cities-clean.RDS
#    measles_data: Timestamped case reports for each week, by city.
#      file: ../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)


# Read in data and format -------------------------------------------------

# From the AERO data repo
file_name <- "../data/raw-data/niger_regional_1995_2005.csv"

niger_measles_raw <- read_csv(
  file_name, 
  col_types = cols(X1 = col_character(), .default = col_integer()),
  col_names = FALSE)


# Generate some information based on the data
num_regions <- nrow(niger_measles_raw)  # one region per row (columns are weeks)
num_weeks <- ncol(niger_measles_raw) - 1  # subtract 1 from ncol() because 
                                          # first column are regions
weeks_per_year <- num_weeks/11  # we know there are 11 years
weeks <- rep(1:52, times = 11)  # data preformatted to all years having 52 weeks

# Create a vector of years for all num_weeks
years <- rep(1995:2005, each = weeks_per_year)

# Function for calculating start of week based on week number and year
calculate_start_of_week = function(week, year){
  date <- ymd(paste(year, 1, 1, sep="-"))
  week(date) <- week
  return(date)
}

# Clean up the data frame
measles_data <- niger_measles_raw %>%
  gather(key = week, value = cases, -X1) %>%
  mutate(
    obs_week = rep(1:num_weeks, each = num_regions),
    obs_week2 = stringr::str_replace(week, "^X", ""),
    obs_week2 = as.integer(obs_week2) - 1,
    year = rep(years, each  = num_regions),
    year2 = (obs_week2 - 1) %/% 52 + 1995,
    week_of_year = rep(weeks, each = num_regions),
    week_of_year2 = (obs_week2 - 1) %% 52 + 1,
    date = calculate_start_of_week(week_of_year, year)
  ) %>%
  dplyr::rename(region = X1) %>%
  filter(grepl("City", region)) %>%
  dplyr::select(region, date, year, year2, week_of_year, week_of_year2, 
                obs_week, obs_week2, cases) %>%
  arrange(region, date) %>%
  mutate(
    time = decimal_date(date)
  )

stopifnot(all(measles_data$year == measles_data$year2))
stopifnot(all(measles_data$obs_week == measles_data$obs_week2))
stopifnot(all(measles_data$week_of_year == measles_data$week_of_year2))
measles_data[c("year2", "obs_week2", "week_of_year2")] <- NULL

# Make tibble with initial conditions NA row, 1 week before data start
initial_conditions_datarow <- tibble(
  region = unique(measles_data$region),
  date = min(measles_data$date) - 7,
  year = year(min(measles_data$date) - 7),
  week_of_year = 52,
  obs_week = NA,
  cases = NA
) %>%
  mutate(
    time = decimal_date(date)
  )

# Add in initial time NAs for observations
measles_data <- measles_data %>%
  bind_rows(initial_conditions_datarow) %>%
  arrange(region, date)


# Read in demographic data and format ------------------------------------

# Define the time variable for 365 day years
time_tbl <- tibble(
  date = seq(min(measles_data$date), max(measles_data$date), by = "day"),
  time = decimal_date(date),
  year = year(date)
) %>%
  dplyr::select(-date)

# Births
birth_file <- "../data/raw-data/niger_crude_birth_rates.csv"
birth_data <- read_csv(birth_file, col_types = cols()) %>%
  mutate(
    date = mdy(date),  # lubridate prefixes any 2digit year 00-68 with 20, 
                       # not a problem for us though
    year = as.character(year(date)),
    rate_per_person_per_year = births_per_thousand/1000
  ) %>%
  dplyr::select(year, rate_per_person_per_year) %>%
  mutate(year = as.numeric(year)) %>%
  filter(year > 1990 & year < 2010) %>%
  mutate(time = year+0.000)

# Interpolate birth rates
birth_spline <- predict(
  smooth.spline(x = birth_data$time, y = birth_data$rate_per_person_per_year), 
  x = time_tbl$time)$y

birth_rates <- time_tbl %>%
  mutate(
    birth_per_person_per_year = birth_spline
  )

# Population size
pop_file <- "../data/raw-data/district_pops.csv"
city_strings <- str_sub(unique(measles_data$region), start = 1, end = 6)
pop_data <- suppressWarnings(
  read_csv(pop_file, col_types = cols())
  )  %>%
  gather(key = year, value = population, -X1) %>%
  rename(region = X1) %>%
  mutate(
    region = ifelse(region == "Niamey I", "Niamey", region),
    year = as.numeric(year),
    time = year + 0.000  # make time double
  ) %>%
  filter(region %in% city_strings)

population_sizes <- tibble()
all_times <- unique(time_tbl$time)

for(do_city in unique(pop_data$region)){
  tmp <- pop_data %>%
    filter(region == do_city)
  
  # Interpolate population size
  N <- predict(
    smooth.spline(x = tmp$time, y = tmp$population), 
    x = all_times)$y
  
  out <- tibble(
    region = paste(do_city, "(City)"),  # to match up with case dataframe
    time = all_times,
    population_size = N
  )
  
  population_sizes <- bind_rows(
    population_sizes,
    out
  )
}
  


# Save the cleaned data ---------------------------------------------------

demog_data <- population_sizes %>%
  left_join(birth_rates, by = "time") %>%
  dplyr::select(-year)


if(!dir.exists("../data/clean-data")){
  dir.create("../data/clean-data")
}

saveRDS(
  object = demog_data,
  file = "../data/clean-data/annual-demographic-data-niger-cities-clean.RDS"
)

saveRDS(
  object = measles_data,
  file = "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
)

