# plot-data-and-filters.R:
#  Script to plot the spatial locations of each focal city and the observed
#  incidence time series. The time series plot also shows the estimated 
#  case reports from the fitted model via particle filtering.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(ggrepel)
library(cowplot)
library(pomp)
library(sf)


# Load spatial data -------------------------------------------------------

region_boundaries <- st_read("../data/spatial-data/NER_adm2.shp")

mycoords <- as_tibble(st_coordinates(region_boundaries)) %>%
  group_by(L3) %>%
  summarise(minx = min(X)) %>%
  mutate(
    NAME_2 = region_boundaries$NAME_2
  )

region_boundaries <- region_boundaries %>%
  left_join(mycoords)

my_cities <- tibble(
  city = c("Agadez", "Maradi", "Niamey", "Zinder"),
  lat = c(16.9, 13.5, 13.5, 13.8),
  lon = c(8.0, 7.1, 2.1, 8.98),
  population = c("(118,244)", "(267,249)", "(1,027,000)", "(322,935)")
) %>%
  mutate(
    label = paste(city, population, sep = "\n")
  )


# Load case observation data ----------------------------------------------

file_name <- "../data/clean-data/weekly-measles-incidence-niger-cities-clean.RDS"
measles_data <- readRDS(file_name) %>%
  filter(year > 1994) %>%
  dplyr::select(-year, -week_of_year, -obs_week, -time) %>%
  mutate(
    region = str_sub(region, 1, str_length(region)-7)  # remove (City) from name
  )



# Load predictive distributions -------------------------------------------

pred_ids <- grep("predictive-dist-states", list.files("../results/"))
pred_files <- list.files("../results/")[pred_ids]

pred_cases <- tibble()
for(do_file in pred_files){
  do_city <- str_split(do_file, "-")[[1]][4]
  do_city <- str_split(do_city, "[.]")[[1]][1]
  
  tmp_data <- measles_data %>%
    filter(region == do_city)
  
  tmp <- readRDS(paste0("../results/", do_file)) %>%
    unnest() %>%
    slice(2:n()) %>%  # drop first row of unobserved data
    mutate(
      upper_est = mean_cases + sdev_cases*2,
      lower_est = mean_cases - sdev_cases*2,
      lower_est = ifelse(lower_est < 0, 0, lower_est),
      city_name = do_city,
      date = tmp_data$date,
      observed_cases = tmp_data$cases
    )
  pred_cases <- bind_rows(pred_cases, tmp)
}


# Make the plots ----------------------------------------------------------

the_map <- ggplot(region_boundaries)+
  geom_sf(fill = "grey90", col = "grey80", size = 0.3) +
  geom_point(data = my_cities, aes(x = lon, y = lat), size = 2) +
  geom_label_repel(
    data = my_cities,
    aes(x = lon, y = lat, label = label),
    box.padding   = 0.35, 
    point.padding = 0.3,
    size = 2
  ) + 
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()

the_series <- ggplot(pred_cases, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_est, ymax = upper_est), color = "black") +
  geom_line(aes(y = observed_cases), color = ptol_pal()(3)[2], size = 0.4) +
  facet_wrap(~city_name, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Reported cases") +
  theme_minimal()

outplot <- plot_grid(the_map, the_series, ncol = 2, labels = "AUTO")

ggsave(filename = "../figures/map-and-series.pdf", plot = outplot, height = 4, width = 8.5, units = "in")


# Calculate model-data agreement ------------------------------------------

r_squares <- pred_cases %>%
  group_by(city_name) %>%
  mutate(
    mean_observations = mean(observed_cases),
    error_numer = (mean_cases - observed_cases)^2,
    error_denom = (mean_observations - observed_cases)^2
  ) %>%
  summarise(
    summed_numer = sum(error_numer),
    summed_denom = sum(error_denom)
  ) %>%
  mutate(
    R2 = 1 - (summed_numer / summed_denom)
  ) %>%
  dplyr::select(city_name, R2)


scatters <- ggplot(pred_cases, aes(x = log(observed_cases+1), y = log(mean_cases+1))) +
  geom_point(color = ptol_pal()(1), size = 2, alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2, size = 1) +
  geom_label(data = r_squares, aes(x = 1.7, y = 7.3, label = paste0("italic(R)^2 == ", round(R2,2))), label.size = NA, parse = TRUE) +
  labs(y = "Expected log(cases + 1)", x = "Observed log(cases + 1)") +
  facet_wrap(~city_name, nrow = 1) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"))
ggsave(filename = "../figures/pred-obs-scatters.pdf", plot = scatters, width = 8.5, height = 2.5, units = "in")
