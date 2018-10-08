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


# Load filtered states ----------------------------------------------------

filter_ids <- grep("filtered-states", list.files("../results/"))
filter_files <- list.files("../results/")[filter_ids]

filtered_states <- tibble()
for(do_file in filter_files){
  do_city <- str_split(do_file, "-")[[1]][3]
  do_city <- str_split(do_city, "[.]")[[1]][1]
  
  tmp <- readRDS(paste0("../results/", do_file)) %>%
    filter(state == "cases") %>%
    unnest() %>%
    mutate(
      city_name = do_city
    )
  filtered_states <- bind_rows(filtered_states, tmp)
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

the_series <- ggplot(filtered_states, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), color = "black") +
  geom_line(aes(y = observation), color = ptol_pal()(3)[2], size = 0.4) +
  facet_wrap(~city_name, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Reported cases") +
  theme_minimal()

outplot <- plot_grid(the_map, the_series, ncol = 2)

ggsave(filename = "../figures/map-and-series.pdf", plot = outplot, height = 3, width = 8.5, units = "in")
