# plot-ews.R
#  Script to plot the time series of early warning signals for each city.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(lubridate)


# Load EWS and data -------------------------------------------------------

ews_file <- "../results/ews-niger-cities.RDS"
ews_data <- readRDS(ews_file) %>%
  filter(year(date) > 1996)  # remove first two years before a full window


# Plot the EWS time series ------------------------------------------------

##  First scale to 0-1 range so easy to see
scale_it <- function(x){
  (x-min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}

ews_data <- ews_data %>%
  group_by(city, ews) %>%
  mutate(
    scaled_value = scale_it(value)
  )

##  Now plot the scaled EWS values
ews_plot <- ggplot(ews_data, aes(x = date, y = scaled_value)) +
  geom_line(aes(color = city)) +
  facet_grid(city~ews) +
  labs(x = "Date", y = "EWS value (scaled to [0,1])") +
  scale_color_manual(values = ptol_pal()(4)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE)

##  Save the plot
ggsave(
  filename = "../figures/ews_ts_plot.pdf", 
  plot = ews_plot, 
  width = 14,
  height = 5, 
  units = "in"
)
