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

outplot <- plot_grid(the_map, the_series, ncol = 2, labels = "AUTO", label_size = 12)

ggsave(filename = "../figures/map-and-series.pdf", plot = outplot, height = 4, width = 8.5, units = "in")


# Calculate model-data agreement ------------------------------------------

for(do_city in unique(pred_cases$city_name)){
  tmp_data <- pred_cases %>%
    filter(city_name == do_city)
  
  mod <- lm(mean_cases ~ observed_cases, data = tmp_data)
}

# r_squares <- pred_cases %>%
#   group_by(city_name) %>%
#   summarise(
#     rmse = sqrt(mean((mean_cases - observed_cases)^2)),
#     corr = cor(mean_cases, observed_cases)
#   )
  
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
  geom_point(aes(color = city_name), size = 1, alpha = 0.3) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_label(data = r_squares, 
             aes(x = 1.7, y = 7.3, label = paste0("italic(R)^2 == ", 
                                                  round(R2,2))), 
             label.size = NA, parse = TRUE, size = 3) +
  labs(y = "Expected log(cases + 1)", x = "Observed log(cases + 1)") +
  facet_wrap(~city_name, nrow = 2, scales = "free") +
  scale_y_continuous(limits = c(0,8)) +
  scale_x_continuous(limits = c(0,8)) +
  scale_color_colorblind()+
  guides(color = FALSE)+
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  ggtitle("A. Model-data agreement")
# ggsave(filename = "../figures/pred-obs-scatters.pdf", plot = scatters, width = 8.5, height = 2.5, units = "in")


# Plot R0 curves ----------------------------------------------------------

# Define function to calculate R0 from seasonal params
calc_R0 <- function(beta = NULL, B = NULL, qis = NULL, season = NULL, 
                    eta = (365/8), mu = 0.05, nu = 0.05, gamma = (365/5)){
  if(is.null(beta)){
    R0 <- (eta*B*mu) / (nu*(eta + nu)*(gamma + nu))
  } else{
    B <- as.numeric((1 + exp(season %*% qis)) * beta)
    R0 <- (eta*B*mu) / (nu*(eta + nu)*(gamma + nu))
  }
  return(R0)
}

# Define computation grid for bootstraps

nboots <- 100
nmifs <- 50
comp_grid <- expand.grid(1:nboots, 1:50)
colnames(comp_grid) <- c("boot_series", "param_set")
comp_grid$do_grid <- 1:nrow(comp_grid)

# Load example pomp file for basis function
pomp_file <- "../code/measles-pomp-object-Agadez.RDS"
measles_pomp <- readRDS(pomp_file)  # exemplar bases

bases <- as_tibble(measles_pomp@covar) %>%
  dplyr::select(starts_with("x")) %>%
  dplyr::slice(1:365) %>%
  mutate(
    day = 1:365
  ) %>%
  gather(key = base, value = value, -day)

season <- bases %>%
  spread(key = base, value = value) %>%
  dplyr::select(-day) %>%
  as.matrix()

seasonal_functions <- tibble()
for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  boots <- read.csv(paste0("../results/bootstrap-mif-lls-", do_city, ".csv"))
  
  b_splines <- boots %>%
    slice(2:n()) %>%  # ignore first row of NAs
    left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
    group_by(boot_series) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    dplyr::select(b1, b2, b3, b4, b5, b6)
  
  betas <- boots %>%
    slice(2:n()) %>%  # ignore first row of NAs
    left_join(comp_grid, by = "do_grid") %>%  # merge in comp grid info
    group_by(boot_series) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    dplyr::select(beta_mu)
  
  seasonal_betas <- tibble()
  for(i in 1:nrow(b_splines)){
    qis <- as.numeric(b_splines[i, ])
    beta_tmp <- as.numeric(betas[i, "beta_mu"])
    
    seasonal_tmp <- tibble(
      beta = as.numeric((1+exp(season %*% qis)) * beta_tmp),
      day = 1:365,
      boot = i
    )
    seasonal_betas <- bind_rows(seasonal_betas, seasonal_tmp)
  }  # end bootstrap loop
  
  seasonal_betas <- seasonal_betas %>%
    mutate(
      city = do_city,
      rnaught = calc_R0(B = beta)
    ) %>%
    group_by(city, day) %>%
    summarise(
      upper_R0 = quantile(rnaught, 0.975),
      lower_R0 = quantile(rnaught, 0.025)
    ) %>%
    ungroup()
  
  seasonal_functions <- bind_rows(seasonal_functions, seasonal_betas)
  
}  # end city loop

seasonal_functions <- seasonal_functions %>%
  mutate(
    date = as.Date(day, origin = "2016-12-31",tz = "UTC")
  )

# Calculate MLE R0 curves for each city

R_0s <- tibble()

for(do_city in c("Agadez", "Maradi", "Niamey", "Zinder")){
  
  # Load fitted parameters and pomp model 
  mle_file <- paste0("../results/initial-mif-lls-", do_city, ".csv")
  
  mles <- read.csv(mle_file) %>% 
    slice(2:n()) %>%  # ignore first row of storage NAs
    filter(loglik == max(loglik, na.rm = TRUE)) %>%
    dplyr::select(-do_grid, -loglik, -loglik_se)
  
  pomp_file <- paste0("./measles-pomp-object-", do_city, ".RDS")
  fitted_pomp <- readRDS(pomp_file)
  
  # Calculte R0 
  qis <- mles %>%
    dplyr::select(b1, b2, b3, b4, b5, b6) %>%
    as.numeric()
  
  beta <- mles %>%
    pull(beta_mu)
  
  bases <- as_tibble(fitted_pomp@covar) %>%
    dplyr::select(starts_with("x")) %>%
    dplyr::slice(1:365) %>%
    mutate(
      day = 1:365
    ) %>%
    gather(key = base, value = value, -day)
  
  season <- bases %>%
    spread(key = base, value = value) %>%
    dplyr::select(-day) %>%
    as.matrix()
  
  N <- round(mean(fitted_pomp@covar[, "N"]))
  
  R0 <- calc_R0(beta = beta, qis = qis, season = season)
  
  tmp_out <- tibble(
    city = do_city,
    day = 1:length(R0),
    R0 = R0
  )
  
  R_0s <- bind_rows(R_0s, tmp_out)
}

R_0s <- R_0s %>%
  mutate(
    date = as.Date(day, origin = "2016-12-31",tz = "UTC")
  )

all_R_0s <- left_join(R_0s, seasonal_functions, by = c("city", "day", "date"))

rnaughts <- ggplot(all_R_0s, aes(x = date)) +
  geom_ribbon(aes(ymin = lower_R0, ymax = upper_R0, fill = city), alpha = 0.6) +
  geom_line(aes(y = R0), color = "white") +
  scale_fill_colorblind(name = NULL) +
  scale_color_colorblind(name = NULL) +
  guides(fill = FALSE) +
  facet_wrap(~city,nrow = 2, scales = "free") +
  labs(x = "Time of year", y = expression(R[0])) +
  scale_x_date(date_labels = "%b", date_breaks = "2 months") +
  scale_y_continuous(limits = c(0,35))+
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
       plot.title = element_text(face = "bold")) +
  ggtitle("B. Seasonal basic reproduction number")


# rnaughts <- ggplot(R_0s, aes(x = date, y = R0, color = city)) +
#   geom_line() +
#   scale_color_colorblind(name = NULL) +
#   scale_x_date(date_labels = "%b", date_breaks = "2 months") +
#   labs(x = "Time of year", y = expression(R[0])) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.9), 
#         legend.box.background = element_blank(),
#         legend.key.size = unit(0.8, 'lines'))

scatters_and_season <- cowplot::plot_grid(scatters, rnaughts, labels = NULL, scale = 0.9, label_size = 12)
ggsave(filename = "../figures/pred-obs-scatters.pdf", plot = scatters_and_season, width = 9.5, height = 5, units = "in")


