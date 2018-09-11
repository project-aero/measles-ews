# plot-global-search.R
#  Script to plot the results from the initial global parameter search using
#  MIF. The 500 parameter sets with the highest likelihood are
#  saved for subsequent analysis.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)


# Load data ---------------------------------------------------------------

mif_traces <- read_csv("../results/initial-mif-traces.csv")
mif_finals <- read_csv("../results/initial-mif-lls.csv")


# Plot MIF traces ---------------------------------------------------------

best_grids <- mif_finals %>%
  arrange(-loglik)

best_grids <- best_grids[1:200, ] %>%
  pull(do_grid)

mif_traces_long <- mif_traces %>%
  gather(key = parameter, value = value, -do_grid, -iteration) %>%
  filter(do_grid %in% best_grids)

ggplot(mif_traces_long, aes(x = iteration, y = value, group = do_grid)) +
  geom_line(alpha = 0.1) +
  facet_wrap(~parameter, scales = "free_y")


# Simulate model at MLEs --------------------------------------------------

mles <- mif_finals %>%
  filter(loglik == max(loglik)) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se)

measles_pomp <- readRDS("measles-pomp-object.RDS")

simulate(
  measles_pomp,
  params = unlist(mles),
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = reports, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())

simulate(
  measles_pomp,
  params = unlist(mles),
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = S, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())
