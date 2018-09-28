# plot-global-search.R
#  Script to plot the results from the initial global parameter search using
#  MIF. The 500 parameter sets with the highest likelihood are
#  saved for subsequent analysis.
#
# Author:
#  Andrew Tredennick

DO_CITY <- "Niamey"

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(pomp)

# Load data ---------------------------------------------------------------

mif_traces <- read.csv(paste0("../results/initial-mif-traces-", DO_CITY, ".csv")) %>%
  as_tibble() %>%
  slice(2:n())
mif_finals <- read.csv(paste0("../results/initial-mif-lls-", DO_CITY, ".csv")) %>%
  as_tibble() %>%
  slice(2:n())

mif_finals <- read.csv("~/Desktop/initial-mif-lls.csv") %>%
  as_tibble() %>%
  slice(2:n())

# Plot MIF traces ---------------------------------------------------------

best_grids <- mif_finals %>%
  filter(loglik < 0) %>%
  arrange(-loglik) %>%
  slice(1:200) %>%
  pull(do_grid)

mif_traces_long <- mif_traces %>%
  gather(key = parameter, value = value, -do_grid, -iteration) %>%
  filter(do_grid %in% best_grids)

ggplot(mif_traces_long, aes(x = iteration, y = value, group = do_grid)) +
  geom_line(alpha = 0.1) +
  facet_wrap(~parameter, scales = "free_y")

mifs_last <- mif_traces_long %>% filter(iteration == 50) %>%
  filter(parameter %in% c("loglik", "S_0")) %>%
  spread(parameter, value)

plot(mifs_last$loglik, mifs_last$S_0)

# Simulate model at MLEs --------------------------------------------------

mles <- mif_finals %>%
  filter(loglik < 0) %>%
  filter(loglik == max(loglik, na.rm = T)) %>%
  # filter(do_grid == 584) %>%
  dplyr::select(-do_grid, -loglik, -loglik_se)

measles_pomp <- readRDS(paste0("measles-pomp-object-", DO_CITY, ".RDS"))

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
  theme(strip.text=element_blank())

simulate(
  measles_pomp,
  params = unlist(mles),
  nsim = 1,
  as.data.frame = TRUE,
  include.data = FALSE) %>%
  ggplot(aes(x = time, y = RE, group = sim, color = (sim == "data"))) +
  geom_line() +
  geom_hline(aes(yintercept = 1), color = "red") +
  labs(y = expression(paste( R[E])))+
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "black"))+
  guides(color = FALSE) 


# Plot basis functions ----------------------------------------------------

bases <- as_tibble(measles_pomp@covar) %>%
  dplyr::select(starts_with("x")) %>%
  dplyr::slice(1:365) %>%
  mutate(
    day = 1:365
  ) %>%
  gather(key = base, value = value, -day)

basis_plot <- ggplot(bases, aes(x = day, y = value, color = base)) +
  geom_line() +
  scale_color_manual(values = ggthemes::ptol_pal()(6)) +
  guides(color = FALSE) +
  labs(y = "Basis value", x = "Day") +
  theme_minimal()

season <- bases %>%
  spread(key = base, value = value) %>%
  dplyr::select(-day) %>%
  as.matrix()

qis <- mles %>%
  dplyr::select(starts_with("b")) %>%
  dplyr::select(-beta_mu, -beta_sd) %>%
  t()

seasonal_beta <- tibble(
  beta = as.numeric((1+exp(season %*% qis)) * mles$beta_mu),
  day = 1:365
)

beta_plot <- ggplot(seasonal_beta, aes(x = day, y = beta)) +
  geom_line() +
  labs(x = "Day", y = expression(paste("Tranmission rate, ",beta, " (", yr^-1, ")"))) +
  # scale_y_continuous(breaks = seq(5,30,5)) +
  theme_minimal()

cowplot::plot_grid(basis_plot, beta_plot, ncol = 1, align = "v")

# Plot likelihood surface -------------------------------------------------

final_mles <- mif_finals %>%
  filter(loglik < 0) %>%
  dplyr::select(-do_grid, -loglik_se)

pairs(final_mles, pch = ".")

ggplot(final_mles, aes(x = S_0, y = rho, color = loglik)) +
  geom_point()
