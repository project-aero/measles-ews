# analyze-likelihood-profiles.R


# Load packages -----------------------------------------------------------

library(tidyverse)


# Load profile results ----------------------------------------------------

profile_data <- read.csv("../results/rho-profile-Niamey.csv") %>%
  slice(2:n()) %>%  # chop off first NA row
  filter(loglik > max(loglik) - 20)

preds <- predict(
  mgcv::gam(
    loglik ~ s(rho_value, bs = "cs"), data = profile_data
  ),
  newdata = data.frame(
    rho_value = seq(0.2, 0.54, by = 0.0001)
  )
)

pred_df <- tibble(
  rho_value = seq(0.2, 0.54, by = 0.0001),
  loglik = preds
)

vspots <- pred_df %>%
  filter(loglik > max(loglik, na.rm = T) - 1.92) %>%
  pull(rho_value)

vspots2 <- profile_data %>%
  filter(loglik > max(pred_df$loglik, na.rm = T) - 1.92) %>%
  pull(rho_value)

cspot <- pred_df %>%
  filter(loglik == max(loglik, na.rm = T)) %>%
  mutate(
    loglik - 1.92
  ) %>%
  pull(loglik)


ggplot(profile_data, aes(x = rho_value, y = loglik)) +
  geom_point(shape = 1) +
  geom_hline(aes(yintercept = max(preds) - 1.92), color = "blue") +
  geom_line(data = pred_df, aes(x = rho_value, y = loglik), color = "red")+
  geom_vline(aes(xintercept = min(vspots)), color = "blue") +
  geom_vline(aes(xintercept = max(vspots)), color = "blue") +
  geom_vline(aes(xintercept = min(vspots2)), color = "blue", linetype = 2) +
  geom_vline(aes(xintercept = max(vspots2)), color = "blue", linetype = 2) +
  labs(x = expression(rho), y = "profile log-likelihood") +
  theme_minimal()
