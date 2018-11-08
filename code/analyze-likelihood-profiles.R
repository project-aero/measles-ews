# analyze-likelihood-profiles.R

mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000)
{
  smooth_fit <- loess(lp ~ parameter, family = "s", span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2))
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  var_b <- m["b","b"]
  var_a <- m["a","a"]
  cov_ab <- m["a","b"]
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
  delta_line <- max(smoothed_loglik) - delta
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp,parameter=parameter,confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(
         parameter=parameter_grid,
         smoothed=smoothed_loglik,
         quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
       ),
       mle=smooth_arg_max, ci=ci, delta=delta, delta_line=delta_line,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
  )
}


# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggthemes)



# Niamey, beta profile ----------------------------------------------------
profile_data <- read.csv("../results/loglik-profile-beta-Niamey.csv") %>%
  slice(2:n()) %>%
  drop_na() %>%
  group_by(value, parameter) %>%
  summarise(loglik = pomp::logmeanexp(loglik))

do_param <- unique(profile_data$parameter)
mcap_out <- mcap(lp = profile_data$loglik, parameter = profile_data$value, lambda = 0.5)

ggplot(mcap_out$fit, aes(x=parameter, shape )) +
  geom_point(data = profile_data, aes(x = value, y = loglik), shape = 19, size = 2, color = "grey50", alpha = 0.1) +
  geom_line(aes(y = smoothed), color = ptol_pal()(2)[2], size = 1) +
  geom_line(aes(y = quadratic), color = ptol_pal()(2)[1], linetype = 2, size = 1) +
  geom_vline(aes(xintercept = mcap_out$ci[1]), color = ptol_pal()(2)[2]) +
  geom_vline(aes(xintercept = mcap_out$ci[2]), color = ptol_pal()(2)[2]) +
  geom_hline(aes(yintercept = mcap_out$delta_line), color = ptol_pal()(2)[2]) +
  labs(x = expression(beta), y = "profile log-likelihood") +
  # coord_cartesian(xlim = c(50, 1000), ylim = c(-1500, -1450)) +
  ggtitle(paste0("95% CI: ", round(mcap_out$ci,2)[1], " - ", round(mcap_out$ci,2)[2]))


# Niamey, rho profile -----------------------------------------------------

profile_data <- read.csv("../results/loglik-profile-rho-Niamey.csv") %>%
  slice(2:n()) %>%
  # drop_na() %>%
  group_by(value, parameter) %>%
  summarise(loglik = pomp::logmeanexp(loglik)) %>%
  ungroup() %>%
  mutate(lik_diff = max(loglik) - loglik) %>%
  filter(value > 0.01)
  # filter(lik_diff < 100)

do_param <- unique(profile_data$parameter)
mcap_out <- mcap(lp = profile_data$loglik, parameter = profile_data$value, lambda = 0.5)

plot(profile_data$value, profile_data$loglik)
abline(v = 0.26)
smooth_fit <- loess(formula = loglik~value, data = profile_data, span = 0.2)
parameter_grid <- seq(min(profile_data$value), max(profile_data$value), length.out = 1000)
smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
delta <- qchisq(0.99, df = 1)
loglik_diff <- max(smoothed_loglik) - smoothed_loglik
ci <- range(parameter_grid[loglik_diff < delta])
lines(parameter_grid, smoothed_loglik, col = "red")


ggplot(mcap_out$fit, aes(x=parameter)) +
  geom_point(data = profile_data, aes(x = value, y = loglik), shape = 1, size = 2, color = "grey50") +
  geom_line(aes(y = smoothed), color = "coral", size = 1) +
  # geom_line(aes(y = quadratic), color = "steelblue", linetype = 2, size = 1) +
  geom_vline(aes(xintercept = mcap_out$ci[1]), color = "coral") +
  geom_vline(aes(xintercept = mcap_out$ci[2]), color = "coral") +
  geom_hline(aes(yintercept = mcap_out$delta_line), color = "coral") +
  labs(x = expression(rho), y = "profile log-likelihood") +
  coord_cartesian(ylim = c(-4000, -1450)) +
  ggtitle(paste0("95% CI: ", round(mcap_out$ci,2)[1], " - ", round(mcap_out$ci,2)[2]))


# Niamey, iota profile ----------------------------------------------------

profile_data <- read.csv("../results/loglik-profile-iota-Niamey.csv") %>%
  slice(2:n()) %>%
  drop_na()

do_param <- unique(profile_data$parameter)
mcap_out <- mcap(lp = profile_data$loglik, parameter = profile_data$value)

ggplot(mcap_out$fit, aes(x=parameter)) +
  geom_point(data = profile_data, aes(x = value, y = loglik), shape = 1, size = 2, color = "grey50") +
  geom_line(aes(y = smoothed), color = "coral", size = 1) +
  geom_line(aes(y = quadratic), color = "steelblue", linetype = 2, size = 1) +
  geom_vline(aes(xintercept = mcap_out$ci[1]), color = "coral") +
  geom_vline(aes(xintercept = mcap_out$ci[2]), color = "coral") +
  geom_hline(aes(yintercept = mcap_out$delta_line), color = "coral") +
  labs(x = expression(rho), y = "profile log-likelihood") +
  # coord_cartesian(ylim = c(-5000, -1450)) +
  ggtitle(paste0("95% CI: ", round(mcap_out$ci,2)[1], " - ", round(mcap_out$ci,2)[2]))




