# analyze-likelihood-profiles.R

mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000)
{
  smooth_fit <- loess(lp ~ parameter,span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter-smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0,length(parameter))
  weight[included] <- (1-(dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp,b=parameter,a=-parameter^2)
  )
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


# Load profile results ----------------------------------------------------

profile_data <- read.csv("../results/rho-profile-Niamey.csv") %>%
  slice(2:n()) %>%  # chop off first NA row
  sample_n(100)

test <- mcap(lp = profile_data$loglik, parameter = profile_data$rho_value)

ggplot(test$fit, aes(x=parameter)) +
  geom_point(data = profile_data, aes(x = rho_value, y = loglik), shape = 1, color = "black") +
  geom_line(aes(y = smoothed), color = "red") +
  geom_line(aes(y = quadratic), color = "blue", linetype = 2) +
  geom_vline(aes(xintercept = test$ci[1]), color = "red") +
  geom_vline(aes(xintercept = test$ci[2]), color = "red") +
  geom_hline(aes(yintercept = test$delta_line), color = "red") +
  labs(x = expression(rho), y = "profile log-likelihood") +
  theme_few() 

