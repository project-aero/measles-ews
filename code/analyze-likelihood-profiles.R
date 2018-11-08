# analyze-likelihood-profiles.R


# MCAP function -----------------------------------------------------------
#  From Ionides and King, Royal Society Interface 2017
#  DOI: 10.1098/rsif.2017.0126

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


# Function for plotting results -------------------------------------------

plot_ci <- function(full_profile, mcap_profile, mcap_out, xlab, line_color = "steelblue", lambda_num = 0.75){
  smooth_fit <- mcap_out$fit[,c(1,2)]
  quad_fit <- mcap_out$fit[,c(1,3)] %>%
    filter(quadratic >= min(full_profile$loglik))
  
  ggplot() +
    geom_point(data = full_profile, aes(x = value, y = loglik), shape = 1, size = 1, color = "grey50", alpha = 0.5) +
    geom_point(data = mcap_profile, aes(x = value, y = loglik), shape = 19, size = 1, alpha = 1) +
    geom_line(data = smooth_fit, aes(x = parameter, y = smoothed), color = line_color) +
    geom_line(data = quad_fit, aes(x = parameter, y = quadratic), color = line_color, linetype = 2) +
    geom_vline(aes(xintercept = mcap_out$ci[1]), color = line_color) +
    geom_vline(aes(xintercept = mcap_out$ci[2]), color = line_color) +
    geom_hline(aes(yintercept = mcap_out$delta_line), color = line_color) +
    labs(x = xlab, y = "profile log-likelihood") +
    ggtitle(paste0("95% C.I. = (", round(mcap_out$ci,2)[1], ", ", round(mcap_out$ci,2)[2], "); span = ", lambda_num)) +
    theme(plot.title = element_text(size = 10))
}


# Function for running MCAP given profile file ----------------------------

run_mcap <- function(file_name, lambda = 0.75){
  profile_data <- read.csv(file_name) %>%
    slice(2:n()) %>%
    drop_na()
  
  mcap_data <- profile_data %>%
    group_by(value, parameter) %>%
    summarise(loglik = pomp::logmeanexp(loglik))
    
  
  do_param <- unique(mcap_data$parameter)
  mcap_out <- mcap(lp = mcap_data$loglik, parameter = mcap_data$value, lambda = lambda, confidence = 0.99)
  
  return(list(profile_data, mcap_data, mcap_out))
}


# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggthemes)


# Niamey, beta profile ----------------------------------------------------

fname <- ("../results/loglik-profile-beta-Niamey.csv")
mcap_all <- run_mcap(fname)
beta_plot <- plot_ci(
  full_profile = mcap_all[[1]], 
  mcap_profile = mcap_all[[2]], 
  mcap_out = mcap_all[[3]], 
  xlab = expression(beta)
)


# Niamey, rho profile -----------------------------------------------------

fname <- ("../results/loglik-profile-rho-Niamey.csv")
mcap_all <- run_mcap(fname, lambda = 0.2)
rho_plot <- plot_ci(
  full_profile = mcap_all[[1]], 
  mcap_profile = mcap_all[[2]], 
  mcap_out = mcap_all[[3]], 
  xlab = expression(rho),
  lambda = 0.2
)


# Combine Niamey ----------------------------------------------------------

cowplot::plot_grid(beta_plot, rho_plot, ncol = 2)
