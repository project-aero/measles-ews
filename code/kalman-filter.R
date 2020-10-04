#!/usr/bin/env R

library(deSolve)
library(dplyr)
library(pomp)

pob <- readRDS("measles-pomp-object-Niamey.RDS")

case_data <- tibble(time = pob@times, reports = pob@data["reports",])
cov_data <- bind_cols(tibble(time = pob@tcovar), as_tibble(pob@covar))

PsystemSEIR <- function(pvec, covf,
                        init.vars = list(xvec, Pmat),
                        time.steps = c(1995, 1995 + 1 / 52)){
  
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      N_t <- covf$N(t)
      mu_t <- covf$mu(t)
      xi1t <- covf$xi1(t)
      xi2t <- covf$xi2(t)
      xi3t <- covf$xi3(t)
      xi4t <- covf$xi4(t)
      xi5t <- covf$xi5(t)
      xi6t <- covf$xi6(t)
      beta_t <- beta_mu * (1 + exp(xi1t * b1 + xi2t * b2 + xi3t * b3 + xi4t * b4 + xi5t * b5 + xi6t * b6))
      eta <- 365 / 8  
      gamma <- 365 / 5
      
      dS <- N_t * mu_t * 0.3 - beta_t * S * I / N_t
      dE <- beta_t * S * I / N_t - eta * E
      dI <- iota + eta * E -  gamma * I
      dC <- gamma * I
      
      list(c(dS=dS, dE=dE, dI=dI, dC = dC))
    })
  }
  iv <- sapply(init.vars, unlist)
  deSolve::lsoda(iv$xvec, time.steps, PModel, pvec)
}

genfun <- function(y) {
  approxfun(cov_data$time, y)
}
covf <- apply(cov_data, 2, genfun)
pvec <- coef(pob)
xvec <- c(S=3e-2, E=1.6e-4, I=1.6e-4, C = 0)
Pmat <- matrix(0) 

init.vars <- list(xvec = xvec, Pmat = Pmat)             
PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars)
