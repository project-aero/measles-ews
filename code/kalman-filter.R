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
      
      F <- rbind(c(-beta_t * I / N_t,    0, -beta_t * S / N_t, 0),
                 c( beta_t * I / N_t, -eta,  beta_t * S / N_t, 0),
                 c(                0,  eta,            -gamma, 0),
                 c(                0,    0,             gamma, 0))
      
      
      f <- c(mu_t * 0.3, beta_t * (S / N_t) * (I / N_t), eta * E / N_t, gamma * I / N_t)
      Q <- rbind(c(f[1] + f[2],       -f[2],           0,     0),
                 c(      -f[2], f[2] + f[3],       -f[3],     0),
                 c(          0,       -f[3], f[3] + f[4], -f[4]),
                 c(          0,           0,       -f[4],  f[4]))
      
      P <- rbind(c(Pss, Pse, Psi, Psc),
                 c(Pse, Pee, Pei, Pec),
                 c(Psi, Pei, Pii, Pic),
                 c(Psc, Pec, Pic, Pcc))
      
      dS <- N_t * mu_t * 0.3 - beta_t * S * I / N_t
      dE <- beta_t * S * I / N_t - eta * E
      dI <- iota + eta * E -  gamma * I
      dC <- gamma * I
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(c(dS=dS, dE=dE, dI=dI, dC = dC, dPss = dP[1,1], dPse = dP[1,2], 
             dPsi = dP[1,3], dPsc = dP[1,4], dPee = dP[2,2], dPei = dP[2,3], 
             dPec = dP[2,4], dPii = dP[3,3], dPic = dP[3,4], dPcc = dP[4,4]))
    })
  }
  deSolve::lsoda(init.vars, time.steps, PModel, pvec)
}

genfun <- function(y) {
  approxfun(cov_data$time, y)
}
covf <- apply(cov_data, 2, genfun)
pvec <- coef(pob)
init.vars <- c(S=3e-2, E=1.6e-4, I=1.6e-4, C = 0,
               Pss = 1, Pse = 0, Psi = 0, Psc = 0, Pee = 1, Pei = 0, Pec = 0, Pii = 1, Pic = 0, Pcc = 1)

out <- PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars, time.steps = case_data$time)

iterate_f_and_P <- function(xhat, PN, pvec, covf, time.steps){
  P <- P / covf$N(time.steps[1])
  init.vars <- c(xhat, Pss = P[1,1], Pse = P[1,2], 
  Psi = P[1,3], Psc = P[1,4], Pee = P[2,2], Pei = P[2,3], Pec = P[2,4], 
  Pii = P[3,3], Pic = P[3,4], Pcc = P[4,4])
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, covf, time.steps)[2, ]
  xhat <- ret[c("S", "E", "I", "C")]
  P <- with(as.list(init.vars),        
            rbind(c(Pss, Pse, Psi, Psc),
                  c(Pse, Pee, Pei, Pec),
                  c(Psi, Pei, Pii, Pic),
                  c(Psc, Pec, Pic, Pcc)))
  PN <- P * covf$N(time.steps[2])
  list(xhat = xhat, PN = PN)
}

xhat <- init.vars[1:4]
P <- with(as.list(init.vars[-c(1:4)]),        
          rbind(c(Pss, Pse, Psi, Psc),
                c(Pse, Pee, Pei, Pec),
                c(Psi, Pei, Pii, Pic),
                c(Psc, Pec, Pic, Pcc)))

iterate_f_and_P(xhat = xhat, PN = P, pvec = pvec, covf = covf, 
                time.steps = c(1996, 1996 + 1 / 52))

# Initialize
xhat0 <- matrix(xhat, ncol = 1)
rownames(xhat0) <- names(xhat)
Phat0 <- P
z_1 <-  case_data$reports[2]
H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
R <- 1

# Predict
XP_1_0 <- iterate_f_and_P(xhat0[, 1], PN = Phat0, pvec = pvec, covf = covf,
                         time.steps = case_data$time[c(2, 3)])
xhat_1_0 <- XP_1_0$xhat
P_1_0 <- XP_1_0$PN
# Update

K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
ytilde_1 <- z_1 - H %*% xhat_1_0
xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
P_1_1 <- (diag(4) - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(case_data)
z <- case_data$reports

ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
rownames(xhat_kk) <- names(xhat)
P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))

K[, 2] <- K_1
xhat_kkmo[, 2] <- xhat_1_0
xhat_kk[, 2] <- xhat_1_1
P_kk[, , 2] <- P_1_1
P_kkmo[, , 2] <- P_1_0
S[, 2] <- H %*% P_kkmo[, , 1] %*% t(H) + R
ytilde_kk[, 2] <- z[2] - H %*% xhat_kk[, 2]
ytilde_k[, 1] <- ytilde_1

for (i in seq(3, T)){
  xhat_init <- xhat_kk[, i - 1]
  xhat_init["C"] <- 0
  PNinit <- P_kk[,,i - 1]
  PNinit[, 4] <- PNinit[4, ] <- 0
  XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                            time.steps = case_data$time[c(i - 1, i)])
  xhat_kkmo[, i] <- XP$xhat
  P_kkmo[, , i] <- XP$PN
  
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
}
