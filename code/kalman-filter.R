#!/usr/bin/env R

library(deSolve)
library(dplyr)
library(pomp)

pob <- readRDS("measles-pomp-object-Niamey.RDS")

case_data <- tibble(time = pob@times, reports = pob@data["reports",])
cov_data <- bind_cols(tibble(time = pob@tcovar), as_tibble(pob@covar))

PsystemSEIR <- function(pvec, covf,
                        init.vars,
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
      S <- exp(lS)
      E <- exp(lE)
      I <- exp(lI)
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
      
      dlS <- (N_t * mu_t * 0.3 - beta_t * S * I / N_t) / S
      dlE <- (beta_t * S * I / N_t - eta * E) / E
      dlI <- (iota + eta * E -  gamma * I) / I
      dC <- (gamma * I)
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(c(dlS = dlS, dlE = dlE, dlI = dlI, dC = dC, dPss = dP[1,1], dPse = dP[1,2], 
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
init.vars <- c(lS=log(3e-2 * 6.2e5), lE=log(1.6e-4 * 6.2e5), lI=log(1.6e-4 * 6.2e5), C = 0,
               Pss = 1, Pse = 0, Psi = 0, Psc = 0, Pee = 1, Pei = 0, Pec = 0, Pii = 1, Pic = 0, Pcc = 1)

out <- PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars, time.steps = case_data$time)

iterate_f_and_P <- function(xhat, PN, pvec, covf, time.steps){
  P <- PN / covf$N(time.steps[1])
  xhat_trans <- c(log(xhat[c("S", "E", "I")]), xhat["C"])
  if(!all(is.finite(xhat_trans))) {
    browser()
  }
  names(xhat_trans)[1:3] <- c("lS", "lE", "lI")
  init.vars <- c(xhat_trans, Pss = P[1,1], Pse = P[1,2], 
  Psi = P[1,3], Psc = P[1,4], Pee = P[2,2], Pei = P[2,3], Pec = P[2,4], 
  Pii = P[3,3], Pic = P[3,4], Pcc = P[4,4])
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, covf, time.steps)[2, ]
  xhat_new <- c(exp(ret[c("lS", "lE", "lI")]), ret["C"])
  names(xhat_new)[1:3] <- c("S", "E", "I")
  P_new  <- with(as.list(ret),        
            rbind(c(Pss, Pse, Psi, Psc),
                  c(Pse, Pee, Pei, Pec),
                  c(Psi, Pei, Pii, Pic),
                  c(Psc, Pec, Pic, Pcc)))
  PN_new <- P_new * covf$N(time.steps[2])
  list(xhat = xhat_new, PN = PN_new)
}

xhat <-  c(S=(3e-2 * 6.2e5), E=(1.6e-4 * 6.2e5), I=(1.6e-4 * 6.2e5), C = 0)
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
z_1 <-  case_data$reports[-1][1]
H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
R <- z_1 * (1 - pvec["rho"])

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


T <- nrow(case_data[-1,])
z <- case_data$reports[-1]

ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat)
P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))

K[, 1] <- K_1
xhat_kkmo[, 1] <- xhat_1_0
xhat_kk[, 1] <- xhat_1_1
P_kk[, , 1] <- P_1_1
P_kkmo[, , 1] <- P_1_0
S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
ytilde_k[, 1] <- ytilde_1

for (i in seq(2, T)){
  xhat_init <- xhat_kk[, i - 1]
  xhat_init["C"] <- 0
  PNinit <- P_kk[,,i - 1]
  PNinit[, 4] <- PNinit[4, ] <- 0
  XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                            time.steps = case_data$time[c(i - 1, i)])
  xhat_kkmo[, i] <- XP$xhat
  P_kkmo[, , i] <- XP$PN
  #R <- xhat_kkmo["C", i] * pvec["rho"] * (1 - pvec["rho"])
  R <- max(5, z[i - 1] * (1 - pvec["rho"]))
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
  P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
}

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log(2 * pi))
}

log_lik(S, ytilde_k)

kfnll <-
  function(cdata,
           pvec,
           logit_beta_mu,
           logit_S0,
           logit_b2,
           logit_b3,
           xhat0 = structure(c(18600, 99.2, 99.2, 0), .Dim = c(4L, 1L), 
                             .Dimnames = list(c("S", "E", "I", "C"), NULL)),
           Phat0 = diag(c(1, 1, 1, 0)),
           just_nll = TRUE) {

    pvec["beta_mu"] <- scaled_expit(logit_beta_mu, a_beta_mu, b_beta_mu)
    pvec["b2"] <- scaled_expit(logit_b2, a_bpar, b_bpar)
    pvec["b3"] <- scaled_expit(logit_b3, a_bpar, b_bpar)
    xhat0["S", 1] <- scaled_expit(logit_S0, a_S0, b_S0)
    
    #print(c("                                     ", pvec["beta_mu"], xhat0["S", 1] / b_S0, pvec["b2"]))
    
    # Initialize
    z_1 <-  cdata$reports[-1][1]
    H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
    #R <- max(5, z[1] * pvec["tau"])
    R <- max(5, z_1 * (1 - pvec["rho"]))
    
    # Predict
    XP_1_0 <- iterate_f_and_P(xhat0[, 1], PN = Phat0, pvec = pvec, covf = covf,
                              time.steps = cdata$time[c(2, 3)])
    xhat_1_0 <- XP_1_0$xhat
    P_1_0 <- XP_1_0$PN
    # Update
    
    K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R)
    ytilde_1 <- z_1 - H %*% xhat_1_0
    xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
    P_1_1 <- (diag(4) - K_1 %*% H) %*% P_1_0
    
    ## Now calculate for each step in simulation
    
    T <- nrow(case_data[-1,])
    z <- cdata$reports[-1]
    
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    
    K[, 1] <- K_1
    xhat_kkmo[, 1] <- xhat_1_0
    xhat_kk[, 1] <- xhat_1_1
    P_kk[, , 1] <- P_1_1
    P_kkmo[, , 1] <- P_1_0
    S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
    ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
    ytilde_k[, 1] <- ytilde_1
    
    for (i in seq(2, T)){
      xhat_init <- xhat_kk[, i - 1]
      xhat_init["C"] <- 0
      PNinit <- P_kk[,,i - 1]
      PNinit[, 4] <- PNinit[4, ] <- 0
      XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                            time.steps = cdata$time[c(i - 1, i)])
      xhat_kkmo[, i] <- XP$xhat
      P_kkmo[, , i] <- XP$PN
      R <- max(5, z[i - 1] * pvec["tau"])
      R <- max(5, z[i - 1] * (1 - pvec["rho"]))
      #R <- max(5, xhat_kk["C", i - 1] * pvec["rho"] * (1 - pvec["rho"]))
      S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
      K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
      ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
      xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
    }
    
    nll <- 0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi))
    if (!just_nll){
      list(nll = nll, xhat_kkmo = xhat_kkmo, xhat_kk = xhat_kk, P_kk = P_kk, ytilde_k = ytilde_k)
    } else {
      nll
    }
  }

library(bbmle)

scaled_expit <- function(y, a, b){
  (b - a) * exp(y)/ (1 + exp(y)) + a
} 

scaled_logit <- function(x, a, b){
  log((x - a) / (b - x))
}
a_beta_mu <- 1
b_beta_mu <- 1000
a_S0 <- 0
b_S0 <- 628e3

a_bpar <- 0
b_bpar <- 4

pvec2 <- pvec
pvec2["rho"] <- 0.1
pvec2["b1"] <- 0.3
pvec2["b2"] <- 0.3
pvec2["b3"] <- 0.6
pvec2["b4"] <- 0.5
pvec2["b5"] <- 0.2
pvec2["b6"] <- 0

m0 <- mle2(minuslogl = kfnll, 
           start = list(logit_beta_mu = scaled_logit(516, a_beta_mu, b_beta_mu), 
                        logit_S0 = scaled_logit(70e3, a_S0, b_S0),
                        logit_b2 = scaled_logit(1.6, a_bpar, b_bpar),
                        logit_b3 = scaled_logit(1, a_bpar, b_bpar)), 
           method = "Nelder-Mead",
           skip.hessian = TRUE,
           control = list(reltol = 1e-4, trace = 2),
           data = list(cdata = case_data, pvec = pvec2))

p0 <- profile(m0)
confint(p0)
plot(p0, absVal = FALSE)
confint(p0)[2,]
scaled_expit(confint(p0)[2,], a_S0, b_S0) / b_S0

kfret <- kfnll(cdata = case_data, pvec = pvec2, logit_beta_mu = 0.1326351, 
               logit_S0 = -2.0008322, logit_b2 = -0.5584351, 
               logit_b3 = -1.0823994, just_nll = FALSE)



plot(case_data$time[-1], kfret$xhat_kkmo["C",] * pvec["rho"])
points(case_data$time[-1], kfret$xhat_kk["C",] * pvec["rho"], col = 2, pch = 2)
lines(case_data$time[-1], case_data$reports[-1])

plot(case_data$time[-1], kfret$xhat_kkmo["I",])
points(case_data$time[-1], kfret$xhat_kk["I",], col = 2, pch = 2)


tgrid <- seq(0, 1, length.out = 100) + 1995
ximat <- cbind(covf$xi1(tgrid),
               covf$xi2(tgrid),
               covf$xi3(tgrid),
               covf$xi4(tgrid),
               covf$xi5(tgrid),
               covf$xi6(tgrid))

matplot(tgrid, ximat)


seasgrid <- 1 + exp(ximat %*% pvec2[c("b1", "b2", "b3", "b4", "b5", "b6")])
R0grid <- scaled_expit(1.1434, a_beta_mu, b_beta_mu) * seasgrid / (365 / 5)
plot(tgrid, R0grid)

par(mfrow = c(2, 1))

plot(case_data$time[-1], kfret$xhat_kkmo["C",] * pvec["rho"], xlim = c(1995, 1997))
points(case_data$time[-1], kfret$xhat_kk["C",] * pvec["rho"], col = 2, pch = 2, xlim = c(1995, 1997))
lines(case_data$time[-1], case_data$reports[-1])

plot(case_data$time[-1], kfret$xhat_kkmo["I",], xlim = c(1995, 1997))
points(case_data$time[-1], kfret$xhat_kk["I",], col = 2, pch = 2, xlim = c(1995, 1997))


plot(kfret$xhat_kk, kfret$xhat_kkmo)
