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
      P_t <- covf$N(t)
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
      S <- S_0 + exp(Nt_birth) - exp(Nt_trans)
      C <- exp(Nt_cases) - 1
      E <- E_0 + exp(Nt_trans) - exp(Nt_prog)
      I <- I_0 + exp(Nt_prog) - exp(Nt_recov)
      if (any(c(S, C, E, I) < 0)){
        browser()
      }
      lambda <- c(beta_t / P_t * S * I + iota,
                  P_t * mu_t * 0.3,
                  eta * E,
                  gamma * I,
                  rho * gamma * I)
      ito_factors_make <- function(Nt){
        exp(-Nt) - 0.5 * exp(-2 * Nt)
      }
      Ntv <- c(Nt_trans, Nt_birth, Nt_prog, Nt_recov, Nt_cases)
      itof <- sapply(Ntv, ito_factors_make)
      dNt_trans <- itof[1] * lambda[1]
      dNt_birth <- itof[2] * lambda[2]
      dNt_prog <- itof[3] * lambda[3]
      dNt_recov <- itof[4] * lambda[4]
      dNt_cases <- itof[5] * lambda[5]
      
      ito_factors_derv_make <- function(Nt){
        exp(-2 * Nt) - exp(-Nt)
      }
      itofp <- sapply(Ntv, ito_factors_derv_make)
      F <- rbind(c(itofp[1] * lambda[1] + itof[1] * beta_t / P_t * (-1) * exp(Ntv[1]) * I, 
                      itof[1] * beta_t / P_t * exp(Ntv[2]) * I, 
                      itof[1] * beta_t / P_t * S * exp(Ntv[3]), 
                      itof[1] * beta_t / P_t * S * (-1) * exp(Ntv[4]),
                      0),
                 c(0, itofp[2] * lambda[2], 0, 0, 0),
                 c(itof[3] * eta * exp(Ntv[1]), 0, itofp[3] * lambda[3] + itof[3] * eta * -1 * exp(Ntv[3]), 0, 0),
                 c(0, 0, itof[4] * gamma * exp(Ntv[3]), itofp[4] * lambda[4] + itof[4] * gamma * -1 * exp(Ntv[4]), 0),
                 c(0, 0, itof[5] * rho * gamma * exp(Ntv[3]), itof[5] * rho * gamma * -1 * exp(Ntv[4]), itofp[5] * lambda[5])) 

      Q <- diag(exp(-Ntv) * lambda)
      
      P <- rbind(c(P11, P12, P13, P14, P15),
                 c(P12, P22, P23, P24, P25),
                 c(P13, P23, P33, P34, P35),
                 c(P14, P24, P34, P44, P45),
                 c(P15, P25, P35, P45, P55))
      
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(c(dNt_trans = dNt_trans, 
             dNt_birth = dNt_birth,
             dNt_prog = dNt_prog,
             dNt_recov = dNt_recov,
             dNt_cases = dNt_cases,
             dP11 = dP[1,1], 
             dP12 = dP[1,2], 
             dP13 = dP[1,3], 
             dP14 = dP[1,4], 
             dP15 = dP[1,5],
             dP22 = dP[2,2], 
             dP23 = dP[2,3], 
             dP24 = dP[2,4],
             dP25 = dP[2,5],
             dP33 = dP[3,3], 
             dP34 = dP[3,4],
             dP35 = dP[3,5],
             dP44 = dP[4,4],
             dP45 = dP[4,5],
             dP55 = dP[5,5]))
    })
  }
  deSolve::lsoda(init.vars, time.steps, PModel, pvec)
}

genfun <- function(y) {
  approxfun(cov_data$time, y)
}
covf <- apply(cov_data, 2, genfun)
pvec <- coef(pob)

pvec["S_0"] <- 0.11 * 628e3
pvec["E_0"] <- 10
pvec["I_0"] <- 30

init.vars <- c(Nt_trans = log(10 + 1),
               Nt_birth = log(10 + 1),
               Nt_prog = log(10 + 1),
               Nt_recov = log(10 + 1),
               Nt_cases = 0,
               P11 = 1,
               P12 = 0,
               P13 = 0,
               P14 = 0,
               P15 = 0,
               P22 = 1,
               P23 = 0,
               P24 = 0,
               P25 = 0,
               P33 = 1,
               P34 = 0,
               P35 = 0,
               P44 = 1,
               P45 = 0,
               P55 = 1)

out <- PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars, 
                   time.steps = case_data$time)

Ssol <- pvec["S_0"] + (exp(out[, "Nt_birth"]) - 1) - (exp(out[, "Nt_trans"]) - 1)
Esol <- pvec["E_0"] + (exp(out[, "Nt_trans"]) - 1) - (exp(out[, "Nt_prog"]) - 1)
Isol <- pvec["I_0"] + (exp(out[, "Nt_prog"]) - 1) - (exp(out[, "Nt_recov"]) - 1)

plot(case_data$time, Ssol)
plot(case_data$time, Esol, log = "y")
plot(case_data$time, Isol)
plot(Ssol, Isol, log = "xy", type = 'b')
plot(case_data$time, out[, "P11"], type = 'l', log = "y")
plot(case_data$time, out[, "P12"])

iterate_f_and_P <- function(xhat, P, pvec, covf, time.steps){

  init.vars <- c(xhat["Nt_trans"],
                 xhat["Nt_birth"],
                 xhat["Nt_prog"],
                 xhat["Nt_recov"],
                 Nt_cases = 0,
                 P11 = P[1,1],
                 P12 = P[1,2],
                 P13 = P[1,3],
                 P14 = P[1,4],
                 P15 = 0,
                 P22 = P[2,2],
                 P23 = P[2,3],
                 P24 = P[2,4],
                 P25 = 0,
                 P33 = P[3,3],
                 P34 = P[3,4],
                 P35 = 0,
                 P44 = P[4,4],
                 P45 = 0,
                 P55 = 0)
  
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, covf, time.steps)[2, ]
  xhat_new <- ret[c("Nt_trans", "Nt_birth", "Nt_prog", "Nt_recov", "Nt_cases")] 
  P_new <- with(as.list(ret),
            rbind(c(P11, P12, P13, P14, P15),
             c(P12, P22, P23, P24, P25),
             c(P13, P23, P33, P34, P35),
             c(P14, P24, P34, P44, P45),
             c(P15, P25, P35, P45, P55)))

  list(xhat = xhat_new, P = P_new)
}

xhat <-  c(Nt_trans = 2.39789527279837, 
           Nt_birth = 2.39789527279837, 
           Nt_prog = 2.39789527279837, 
           Nt_recov = 2.39789527279837, 
           Nt_cases = 0)

P <- with(as.list(init.vars),
              rbind(c(P11, P12, P13, P14, P15),
                    c(P12, P22, P23, P24, P25),
                    c(P13, P23, P33, P34, P35),
                    c(P14, P24, P34, P44, P45),
                    c(P15, P25, P35, P45, P55)))

iterate_f_and_P(xhat = xhat, P = P, pvec = pvec, covf = covf, 
                time.steps = c(1996, 1996 + 1 / 52))

keep_valid <- function(xh, pr){
  ## constrain to valid S, E, I values
  S <- pr["S_0"] + exp(xh["Nt_birth"]) - exp(xh["Nt_trans"])
  if (S < 0) {
    S <- 0
    xh["Nt_trans"] <- log(pr["S_0"] + exp(xh["Nt_birth"]))
  }
  E <- pr["E_0"] + exp(xh["Nt_trans"]) - exp(xh["Nt_prog"])
  if (E < 0) {
    E <- 0
    xh["Nt_prog"] <- log(pr["E_0"] + exp(xh["Nt_trans"]))
  }
  I <- pr["I_0"] + exp(xh["Nt_prog"]) - exp(xh["Nt_recov"])
  if (I < 0){
    xh["Nt_recov"] <- log(pr["I_0"] + exp(xh["Nt_prog"]))
  }
  xh
}


# Initialize
xhat0 <- matrix(xhat, ncol = 1)
rownames(xhat0) <- names(xhat)
Phat0 <- P
z_1 <-  log(case_data$reports[-1][1] + 1)
H <- matrix(c(0, 0, 0, 0, 1), ncol = 5)
R <- 1#z_1 * pvec["tau"]

# Predict
XP_1_0 <- iterate_f_and_P(xhat0[, 1], 
                          P = Phat0, 
                          pvec = pvec, 
                          covf = covf,
                          time.steps = case_data$time[c(2, 3)])
xhat_1_0 <- keep_valid(XP_1_0$xhat, pvec)
P_1_0 <- XP_1_0$P
# Update

K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
ytilde_1 <- z_1 - H %*% xhat_1_0
xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
P_1_1 <- (diag(5) - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(case_data[-1,])
z <- log(case_data$reports[-1] + 1)

ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(5, T))
rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat)
P_kk <- P_kkmo <- array(NA_real_, dim = c(5, 5, T))

K[, 1] <- K_1
xhat_kkmo[, 1] <- xhat_1_0
xhat_kk[, 1] <- xhat_1_1
P_kk[, , 1] <- P_1_1
P_kkmo[, , 1] <- P_1_0
S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
ytilde_k[, 1] <- ytilde_1

for (i in seq(2, T)){
  XP <- iterate_f_and_P(xhat = xhat_kk[, i - 1], 
                        P = P_kk[,, i -1], 
                        pvec = pvec, 
                        covf = covf,
                        time.steps = case_data$time[c(i - 1, i)])
  xhat_kkmo[, i] <- XP$xhat
  P_kkmo[, , i] <- XP$P
  #R <- xhat_kkmo["C", i] * pvec["rho"] * (1 - pvec["rho"])
  #R <- max(1, z[i - 1] * pvec["tau"])
  R <- 1
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  xhat_kk[, i] <- keep_valid(xhat_kk[, i], pvec)
  P_kk[, , i] <- (diag(5) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]

}

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log(2 * pi))
}

log_lik(S, ytilde_k)

Ssol <- pvec["S_0"] + (exp(xhat_kkmo["Nt_birth",]) - 1) - (exp(xhat_kkmo["Nt_trans",]) - 1)
Esol <- pvec["E_0"] + (exp(xhat_kkmo["Nt_trans",]) - 1) - (exp(xhat_kkmo["Nt_prog",]) - 1)
Isol <- pvec["I_0"] + (exp(xhat_kkmo["Nt_prog",]) - 1) - (exp(xhat_kkmo["Nt_recov",]) - 1)
csol <- 0 + exp(xhat_kkmo["Nt_cases",]) - 1
csol2 <- 0 + exp(xhat_kk["Nt_cases", ]) - 1
Isol2 <- pvec["I_0"] + (exp(xhat_kk["Nt_prog",]) - 1) - (exp(xhat_kk["Nt_recov",]) - 1)

kfnll <-
  function(cdata,
           pvec,
           logit_beta_mu,
           logit_S0,
           logit_I0,
           logit_E0,
           logit_b1,
           logit_b2,
           logit_b3,
           logit_b4, 
           logit_b5,
           logit_b6,
           logit_rho,
           logit_iota,
           logit_tau,
           logit_tau2, 
           xhat0 = structure(c(18600, 99.2, 99.2, 0), .Dim = c(4L, 1L), 
                             .Dimnames = list(c("S", "E", "I", "C"), NULL)),
           Phat0 = diag(c(1, 1, 1, 0)),
           just_nll = TRUE) {

    pvec["beta_mu"] <- scaled_expit(logit_beta_mu, a_beta_mu, b_beta_mu)
    pvec["b1"] <- scaled_expit(logit_b1, a_bpar, b_bpar)
    pvec["b2"] <- scaled_expit(logit_b2, a_bpar, b_bpar)
    pvec["b3"] <- scaled_expit(logit_b3, a_bpar, b_bpar)
    pvec["b4"] <- scaled_expit(logit_b4, a_bpar, b_bpar)
    pvec["b5"] <- scaled_expit(logit_b5, a_bpar, b_bpar)
    pvec["b6"] <- scaled_expit(logit_b6, a_bpar, b_bpar)
    pvec["rho"] <- scaled_expit(logit_rho, a_rho, b_rho)
    pvec["iota"] <- scaled_expit(logit_iota, a_iota, b_iota)
    pvec["tau"] <- scaled_expit(logit_tau, a_tau, b_tau)
    xhat0["S", 1] <- scaled_expit(logit_S0, a_S0, b_S0)
    xhat0["I", 1] <- scaled_expit(logit_I0, a_I0, b_I0)
    xhat0["E", 1] <- scaled_expit(logit_E0, a_E0, b_E0)
    tau2 <- scaled_expit(logit_tau2, a_tau2, b_tau2)
    
    #print(c("                                     ", pvec["beta_mu"], xhat0["S", 1] / b_S0, pvec["b2"]))
    
    # Initialize
    z_1 <-  cdata$reports[-1][1]
    H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
    #R <- max(5, z[1] * pvec["tau"])
    R <- max(tau2, z_1  * (1 + 1 /  pvec["tau"])) 
    
    
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
      #R <- max(5, z[i - 1] * pvec["tau"])
      #R <- max(.5, z[i - 1] * (1 - pvec["rho"]))
      R <- max(tau2, xhat_kkmo["C", i] * pvec["rho"] * (1 + 1 /  pvec["tau"]))
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
      list(nll = nll, xhat_kkmo = xhat_kkmo, xhat_kk = xhat_kk, 
           P_kkmo = P_kkmo, P_kk = P_kk, 
           ytilde_k = ytilde_k, S = S)
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
a_I0 <- 0
b_I0 <- 100
a_E0 <- 0
b_E0 <- 100

a_bpar <- 0
b_bpar <- 5

a_rho <- 0
b_rho <- 1

a_iota <- 0
b_iota <- 100

a_tau <- 0.001
b_tau <- 1000

a_tau2 <- 0
b_tau2 <- 20

pvec2 <- pvec
pvec2["rho"] <- 0.1
pvec2["b1"] <- 0.3
pvec2["b2"] <- 0.3
pvec2["b3"] <- 0.6
pvec2["b4"] <- 0.5
pvec2["b5"] <- 0.2
pvec2["b6"] <- 0

Phat0 <- diag(c(1e4, 1e2, 1e2, 0))


system.time(m0 <- mle2(minuslogl = kfnll, 
           start = list(logit_beta_mu = scaled_logit(234, a_beta_mu, b_beta_mu), 
                        logit_S0 = scaled_logit(46098, a_S0, b_S0),
                        logit_I0 = scaled_logit(50, a_I0, b_I0),
                        logit_E0 = scaled_logit(25, a_E0, b_E0),
                        logit_b1 = scaled_logit(1.56, a_bpar, b_bpar),
                        logit_b2 = scaled_logit(2.81, a_bpar, b_bpar),
                        logit_b3 = scaled_logit(1.03, a_bpar, b_bpar),
                        logit_b4 = scaled_logit(0.107, a_bpar, b_bpar),
                        logit_b5 = scaled_logit(0.0134, a_bpar, b_bpar),
                        logit_b6 = scaled_logit(3.9512, a_bpar, b_bpar),
                        logit_rho = scaled_logit(0.30, a_rho, b_rho),
                        logit_iota = scaled_logit(10, a_iota, b_iota),
                        logit_tau = scaled_logit(0.1, a_tau, b_tau),
                        logit_tau2 = scaled_logit(5, a_tau2, b_tau2)),
           method = "Nelder-Mead",
           skip.hessian = TRUE,
           control = list(reltol = 1e-4, trace = 1, maxit = 1000),
           data = list(cdata = case_data, pvec = pvec2, Phat0 = Phat0)))

#p0 <- profile(m0)
#confint(p0)
#plot(p0, absVal = FALSE)
#confint(p0)[2,]

scaled_expit(coef(m0)["logit_S0"], a_S0, b_S0)
scaled_expit(coef(m0)["logit_I0"], a_I0, b_I0)
scaled_expit(coef(m0)["logit_E0"], a_E0, b_E0)
(rho_hat <- scaled_expit(coef(m0)["logit_rho"], a_rho, b_rho))
scaled_expit(coef(m0)["logit_iota"], a_iota, b_iota)
scaled_expit(coef(m0)["logit_tau"], a_tau, b_tau)
scaled_expit(coef(m0)["logit_tau2"], a_tau2, b_tau2)

kfret <- with(as.list(coef(m0)), 
              kfnll(cdata = case_data, pvec = pvec2, 
               logit_beta_mu = logit_beta_mu, 
               logit_S0 = logit_S0,
               logit_E0 = logit_E0,
               logit_I0 = logit_I0,
               logit_b1 = logit_b1,
               logit_b2 = logit_b2, 
               logit_b3 = logit_b3,
               logit_b4 = logit_b4,
               logit_b5 = logit_b5,
               logit_b6 = logit_b6,
               logit_rho = logit_rho,
               logit_iota = logit_iota,
               logit_tau = logit_tau,
               logit_tau2 = logit_tau2,
               just_nll = FALSE))

par(mfrow = c(1, 1))
test <- case_data$time > 1990
qqnorm(kfret$ytilde_k[test]/ kfret$S[test]) # evalutate departure from normality
abline(0, 1)

par(mfrow = c(4, 1))
plot(case_data$time[-1], kfret$xhat_kkmo["C",] * rho_hat)
points(case_data$time[-1], kfret$xhat_kk["C",] * rho_hat, col = 2, pch = 2)
lines(case_data$time[-1], case_data$reports[-1])
plot(case_data$time[-1], kfret$S, log = "y")
plot(case_data$time[-1], kfret$ytilde_k)
plot(case_data$time[-1], kfret$ytilde_k / kfret$S)

par(mfrow = c(2, 1))
plot(case_data$time[-1], kfret$xhat_kkmo["I",])
points(case_data$time[-1], kfret$xhat_kk["I",], col = 2, pch = 2)

plot(case_data$time[-1], kfret$xhat_kkmo["S",])
points(case_data$time[-1], kfret$xhat_kk["S",], col = 2, pch = 2)

tgrid <- seq(0, 1, length.out = 100) + 1995
ximat <- cbind(covf$xi1(tgrid),
               covf$xi2(tgrid),
               covf$xi3(tgrid),
               covf$xi4(tgrid),
               covf$xi5(tgrid),
               covf$xi6(tgrid))

matplot(tgrid, ximat)

is_spline_par <- grepl("^logit_b[1-6]$", names(coef(m0)))
bhat <- scaled_expit(coef(m0)[is_spline_par], a_bpar, b_bpar)
seasgrid <- 1 + exp(ximat %*% bhat)
beta_mu_hat <- scaled_expit(coef(m0)["logit_beta_mu"], a_beta_mu, b_beta_mu)
R0grid <- beta_mu_hat * seasgrid / (365 / 5)
plot(tgrid, R0grid)

plot(log(case_data$reports[-1] + 1), log(rho_hat * kfret$xhat_kkmo["C",] + 1))
cor(case_data$reports[-1], rho_hat * kfret$xhat_kkmo["C",]) ^ 2

1 - sum((case_data$reports[-1] - rho_hat * kfret$xhat_kkmo["C", ])^2) / sum((case_data$reports[-1] - mean(case_data$reports[-1])) ^ 2)

upper95 <- kfret$xhat_kkmo["C",] + sqrt(kfret$P_kkmo[4, 4, ]) * 1.96
lower95 <- kfret$xhat_kkmo["C",] - sqrt(kfret$P_kkmo[4, 4, ]) * 1.96
lower95 <- ifelse(lower95 < 0, 0, lower95)

plot(case_data$time[-1], (rho_hat * upper95) ^ .5, type = 'l')
lines(case_data$time[-1], (rho_hat * lower95) ^ .5, type = 'l')
lines(case_data$time[-1], case_data$reports[-1] ^ .5, col = 2)
mean(case_data$reports[-1] >= rho_hat * lower95 & case_data$reports[-1] < rho_hat * upper95)
