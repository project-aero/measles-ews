

library(spaero)
library(pomp)
library(deSolve)

sim <- create_simulator(process_model = "SIS",
                        transmission = "frequency-dependent")
coef(sim)["beta_par"] <- 30
sd <- pomp::simulate(sim, times = seq(0, 10, by = 1/52), seed = 2, as.data.frame = TRUE)
plot(I ~ time, data = sd, log = 'y')
lines(cases ~ time, data = sd)

RunDetermSIS <- function(beta=(R0 * gamma), eta=(imports * R0 / pop.size),
                         gamma= 24, imports=0,
                         pop.size=1e5, R0=30 / 24, init.vars=c(I=0, C=0),
                         time.steps=seq(0, 1 / 52, len=1)) {
  ## Solves for the trajectory of the deterministic  SIS model
  ##
  ## Args:
  ##   beta: numeric. The transmission rate.
  ##   eta: numeric. The rate of infection from outside.
  ##   gamma: numeric. The recovery rate.
  ##   imports: numeric. The expected number of imported cases.
  ##   pop.size: numeric. The population size.
  ##   R0: numeric. The basic reproduction number (with no vaccination).
  ##   init.vars: numeric with name 'I'. The initial number infected.
  ##   time.steps: numeric vector. Time points at which to solve for the
  ##     trajectory.
  ##
  ## Returns:
  ##   numeric matrix. The time steps and infected proportion at those times.
  
  R0 <- beta / gamma  # In case user provides beta instead of R0
  stopifnot(c(beta, eta, gamma, imports, pop.size, R0, init.vars) >= 0)
  stopifnot(c(init.vars) <= pop.size)
  stopifnot(diff(time.steps) > 0)
  
  parameters <- c(eta = eta, beta = beta, gamma = gamma, pop.size = pop.size)
  SISModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      dI <- beta * (pop.size - I) * I / pop.size + eta * (pop.size - I) - gamma * I
      dC <- gamma * I
      list(c(dI=dI, dC))
    })
  }
  lsoda(init.vars, time.steps, SISModel, parameters)
}

odesol <- RunDetermSIS(init.vars = c(I=100, C=0), time.steps = sd$time[sd$I > 100])
lines(I~time, data = odesol, col = 2)
lines(C~time, data = odesol, col = 3)

f <- function(x, ...){
  ret <- RunDetermSIS(init.vars = c(I=x[1], C=0), ..., time.steps = c(0, 1 / 52))[2, c("I", "C")]
  matrix(ret, ncol = 1)
}

DriftMatrixSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, eval.pars = c(I = 0)){
  I <- eval.pars['I']
  ret <- beta  - 2 * beta * I / pop.size - eta - gamma
  matrix(ret)
}

DiffusionMatrixSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, eval.pars = c(I = 0)){
  I <- eval.pars['I']
  ret <- beta * I / pop.size * (1 - I / pop.size) + eta * (1 - I / pop.size) + gamma * I / pop.size
  matrix(ret)
}

PsystemSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, 
                       init.vars = c(I = 0, C = 0, Pii = 1, Pic = 1, Pci = 1, Pcc = 1),
                       time.steps = c(0, 1 / 52)){
  
  parameters <- c(eta = eta, beta = beta, gamma = gamma, pop.size = pop.size)
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      dI <- beta * (pop.size - I) * I / pop.size + eta * (1 - I) - gamma * I
      dC <- gamma * I
      F <- rbind(c(beta - 2 * beta * I / pop.size - eta - gamma, 0),
                 c(gamma, 0))
      Pmat <- rbind(c(Pii, Pic),
                    c(Pci, Pcc))
      Qmat <- rbind(c(beta * I / pop.size * (1 - I / pop.size) + eta * (1 - I / pop.size) + gamma * I / pop.size, -gamma * I / pop.size),
                    c(-gamma * I / pop.size, gamma * I / pop.size))
      dPmat <-  F %*% Pmat + Pmat %*% t(F) + Qmat
      list(c(dI=dI, dC = dC, 
             dPii = dPmat[1,1], dPic = dPmat[1,2], 
             dPci = dPmat[2,1], dPcc = dPmat[2,2]))
    })
  }
  lsoda(init.vars, time.steps, PModel, parameters)
}

iterate_P <- function(xhat, P, ...){
  ret <- PsystemSIS(init.vars = c(I = xhat[1], C = 0, 
                           Pii = P[1,1] / sqrt(1e5), Pic = 0, 
                           Pci = 0, Pcc = 0), ...)[2, c("Pii", "Pic", "Pci", "Pcc")]
  rbind(c(ret[1], ret[2]),
        c(ret[3], ret[4])) * sqrt(1e5)
}

iterate_f_and_P <- function(xhat, P, pop.size = 1e5, ...){
  ret <- PsystemSIS(init.vars = c(I = xhat[1], 
                                  C = 0, 
                                  Pii = P[1,1] / sqrt(pop.size), 
                                  Pic = 0, 
                                  Pci = 0, 
                                  Pcc = 0), pop.size = pop.size, ...)[2, c("I", "C", "Pii", "Pic", "Pci", "Pcc")]
  list(xhat = matrix(c(ret[1], ret[2]), ncol = 1), 
       Phat = rbind(c(ret[3], ret[4]),
                    c(ret[5], ret[6])) * sqrt(pop.size))
}


# Initialize
xhat0 <- matrix(c(0, 0), ncol = 1)
Phat0 <- rbind(c(1, 0),
               c(0, 0))
dt <- 1 / 52
z_1 <- sd$reports[1]
H <- matrix(c(0, coef(sim)["rho"]), ncol = 2)
R <- sd$cases * H[2] * (1 - H[2])
R[R < 1] <- 1

# Predict
xhat_1_0 <- f(xhat0, beta = 30)
P_1_0 <- iterate_P(xhat0, Phat0, beta = 30)

# Update

K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
ytilde_1 <- z_1 - H %*% xhat_1_0
xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
P_1_1 <- (diag(2) - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(sd)
z <- sd$reports

ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(2, T))
P_kk <- P_kkmo <- array(NA_real_, dim = c(2, 2, T))

K[, 1] <- K_1
xhat_kkmo[, 1] <- xhat_1_0
xhat_kk[, 1] <- xhat_1_1
P_kk[, , 1] <- P_1_1
P_kkmo[, , 1] <- P_1_0
S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R[1]
ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
ytilde_k[, 1] <- ytilde_1

for (i in seq(2, T)){
  xhat_kkmo[, i] <- f(xhat_kk[, i - 1], beta = 30)
  P_kkmo[, , i] <- iterate_P(xhat_kk[, i - 1], P_kk[, , i - 1], beta = 30)
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R[i]
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  P_kk[, , i] <- (1 - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
}

plot(I ~ time, data = sd)
points(cases ~ time, data = sd, col = "orange")
points(reports ~ time, data = sd, col = "red")
lines(sd$time, xhat_kk[1,], col = "blue")
lines(sd$time, xhat_kk[2,], col = "blue")
legend("topleft", col = c("black", "orange", "red", "blue"),
       legend = c("No. infecteds", "All cases", "Reported cases", "Kalman filter"),
       pch = c(1, 1, 1, NA), lty = c(NA, NA, NA, 1))

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log(2 * pi))
}

log_lik(S, ytilde_k)

# Now create  a higher-level function for the log likihood

kfnll <-
  function(z,
           beta = 30,
           rho = 0.1,
           gamma = 24,
           dt = 1 / 52,
           xhat0 = c(0, 0),
           Phat0 = rbind(c(1, 0), 
                         c(0, 0))) {
    z_1 <- z[1]
    H <- matrix(c(0, rho), ncol = 2)
    
    # Predict
    xP <- iterate_f_and_P(xhat0, Phat0, beta = beta, gamma = gamma)
    xhat_1_0 <- xP$xhat
    PP_1_0 <- xP$Phat
    # Update
    
    K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
    ytilde_1 <- z_1 - H %*% xhat_1_0
    xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
    P_1_1 <- (diag(2) - K_1 %*% H) %*% P_1_0
    
    T <- length(z)
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(2, T))
    P_kk <- P_kkmo <- array(NA_real_, dim = c(2, 2, T))

    K[, 1] <- K_1
    xhat_kkmo[, 1] <- xhat_1_0
    xhat_kk[, 1] <- xhat_1_1
    P_kk[, , 1] <- P_1_1
    P_kkmo[, , 1] <- P_1_0
    Rc <- xhat_kkmo[2, 1] * rho * (1 - rho)
    if(Rc < 1){
      Rc <- 1
    }
    S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + Rc
    ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
    ytilde_k[, 1] <- ytilde_1
    
    for (i in seq(2, T)){
      xP <- iterate_f_and_P(xhat_kk[, i - 1], P_kk[, , i - 1], beta = beta, gamma = gamma)
      xhat_kkmo[, i] <- xP$xhat
      P_kkmo[, , i] <- xP$Phat
      Rc <- xhat_kkmo[2, i] * rho * (1 - rho)
      if(Rc < 1){
        Rc <- 1
      }
      S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + Rc
      K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
      ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
      xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
      P_kk[, , i] <- (1 - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
    }
    
    0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi))
  }

kfnll(sd$reports)

## BBMLE estimates

library(bbmle)


m0 <- mle2(minuslogl = kfnll, start = list(rho = .2), 
           method = "L-BFGS-B", 
           lower = 0.05, upper = .95, #if R0 goes below 1, then the lsoda calls have problems
           trace = TRUE, data = list(z = sd$reports))

p0 <- profile(m0)
confint(p0)
plot(p0, absVal = FALSE)
