

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
                         pop.size=1e5, R0=30 / 24, init.vars=c(I=0),
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
      dI <- beta * (pop.size - I) * I / pop.size + eta * (1 - I) - gamma * I
      list(dI=dI)
    })
  }
  lsoda(init.vars, time.steps, SISModel, parameters)
}

odesol <- RunDetermSIS(init.vars = c(I=100), time.steps = sd$time[sd$I > 100])
lines(I~time, data = odesol, col = 2)

f <- function(x){
  RunDetermSIS(init.vars = c(I=x[1]), time.steps = c(0, 1 / 52))[2, "I"]
}

DriftMatrixSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, eval.pars = c(I = 0)){
  I <- eval.pars['I']
  ret <- beta  - 2 * beta * I / pop.size - eta - gamma
  matrix(ret)
}

DiffusionMatrixSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, eval.pars = c(I = 0)){
  I <- eval.pars['I']
  ret <- beta * I * (pop.size - I) + eta * (pop.size - I) + gamma * I
  matrix(ret)
}

PsystemSIS <- function(beta = 30, eta = 0, gamma = 24, pop.size = 1e5, init.vars = c(I = 0, P = 1),
                       time.steps = c(0, 1 / 52)){
  
  parameters <- c(eta = eta, beta = beta, gamma = gamma, pop.size = pop.size)
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      dI <- beta * (pop.size - I) * I / pop.size + eta * (1 - I) - gamma * I
      F <- (beta - 2 * beta * I / pop.size - eta - gamma)
      Q <- beta * I * (pop.size - I) + eta * (pop.size - I) + gamma * I
      dP <- 2 * F * P + Q
      list(c(dI=dI, dP = dP))
    })
  }
  lsoda(init.vars, time.steps, PModel, parameters)
}

iterate_P <- function(xhat, P){
  PsystemSIS( init.vars = c(I = xhat, P = P))[2, "P"]
}

# Initialize
xhat0 <- 0
Phat0 <- 1
dt <- 1 / 52
z_1 <- sd$reports[1]
H <- coef(sim)["rho"] * dt * coef(sim)["gamma"]
R <- 1

# Predict
xhat_1_0 <- f(xhat0)
P_1_0 <- iterate_P(xhat0, Phat0)

# Update

K_1 = P_1_0 * t(H) %*% solve(H %*% P_1_0 %*% t(H) + R)
xhat_1_1 <- xhat_1_0 + K_1 %*% (z_1 - H * xhat_1_0)
P_1_1 <- (1 - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(sd)
z <- sd$reports
ytilde_kk <- K <- xhat_kk <- xhat_kkmo <- P_kk <- P_kkmo <- S <- numeric(T)

K[1] <- K_1
xhat_kkmo[1] <- xhat_1_0
xhat_kk[1] <- xhat_1_1
P_kk[1] <- P_1_1
P_kkmo[1] <- P_1_0
S[1] <- H %*% P_kkmo[1] %*% t(H) + R
ytilde_kk[1] <- z[1] - H %*% xhat_kk[1]

for (i in seq(2, T)){
  xhat_kkmo[i] <- f(xhat_kk[i - 1])
  P_kkmo[i] <- iterate_P(xhat_kk[i - 1], P_kk[i - 1])
  S[i] <- H %*% P_kkmo[i] %*% t(H) + R
  K[i] <- P_kkmo[i] * t(H) %*% solve(S[i])
  xhat_kk[i] <- xhat_kkmo[i] + K[i] %*% (z[i] - H * xhat_kkmo[i])
  P_kk[i] <- (1 - K[i] %*% H) %*% P_kkmo[i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[i]
}


plot(I ~ time, data = sd)
points(cases ~ time, data = sd, col = "orange")
points(sd$time, xhat_kk, col = "blue")
plot(I ~ time, data = sd, log = 'y')
lines(cases ~ time, data = sd)
legend("topleft", col =c ("black", "orange", "blue"), 
       legend = c("No. infecteds", "Reported cases", "Kalman filter"),
       pch = 1)

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log( 2 * pi))
}

log_lik(S, ytilde_kk)


