

library(spaero)
library(pomp)

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
  ##   beta: numeric. The transmission rate without vaccination.
  ##   eta: numeric. The rate of infection from outside.
  ##   gamma: numeric. The recovery rate.
  ##   imports: numeric. The expected number of imported cases.
  ##   pop.size: numeric. The population size.
  ##   R0: numeric. The basic reproduction number (with no vaccination).
  ##   init.vars: numeric with name 'I'. The initial fraction infected.
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
