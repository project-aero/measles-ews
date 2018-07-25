# Example of fitting a stochastic SIR model using `pomp`


# Load packages -----------------------------------------------------------

library(tidyverse)
library(pomp)


# Define stochastic SIR model ---------------------------------------------

stoch_sir <- function(ntimes, S0, I0, births, beta){
  S <- I <- numeric(ntimes)
  S[1] <- S0
  I[1] <- I0
  
  for(t in 2:ntimes){
    inv_lambda <- exp(-beta * (I[t-1] + 0.1))
    delta <- rbinom(1, S[t-1], inv_lambda)
    S[t] <- delta + births
    I[t] <- S[t-1] - delta
  }
  return(
    tibble(
      ntime = 1:ntimes,
      S = S,
      I = I
    )
  )
}


# Simulate data -----------------------------------------------------------

ntimes <- 26*10
S0 <- 10000
I0 <- 1
births <- 100
beta <- 0.0001
report_fraction <- 0.5

sim <- stoch_sir(ntimes, S0, I0, births, beta) %>%
  gather(key = state, value = abundance, S:I)

ggplot(sim, aes(x = ntime, y = abundance)) +
  geom_line() +
  facet_wrap(~state, scales = "free")


# Define the observation data ---------------------------------------------

case_data <- sim %>%
  filter(state == "I") %>%
  mutate(
    cases = abundance * report_fraction,
    year = rep(1990:(1990-1+ntimes/26), each = 26),
    biweek = rep(1:26, times = ntimes/26),
    week = 1:n(),
    births = births
  ) %>%
  dplyr::select(year, biweek, week, cases, births)


# Define pomp models ------------------------------------------------------

# Process model
sir_proc <- Csnippet(
  "
  double dN[2];
  double inv_lambda = exp(-beta_r * (nearbyint(I) + 0.1));
  double delta = rbinom(nearbyint(S), inv_lambda);
  dN[0] = delta + births;
  dN[1] = S - delta;
  S = nearbyint(dN[0]);
  I = nearbyint(dN[1]);
  "
)

# Measurement model
rmeas <- Csnippet(
  "
  cases = rnbinom_mu(1/od, nearbyint(rho * I));
  "
)

# Likelihood density model
dmeas <- Csnippet(
  "
  lik = dnbinom_mu(cases, 1/od, nearbyint(rho * I), give_log);
  "
)

# Set up parameter transformations
logtrans <- Csnippet(
  "
  Tbeta_r = log(beta_r);
  Trho = logit(rho);
  Tod = log(od);
  "
)

exptrans <- Csnippet(
  "
  Tbeta_r = exp(beta_r);
  Trho = expit(rho);
  Tod = exp(od);
  "
)

obs_data <- case_data %>%
  dplyr::select(week, cases) %>%
  mutate(cases = round(cases, 0))

birth_data <- case_data %>%
  dplyr::select(week, births) %>%
  bind_rows(
    tibble(
      week = 0,
      births = 100
    )
  ) %>%
  arrange(week)

# Generate pomp object
sir_pomp <- pomp(
  data = obs_data,
  time = "week",
  covar = birth_data,
  tcovar = "week",
  t0 = 0,
  rprocess = discrete.time.sim(sir_proc, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  statenames = c("S", "I"),
  toEstimationScale = logtrans,
  fromEstimationScale = exptrans,
  paramnames = c("beta_r", "rho", "od")
)

# Make sure the model simulates reasonable dynamics -- It does!
pomp_test <- simulate(sir_pomp, params = c(beta_r = 0.0002, rho = 0.5, od = 0.1, S.0 = 10000, I.0 = 1))
plot(pomp_test)


# Perform iterated filter to find MLEs ------------------------------------

guesses <- sobolDesign(
  lower = c(beta_r = 0.00001, rho = 0.3, od = 0.01, S.0 = 5000, I.0 = 1),
  upper = c(beta_r = 0.0005, rho = 0.7, od = 0.1, S.0 = 15000, I.0 = 20),
  nseq = 10
)

guesses$S.0 <- round(guesses$S.0, 0)
guesses$I.0 <- round(guesses$I.0, 0)

pomp_test <- simulate(sir_pomp, params = unlist(guesses[10, ]))
plot(pomp_test)

mles_out <- {}
for(i in 1:nrow(guesses)){
  mf <- mif2(
    object = sir_pomp, 
    start = unlist(guesses[i, ]),
    Nmif = 50,
    Np = 1000,
    transform = TRUE,
    cooling.fraction.50 = 0.8,
    cooling.type = "geometric",
    rw.sd = rw.sd(beta_r = 0.02, rho = 0.02, od = 0.01, I.0 = 0.02, S.0 = 0.1)
  ) %>%
    mif2()
  
  ll <- logmeanexp(replicate(10, logLik(pfilter(mf))), se=TRUE)
  mles <- data.frame(loglik=ll[1],loglik.se=ll[2],as.list(coef(mf)))
  mles_out <- rbind(mles_out, mles)
  
  cat(paste("Done with mif run", i, "of", nrow(guesses)))
  cat("\n")
}

test <- pfilter(sir_pomp, params = unlist(guesses[1,]), Np = 1000)
test <- pfilter(sir_pomp, params = c(beta_r = 1.066219e-04, rho = 4.664018e-01, od = 0.01, S.0 = 1.187500e+04, I.0 = 12), Np = 1000)
plot(test)
