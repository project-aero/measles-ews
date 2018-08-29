# Example of fitting a stochastic SIR model with seasonality using `pomp`


# Load packages -----------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
registerDoParallel()


# Define stochastic SIR model ---------------------------------------------

stoch_sir <- function(ntimes, S0, I0, births, beta, amp, phase = 1){
  S <- I <- numeric(ntimes)
  S[1] <- S0
  I[1] <- I0
  
  for(t in 2:ntimes){
    betanow <- beta * (1 + (amp * sin((2*pi*t)/26 + phase)))
    inv_lambda <- exp(-betanow * (I[t-1] + 0.1))
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
births <- 150
beta <- 0.00015
amp <- 0.2
phase <- 1
report_fraction <- 0.5

set.seed(72826491)
sim <- stoch_sir(ntimes, S0, I0, births, beta, amp, phase) %>%
  gather(key = state, value = abundance, S:I)

ggplot(sim, aes(x = ntime, y = abundance)) +
  geom_line() +
  facet_wrap(~state, scales = "free") +
  scale_x_continuous(breaks = seq(1, 26*10, by = 26))


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
  double beta_r = beta_mu * (1 + dot_product(K, &xi1, &b1));
  double inv_lambda = exp(-beta_r * (I + 0.1));
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
  cases = rnorm(rho*I, sqrt( pow(tau*I,2) + rho*I ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
  
  //cases = rnbinom_mu(1/od, nearbyint(rho * I));
  "
)

# Likelihood density model
dmeas <- Csnippet(
  "
  double tol = 1.0e-25;
  double mean_cases = rho*I;
  double sd_cases = sqrt(pow(tau*I, 2) + mean_cases);
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);

  //lik = dnbinom_mu(cases, 1/od, nearbyint(rho * I), give_log);
  "
)

# Map initial value parameters to initial values of states
init <- Csnippet(
  "
  I = nearbyint(I_0);
  S = nearbyint(S_0);
  "
)

# Set up parameter transformations
to_est <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Trho = logit(rho);
  Ttau = log(tau);
  TS_0 = log(S_0);
  TI_0 = log(I_0);
  "
)

from_est <- Csnippet(
  "
  Tbeta_mu = exp(beta_mu);
  Trho = expit(rho);
  Ttau = exp(tau);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  "
)

obs_data <- case_data %>%
  dplyr::select(week, cases) %>%
  mutate(cases = round(cases, 0))

bspline_basis <- periodic.bspline.basis(
  obs_data$week,
  nbasis = 6,
  degree = 3,
  period = 26,
  names = "xi%d"
) %>%
  as_tibble()

covar_data <- case_data %>%
  dplyr::select(week, births) %>%
  arrange(week) %>%
  bind_cols(bspline_basis)

# Set some realistic parameter values for testing via simulation
params <- c(beta_mu = 0.0001, b1=3, b2=0, b3=0, b4=0, b5=0, b6=0, rho = 0.45, tau = 0.1, S_0 = 9000, I_0 = 1)

# Generate pomp object
sir_pomp <- pomp(
  data = obs_data,
  time = "week",
  covar = covar_data,
  tcovar = "week",
  t0 = 0,
  rprocess = discrete.time.sim(sir_proc, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  initializer = init,
  statenames = c("S", "I"),
  toEstimationScale = to_est,
  fromEstimationScale = from_est,
  paramnames = c("beta_mu", "b1", "b2", "b3", "b4", "b5", "b6", "rho", "tau", "S_0", "I_0"),
  params = params,
  globals = "int K = 6;"
)

# Make sure the model simulates reasonable dynamics -- It does!
nsim <- 9
test_sim <- simulate(sir_pomp, nsim=nsim,as.data.frame=TRUE,include.data=TRUE)
ggplot(data = test_sim, aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())

# Make sure the likelihood function is working -- It is!
pf <- pfilter(sir_pomp, Np = 1000)
logLik(pf)

set.seed(493536993,kind="L'Ecuyer")
t1 <- system.time(
  pf1 <- foreach(i=1:10,.packages='pomp',
                 .options.multicore=list(set.seed=TRUE)
  ) %dopar% {
    pfilter(sir_pomp, Np=5000)
  }
)
(L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))

stew(file="local_search_season.rda",{
  w1 <- getDoParWorkers()
  t1 <- system.time({
    m1 <- foreach(i=1:10,
                  .packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      mf <- mif2(sir_pomp,
                 Np=1000,
                 Nmif=50,
                 cooling.type="geometric",
                 cooling.fraction.50=1,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   beta_mu=0.02, rho=0.02, tau=0.02,
                   b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
                   I_0=ivp(0.2), S_0=ivp(0.2)
                 )
      )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  }
  )
},seed=318817883,kind="L'Ecuyer")
  
pairs(~loglik+beta_mu+rho+tau+I_0+S_0,data=m1)
# write.csv(m1, file="sir_params.csv",row.names=FALSE,na="")

mles_local <- m1 %>%
  filter(loglik == max(loglik)) %>%
  dplyr::select(-loglik, -loglik.se) %>%
  unlist()

fitted_sim <- simulate(sir_pomp, nsim=9, as.data.frame=TRUE, include.data=TRUE, params = mles_local)
ggplot(data = fitted_sim, aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1))+
  guides(color = FALSE) +
  facet_wrap(~sim, ncol = 2) +
  scale_y_sqrt() +
  theme(strip.text=element_blank())

mf <- mif2(sir_pomp,
           start=mles_local,
           Np=2000,
           Nmif=300,
           cooling.type="geometric",
           cooling.fraction.50=0.5,
           transform=TRUE,
           rw.sd=rw.sd(
             beta_mu=0.02, rho=0.02, tau=0.02,
             b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
             I_0=ivp(0.2), S_0=ivp(0.2)
           )
          )

# Check to see if b-splines are picking up correct seasonality
betasp <- coef(mf)[c("b1","b2","b3","b4","b5","b6")]
beta_est <- coef(mf)["beta_mu"]
plot(as.matrix(bspline_basis[1:26,]) %*% betasp * beta_est,
     ylim = c(0.00009, 0.0002), col = "blue", xlab = "Biweek", 
     ylab = expression(beta))
lines(beta * (1 + (amp * sin((2*pi*1:26)/26 + phase))), col = "red")


# Global MLE search -------------------------------------------------------

param_box <- rbind(
  beta_mu = c(0.00001, 0.001),
  b1 = c(0, 10),
  b2 = c(0, 10),
  b3 = c(0, 10),
  b4 = c(0, 10),
  b5 = c(0, 10),
  b6 = c(0, 10),
  rho = c(0.3, 0.8),
  tau = c(0.00001, 0.2),
  S_0 = c(5000, 15000),
  I_0 = c(0, 10)
)

stew(file="global_search_seasonal.rda",{
  w2 <- getDoParWorkers()
  t2 <- system.time({
    m2 <- foreach(i=1:10,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- apply(param_box,1,function(x)runif(1,x[1],x[2]))
      mf <- mif2(sir_pomp,
                 start=c(guess),
                 Np=2000,
                 Nmif=300,
                 cooling.type="geometric",
                 cooling.fraction.50=0.5,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   beta_mu=0.02, rho=0.02, tau=0.02,
                   b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
                   I_0=ivp(0.2), S_0=ivp(0.2)
                 )
                )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
},seed=290860873,kind="L'Ecuyer")

mles <- unlist(
  filter(m2, loglik == max(loglik)) %>%
    dplyr::select(-loglik, -loglik.se)
  )

fitted_sim <- simulate(sir_pomp, nsim=100, as.data.frame=TRUE, include.data=TRUE, params = mles)
ggplot(data = fitted_sim, aes(x = time, y = cases, group = sim, color = (sim == "data"), alpha = (sim == "data"))) +
  geom_line() +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "orange"))+
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1))+
  guides(color = FALSE, alpha = FALSE) +
  scale_y_sqrt()


# Perform particle MCMC for final inference -------------------------------

stew(
  file = "sir_mcmc_seasonal.rda",
  {
    outmcmc <- pmcmc(sir_pomp, Nmcmc=10000,Np=200,start=mles,
                        proposal=mvn.rw.adaptive(rw.sd=c(beta_mu=0.00001,b1=0.01, b2 = 0.01, b3 = 0.01, b4 = 0.01, b5 = 0.01, b6 = 0.01,
                                                         S_0=100,tau=0.01,I_0 = 0.01, rho = 0.1),
                                                 scale.start=10,shape.start=10))
  },
  seed=29086087,kind="L'Ecuyer"
)

traces <- conv.rec(outmcmc)
traces <- window(traces,thin=1,start=1000)
plot(traces[,"rho"])
hist(traces[,"beta_mu"], col = "grey", border = "white")
abline(v = beta, col = "red", lwd = 2)
mle_rho = mean(traces[,"rho"])

test <- outmcmc %>% 
  filter.traj() %>%
  reshape2::melt() %>%
  as_tibble() %>%
  filter(rep > 1000) %>%
  dplyr::group_by(variable, time) %>%
  dplyr::summarise(
    median_value = quantile(value, 0.5),
    upper_value = quantile(value, 0.975),
    lower_value = quantile(value, 0.025)
  ) %>%
  filter(time > 0) %>%
  mutate(
    observation = filter(fitted_sim, sim == "data") %>% pull(cases)
  )

ggplot(test, aes(x = time)) +
  geom_ribbon(aes(ymin = lower_value, ymax = upper_value), alpha = 0.3) +
  geom_line(aes(y = median_value)) +
  facet_wrap(~variable, scales = "free", ncol = 1) +
  scale_y_sqrt()

ggplot(filter(test, variable == "I"), aes(x = time)) +
  geom_ribbon(aes(ymin = lower_value, ymax = upper_value), alpha = 0.3, fill = "coral") +
  geom_line(aes(y = median_value), color = "coral") +
  geom_point(aes(y = observation/mle_rho), color = "grey35", size = 1, alpha = 0.6) +
  scale_y_sqrt() 
  