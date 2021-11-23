## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, message = FALSE, warning = FALSE)

library(tidyverse)   # data wrangling
library(lubridate)   # time and date functions
library(pomp)        # fitting state-space models with particle filtering
library(foreach)     # functions for parallel computing
library(doParallel)  # functions for parallel computing

registerDoParallel()
theme_set(theme_minimal())

## ----load-data-----------------------------------------------------------




## ----extras, eval = FALSE, echo=FALSE------------------------------------
## // double beta_r = dot_product(K, &xi1, &b1);
##  //inv_lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
## // cases = rnorm(rho*I, sqrt( pow(tau*I,2) + rho*I ) );
##  // double tol = 1.0e-25;
##  // double sd_cases = sqrt(pow(tau*I, 2) + mean_cases);
##  // if (cases > 0.0) {
##  //  lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol;
##  // } else{
##  //  lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
##  // }
##  // if (give_log) lik = exp(lik);
## 
## # Define parameter search space based on probe.match results
## # param_box <- rbind(
## #   beta_mu = c(0.0000001, 0.0001),
## #   b1 = c(-5, 5),
## #   b2 = c(-5, 5),
## #   b3 = c(-5, 5),
## #   b4 = c(-5, 5),
## #   b5 = c(-5, 5),
## #   b6 = c(-5, 5),
## #   rho = c(0.3, 0.7),
## #   tau = c(0.00001, 0.2),
## #   S_0 = c(10000, 50000),
## #   I_0 = c(5, 100),
## #   psi = c(0, 10)
## # )
## 
## 
## # bake(file="mif1.rds",seed=900242057,kind="L'Ecuyer",{
## #   foreach (i=1:20,.packages='pomp',.combine=c,.options.multicore=list(set.seed=TRUE)) %dopar% {
## #     guess <- make_guess(probe_ests, sd_perc = 0.2)
## #     mf <- mif2(measles_pomp,
## #          start=c(guess),
## #          Np=1000,
## #          Nmif=50,
## #          cooling.type="geometric",
## #          cooling.fraction.50=1,
## #          transform=TRUE,
## #          rw.sd=rw.sd(
## #            beta_mu=0.02, rho=0.02, tau=0.02, sigma_env = 0.02,
## #            b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
## #            I_0=ivp(0.1), S_0=ivp(0.1), psi = 0.02
## #          )
## #         )
## #     ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
## #     data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
## #   }
## # }) -> mifs1
## 
## # mifs1 %>%
## #   conv.rec(c("loglik", "nfail", names(params))) %>%
## #   reshape2::melt() %>%
## #   ggplot(aes(x=iteration,y=value,color=as.factor(L1),group=L1))+
## #   geom_line()+
## #   guides(color=FALSE)+
## #   labs(x="MIF2 Iteration",y="")+
## #   facet_wrap(~variable,scales="free_y",ncol=2)+
## #   theme_bw()
## #
## # iter1_coef <- mifs1 %>%
## #   conv.rec(c("loglik", names(params))) %>%
## #   reshape2::melt() %>%
## #   filter(iteration == 49) %>%
## #   spread(key = variable, value = value) %>%
## #   filter(loglik == max(loglik, na.rm = TRUE)) %>%
## #   dplyr::select(-iteration, -L1, -loglik) %>%
## #   unlist()
## #
## # bake(file="mif2.rds",seed=902057,kind="L'Ecuyer",{
## #   foreach (i=1:20,.packages='pomp',.combine=c,.options.multicore=list(set.seed=TRUE)) %dopar% {
## #     guess <- make_guess(iter1_coef, sd_perc = 0.2)
## #     mif2(measles_pomp,
## #          start=c(guess),
## #          Np=1000,
## #          Nmif=50,
## #          cooling.type="geometric",
## #          cooling.fraction.50=0.5,
## #          transform=TRUE,
## #          rw.sd=rw.sd(
## #            beta_mu=0.02, rho=0.02, tau=0.02, sigma_env = 0.02,
## #            b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
## #            I_0=ivp(0.1), S_0=ivp(0.1), psi = 0.02
## #          )
## #         )
## #   }
## # }) -> mifs2
## #
## # mifs2 %>%
## #   conv.rec(c("loglik", "nfail", names(params))) %>%
## #   reshape2::melt() %>%
## #   ggplot(aes(x=iteration,y=value,color=as.factor(L1),group=L1))+
## #   geom_line()+
## #   guides(color=FALSE)+
## #   labs(x="MIF2 Iteration",y="")+
## #   facet_wrap(~variable,scales="free_y",ncol=2)+
## #   theme_bw()
## #
## # iter2_coef <- mifs2 %>%
## #   conv.rec(c("loglik", names(params))) %>%
## #   reshape2::melt() %>%
## #   filter(iteration == 49) %>%
## #   spread(key = variable, value = value) %>%
## #   filter(loglik == max(loglik, na.rm = TRUE)) %>%
## #   dplyr::select(-iteration, -L1, -loglik) %>%
## #   unlist()
## 
## 
## # # Simulate from initial conditions
## # model_sims <- simulate(
## #   measles_pomp,
## #   params = iter2_coef,
## #   nsim = 1000,
## #   as.data.frame = TRUE,
## #   include.data = TRUE) %>%
## #   as_tibble()
## #
## # model_data <- model_sims %>%
## #   filter(sim == "data")
## #
## # simulated_data <- model_sims %>%
## #   filter(sim != "data") %>%
## #   group_by(time) %>%
## #   summarise(
## #     med_estimate = median(cases),
## #     upper_estimate = quantile(cases, 0.95),
## #     lower_estimate = quantile(cases, 0.05)
## #   )
## #
## # ggplot(simulated_data, aes(x = time)) +
## #   geom_ribbon(aes(ymin = lower_estimate, ymax = upper_estimate), alpha = 0.2, fill = "blue") +
## #   geom_line(aes(y = med_estimate), color = "blue") +
## #   geom_line(data = model_data, aes(x = time, y = cases))
## #
## # ggplot(filter(model_sims, sim != "data")) +
## #   geom_line(aes(x = time, y = S, group= sim), alpha = 0.2)
## 

## ----process-model-------------------------------------------------------
measles_process <- Csnippet(
  "
  double dN[2];

  double beta_r = beta_mu * (1 + dot_product(K, &xi1, &b1));
  double var_epsilon = sigma_env;
  beta_r *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon, var_epsilon);
  double inv_lambda = exp(-beta_r * (I + psi));
 
  if(inv_lambda > 1){
    inv_lambda = 1;
  }
  double delta = rbinom(nearbyint(S), inv_lambda);
  dN[0] = delta + births;
  dN[1] = S - delta;
  S = nearbyint(dN[0]);
  I = nearbyint(dN[1]);
  beta_t = beta_r;
  "
)

## ----measurement-model---------------------------------------------------
rmeas <- Csnippet(
  "
  cases = rnbinom_mu(1/tau, rho*I);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
  "
)

## ----likelihood-model----------------------------------------------------
dmeas <- Csnippet(
  "
  double mean_cases = nearbyint(rho*I);
  lik = dnbinom_mu(cases, 1/tau, mean_cases, give_log);
  "
)

## ----mapping-------------------------------------------------------------
# Map initial value parameters to initial values of states
init <- Csnippet(
  "
  I = nearbyint(I_0);
  S = nearbyint(S_0);
  beta_t = 0;
  "
)

# Set up parameter transformations
to_est <- Csnippet(
  "
  Tbeta_mu = log(beta_mu);
  Trho = logit(rho);
  Ttau = log(tau);
  Tsigma_env = log(sigma_env);
  TS_0 = log(S_0);
  TI_0 = log(I_0);
  Tpsi = log(psi);
  "
)

from_est <- Csnippet(
  "
  Tbeta_mu = exp(beta_mu);
  Trho = expit(rho);
  Ttau = exp(tau);
  Tsigma_env = exp(sigma_env);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  Tpsi = exp(psi);
  "
)

## ----data-and-covar------------------------------------------------------
# Fetch observations
obs_data <- biweek_data %>%
  dplyr::select(biweek, cases) %>%
  mutate(cases = round(cases, 0))

# Generate basis functions for seasonality
bspline_basis <- periodic.bspline.basis(
  obs_data$biweek,
  nbasis = 6,
  degree = 3,
  period = 26,
  names = "xi%d"
) %>%
  as_tibble()

# Combine with births data for covariate table
covar_data <- biweek_data %>%
  dplyr::select(biweek, births) %>%
  arrange(biweek) %>%
  bind_cols(bspline_basis)

## ----pomp-container------------------------------------------------------
# Set some realistic parameter values for testing via simulation
params <- c(
  beta_mu = 0.00001,
  b1 = 3,
  b2 = 0,
  b3 = 1.5,
  b4 = 6,
  b5 = 5,
  b6 = 3,
  psi = 1,
  rho = 0.1,
  tau = 1,
  sigma_env = 0.01,
  S_0 = 100000, 
  I_0 = 100
)

# Generate pomp object
measles_pomp <- pomp(
  data = obs_data,
  time = "biweek",
  covar = covar_data,
  tcovar = "biweek",
  t0 = 1,
  rprocess = discrete.time.sim(measles_process, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  initializer = init,
  statenames = c("S", "I", "beta_t"),
  toEstimationScale = to_est,
  fromEstimationScale = from_est,
  paramnames = names(params),
  params = params,
  globals = "int K = 6;"
)

## ----simulation-test-----------------------------------------------------
simulate(
  measles_pomp, 
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
    geom_line() +
    scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
    guides(color = FALSE) +
    facet_wrap(~sim, ncol = 2) +
    scale_y_sqrt() +
    theme(strip.text=element_blank())

## ----probe---------------------------------------------------------------
probe.zeroes <- function(y){
  xy <- y["cases", ]
  as.numeric(length(which(xy == 0)))
}
probe.max <- function(y){
  max(y["cases", ], na.rm = TRUE)
}
probe.cumsum <- function(y){
  cases <- y["cases", ]
  max(cumsum(cases))
}

# probe(measles_pomp, probes = list(probe.zeroes), nsim = 500)

plist <- list(
  probe.zeroes,
  probe.max,
  probe.cumsum,
  probe.acf("cases", lags = 1,transform = sqrt, type = "correlation")
)


stew(
  file = "probe_match_results.rda", {
    pm <- probe.match(
      measles_pomp,
      probes = plist,
      est = names(params),
      nsim = 500,
      transform = TRUE,
      start = params,
      method = "Nelder-Mead",
      maxit = 10000
    )
  },
  seed = 290808325,
  kind = "L'Ecuyer"
)

pm_ll <- logLik(pfilter(measles_pomp, Np = 1000, params = coef(pm)))

simulate(
  measles_pomp, 
  params = coef(pm),
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
    geom_line() +
    scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
    guides(color = FALSE) +
    facet_wrap(~sim, ncol = 2) +
    scale_y_sqrt() +
    theme(strip.text=element_blank())

## ----global-mif----------------------------------------------------------
probe_ests <- coef(pm)

make_guess <- function(params, sd_perc = 0.1){
  n <- length(params)
  out <- numeric(n)
  for(i in 1:n){
    out[i] <- rnorm(1, params[i], abs(params[i])*sd_perc)
  }
  names(out) <- names(params)
  return(out)
}

stew(file="mif1.rda",{
  w1 <- getDoParWorkers()
  t1 <- system.time({
    m1 <- foreach(i=1:20,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- make_guess(probe_ests, sd_perc = 0.2)
      mf <- mif2(measles_pomp,
                 start=c(guess),
                 Np=1000,
                 Nmif=50,
                 cooling.type="geometric",
                 cooling.fraction.50=1,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   beta_mu=0.02, rho=0.02, tau=0.01, sigma_env = 0.02,
                   b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
                   I_0=ivp(0.1), S_0=ivp(0.1), psi = 0.02
                 )
                )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
},seed=290860873,kind="L'Ecuyer")

m1_ests <- m1 %>%
  as_tibble() %>%
  arrange(-loglik) %>%
  mutate(
    rownum = 1:n()
  ) %>%
  filter(rownum <= 10) %>%
  select(-loglik, -loglik.se, -rownum)



stew(file="mif2.rda",{
  w2 <- getDoParWorkers()
  t2 <- system.time({
    m2 <- foreach(i=1:10,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- make_guess(unlist(m1_ests[i,]), sd_perc = 0.2)
      mf <- mif2(measles_pomp,
                 start=c(guess),
                 Np=1000,
                 Nmif=50,
                 cooling.type="geometric",
                 cooling.fraction.50=0.5,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   beta_mu=0.02, rho=0.02, tau=0.01, sigma_env = 0.02,
                   b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
                   I_0=ivp(0.1), S_0=ivp(0.1), psi = 0.02
                 )
                )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
},seed=2908873,kind="L'Ecuyer")



iter2_coef <- m2 %>%
  as_tibble() %>%
  filter(loglik == max(loglik)) %>%
  select(-loglik, -loglik.se) %>%
  unlist()

(ll_end <- logmeanexp(replicate(50, logLik(pfilter(measles_pomp, params = iter2_coef, Np=5000))),se=TRUE))

simulate(
  measles_pomp, 
  params = iter2_coef,
  nsim = 9,
  as.data.frame = TRUE,
  include.data = TRUE) %>%
  ggplot(aes(x = time, y = cases, group = sim, color = (sim == "data"))) +
    geom_line() +
    scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "red"))+
    guides(color = FALSE) +
    facet_wrap(~sim, ncol = 2) +
    scale_y_sqrt() +
    theme(strip.text=element_blank())

bs <- iter2_coef[c("b1","b2","b3","b4","b5","b6")]
beta_season <- (as.matrix(bspline_basis[1:26,]) %*% bs) * iter2_coef["beta_mu"]
betas <- tibble(
  beta_value = as.numeric(beta_season),
  biweek = 1:26
)

ggplot(betas, aes(x = biweek, y = beta_value)) +
  geom_line(color = "blue") +
  geom_point(color = "white", size = 4) +
  geom_point(color = "blue", size = 2) +
  ylab(expression(beta)) +
  xlab("Biweek") +
  ggtitle("Estimate seasonality in transmission")

## ----pmcmc---------------------------------------------------------------
hyperparams <- list(min = iter2_coef/10, max = iter2_coef * 10)

measles_dprior <- function (params, ..., log) {
 f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max, log = TRUE))
 if (log) f else exp(f)
}



stew(
  file = "pmcmc_results.rda",
  {
    chain <- pmcmc(
      pomp(measles_pomp, dprior = measles_dprior), 
      Nmcmc=2000, 
      Np=100,
      start= iter2_coef,
      max.fail=Inf,
      proposal = mvn.diag.rw(
        rw.sd = c(
          beta_mu=0.00001,
          b1 = 1,
          b2 = 1,
          b3 = 1,
          b4 = 1,
          b5 = 1,
          b6 = 1,
          S_0 = 1000,
          psi = 0.1,
          sigma_env = 0.1,
          tau = 0.001,
          I_0 =10,
          rho =0.01
        )
      )
    )
    
  mcmcout <- chain %>% pmcmc(Nmcmc=10000,proposal=mvn.rw(covmat(chain))) 
  },
  seed=29012387, kind="L'Ecuyer"
)

traces <- conv.rec(mcmcout)
traces <- window(traces,thin=1,start=1)
plot(traces[,"S_0"])
hist(traces[,"beta_mu"], col = "grey", border = "white")
mle_rho = mean(traces[,"rho"])

test <- mcmcout %>% 
  filter.traj() %>%
  reshape2::melt() %>%
  as_tibble() %>%
  filter(variable != "beta_t") %>%
  # filter(rep > 1000) %>%
  dplyr::group_by(variable, time) %>%
  dplyr::summarise(
    median_value = quantile(value, 0.5),
    upper_value = quantile(value, 0.975),
    lower_value = quantile(value, 0.025)
  ) 

ggplot(test, aes(x = time)) +
  geom_ribbon(aes(ymin = lower_value, ymax = upper_value), alpha = 0.3) +
  geom_line(aes(y = median_value)) +
  facet_wrap(~variable, scales = "free", ncol = 1) +
  ggtitle("Estimates of latent states")

pred_obs <- test %>%
  filter(variable == "I") %>%
  mutate(
    median_value = median_value*mle_rho,
    upper_value = upper_value*mle_rho,
    lower_value = lower_value*mle_rho,
    cases = biweek_data$cases,
    date = biweek_data$date
  )

ggplot(pred_obs, aes(x = date)) +
  geom_point(aes(y = cases), color = "grey35") +
  geom_ribbon(aes(ymin = lower_value, ymax = upper_value), fill = "coral", alpha = 0.4) +
  geom_line(aes(y = median_value), color = "coral") +
  labs(x = "Date", y = "Reported cases") +
  ggtitle("Estimates of infectious state and data")

popsmooth <- predict(lm(biweek_data$population~biweek_data$biweek), x=biweek_data$biweek)

reff_time <- mcmcout %>% 
  filter.traj() %>%
  reshape2::melt() %>%
  as_tibble() %>%
  filter(variable == "beta_t") %>%
  mutate(
    date = rep(biweek_data$date, each = length(unique(rep))),
    population = rep(popsmooth, each = length(unique(rep)))
  ) %>%
  dplyr::group_by(variable, date) %>%
  dplyr::summarise(
    median_value = median(value*population),
    upper_value = quantile(value*population, 0.975),
    lower_value = quantile(value*population, 0.025)
  )

ggplot(reff_time, aes(x = week(date), y = median_value)) +
  geom_ribbon(aes(ymin = lower_value, ymax = upper_value), alpha = 0.2) +
  geom_line(size = 0.3) +
  geom_point(size = 0.7) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = 2) +
  facet_wrap(~year(date)) +
  xlab("Week") +
  ylab(expression(R[t])) +
  ggtitle("Effective reproduction number over time")

