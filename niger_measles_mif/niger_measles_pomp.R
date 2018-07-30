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
file_name <- "../niger_measles/niger_regional_1995_2005.csv"
niger_measles_raw <- read_csv(file_name, col_types = cols())

num_regions <- nrow(niger_measles_raw)
num_weeks <- ncol(niger_measles_raw) - 1  # subtract 1 from ncol() because first column are regions
weeks_per_year <- num_weeks/11
weeks <- rep(1:52, times = 11)

# Create a vector of years for all num_weeks
years <- rep(1995:2005, each = weeks_per_year)

# Function for calculating start of week based on week number and year
calculate_start_of_week = function(week, year) {
  date <- ymd(paste(year, 1, 1, sep="-"))
  week(date) = week
  return(date)
}

# Clean up the data frame
measles_data <- niger_measles_raw %>%
  gather(key = week, value = cases, -X1) %>%
  mutate(
    week_num = rep(1:num_weeks, each = num_regions),
    year = rep(years, each  = num_regions),
    week = rep(weeks, each = num_regions),
    date = calculate_start_of_week(week, year)
  ) %>%
  dplyr::rename(region = X1) %>%
  filter(region == "Niamey (City)")

# Read in population data and estimate births
population <- read_csv("../niger_measles/district_pops.csv") %>%
  gather(key = year, value = population, -X1) %>%
  rename(district = X1) %>%
  filter(district == "Niamey I")

birth_rates <- read_csv("../niger_measles/niger_crude_birth_rates.csv") %>%
  mutate(
    date = mdy(date),  # lubridate prefixes any 2digit year 00-68 with 20, not a problem for us though
    year = as.character(year(date)),
    rate_per_person = births_per_thousand/1000
  ) %>%
  select(year, rate_per_person)

births <- population %>%
  left_join(birth_rates, by = "year") %>%
  mutate(
    births_per_year = population * rate_per_person,
    births_per_week = births_per_year / 52
  )

# Aggregate weekly data to biweeks
biweek_data <- measles_data %>%
  mutate(biweek = rep(1:(n()/2), each = 2)) %>%
  group_by(biweek) %>%
  summarise(
    cases = sum(cases),
    date = min(date)
  ) %>%
  mutate(
    births = round(rep(births$births_per_week*2, each = 26)),
    population = round(rep(population$population, each = 26))
  )

ggplot(biweek_data, aes(x = date, y = cases)) +
  geom_point() +
  geom_line() +
  scale_y_sqrt() +
  labs(x = "Date", y = expression(paste(sqrt(Reported~cases))))

## ----process-model-------------------------------------------------------
measles_process <- Csnippet(
  "
  double dN[2];
  double beta_r = beta_mu * (1 + dot_product(K, &xi1, &b1));
  double inv_lambda = exp(-beta_r * (I + psi));
  double delta = rbinom(nearbyint(S), inv_lambda);
  dN[0] = delta + births;
  dN[1] = S - delta;
  S = nearbyint(dN[0]);
  I = nearbyint(dN[1]);
  "
)

## ----measurement-model---------------------------------------------------
rmeas <- Csnippet(
  "
  cases = rnorm(rho*I, sqrt( pow(tau*I,2) + rho*I ) );
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
  double tol = 1.0e-25;
  double mean_cases = rho*I;
  double sd_cases = sqrt(pow(tau*I, 2) + mean_cases);
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
  "
)

## ----mapping-------------------------------------------------------------
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
  beta_mu = 0.000003,
  b1 = 3,
  b2 = 0,
  b3 = 1.5,
  b4 = 6,
  b5 = 5,
  b6 = 3,
  psi = 1,
  rho = 0.5,
  tau = 0.1,
  S_0 = 20000, 
  I_0 = 20
)

# Generate pomp object
measles_pomp <- pomp(
  data = obs_data,
  time = "biweek",
  covar = covar_data,
  tcovar = "biweek",
  t0 = 0,
  rprocess = discrete.time.sim(measles_process, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  initializer = init,
  statenames = c("S", "I"),
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

## ----global-mif----------------------------------------------------------
# Define parameter search space
param_box <- rbind(
  beta_mu = c(0.0000001, 0.0001),
  b1 = c(-5, 5),
  b2 = c(-5, 5),
  b3 = c(-5, 5),
  b4 = c(-5, 5),
  b5 = c(-5, 5),
  b6 = c(-5, 5),
  rho = c(0.3, 0.7),
  tau = c(0.00001, 0.2),
  S_0 = c(10000, 50000),
  I_0 = c(5, 100),
  psi = c(0, 10)
)

stew(file="global_search_measles.rda",{
  w1 <- getDoParWorkers()
  t1 <- system.time({
    m1 <- foreach(i=1:2,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- apply(param_box,1,function(x)runif(1,x[1],x[2]))
      mf <- mif2(measles_pomp,
                 start=c(guess),
                 Np=2000,
                 Nmif=300,
                 cooling.type="geometric",
                 cooling.fraction.50=0.5,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   beta_mu=0.02, rho=0.02, tau=0.02,
                   b1 = 0.02, b2 = 0.02, b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02,
                   I_0=ivp(0.2), S_0=ivp(0.2), psi = 0.02
                 )
                )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
},seed=290860873,kind="L'Ecuyer")



