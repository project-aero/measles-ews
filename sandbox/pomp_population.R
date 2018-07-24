# https://kingaa.github.io/pomp/vignettes/getting_started.html#specifying-a-discrete-time-process-model-and-skeleton

library(tidyverse)
library(pomp)

parus_dat <- read.csv(text="
                      year,P
                      1960,148
                      1961,258
                      1962,185
                      1963,170
                      1964,267
                      1965,239
                      1966,196
                      1967,132
                      1968,167
                      1969,186
                      1970,128
                      1971,227
                      1972,174
                      1973,177
                      1974,137
                      1975,172
                      1976,119
                      1977,226
                      1978,166
                      1979,161
                      1980,199
                      1981,306
                      1982,206
                      1983,350
                      1984,214
                      1985,175
                      1986,211"
)

ggplot(parus_dat, aes(x = year, y = P)) +
  geom_line()

# Define the process model
bev_holt_proc <- Csnippet(
  "
  double eps = rlnorm(-sigma*sigma/2,sigma);
  N = a*N/(1+b*N)*eps;
  "
)

# Define measurement model
rmeas <- Csnippet(
  "
  P = rpois(N);
  "
)

# Define likelihood density
dmeas <- Csnippet(
  "
  lik = dpois(P,N,give_log);
  "
)

# Set up parameter transformations
logtrans <- Csnippet(
  "
  Ta = log(a);
  Tb = log(b);
  Tsigma = log(sigma);
  "
)

exptrans <- Csnippet(
  "
  Ta = exp(a);
  Tb = exp(b);
  Tsigma = exp(sigma);
  "
)

# Generate pomp object
parus_bh <- pomp(
  data = parus_dat,
  time = "year",
  t0 = 1959,
  rprocess = discrete.time.sim(bev_holt_proc, delta.t = 1),
  rmeasure = rmeas,
  dmeasure = dmeas,
  statenames = "N",
  toEstimationScale = logtrans,
  fromEstimationScale = exptrans,
  paramnames = c("a", "b", "sigma")
)

# Perform iterated filtering
guesses <- sobolDesign(lower=c(a=0.1,b=0.00001,sigma=0,N.0=100),
                       upper=c(a=20,b=0.01,sigma=5,N.0=400),
                       nseq=10)

mles_out <- {}
for(i in 1:nrow(guesses)){
  mf <- mif2(object = parus_bh, start=c(a = guesses$a[i], b = guesses$b[i], sigma = guesses$sigma[i], N.0 = guesses$N.0[i]),
             Nmif=50,Np=1000,transform=TRUE,
             cooling.fraction.50=0.8,cooling.type="geometric",
             rw.sd=rw.sd(a=0.02,b=0.02,sigma=0.2))
  ll <- logmeanexp(replicate(5,logLik(pfilter(mf))),se=TRUE)
  mles <- data.frame(loglik=ll[1],loglik.se=ll[2],as.list(coef(mf)))
  mles_out <- rbind(mles_out, mles)
}

# PMCMC
parus_bh %<>%
  pomp(dprior=Csnippet("
    lik = dunif(a,0,50,1)+dunif(b,0,1,1)+dunif(sigma,0,2,1);
    lik = (give_log) ? lik : exp(lik);
  "),paramnames=c("a","b","sigma"))


starts <- mles_out %>%
  filter(loglik == max(loglik)) %>%
  filter(a > 0 & a < 50 & b > 0 & b < 1) %>%
  dplyr::select(-loglik, -loglik.se)

chain <- parus_bh %>%
  pmcmc(Nmcmc=2000,Np=200,start=starts,
        proposal=mvn.rw.adaptive(rw.sd=c(a=0.01,b=0.01,sigma=0.01),
                                 scale.start=100,shape.start=100))
        
chain_final <- pmcmc(chain, Nmcmc=10000, proposal=mvn.rw(covmat(chain)))

posts <- chain_final %>% 
  filter.traj() %>% 
  reshape2::melt() %>%
  filter(rep > 1000 & rep %% 50 == 0) %>%
  group_by(time) %>%
  summarise(
    median_N = quantile(value, 0.5),
    upper_N = quantile(value, 0.975),
    lower_N = quantile(value, 0.025)
  ) %>%
  filter(time > min(time)) %>%
  mutate(
    observation = parus_dat$P
  )

ggplot(posts, aes(x = time)) +
  geom_line(aes(y = median_N)) +
  geom_line(aes(y = upper_N), linetype = 2) +
  geom_line(aes(y = lower_N), linetype = 2) +
  geom_point(aes(y = observation), size = 2) +
  theme_minimal()
