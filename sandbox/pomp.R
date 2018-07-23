# pomp_test.R

library(pomp)



# Simple SIR --------------------------------------------------------------

rmeas <- "cases = rnbinom_mu(theta, rho * H);"
dmeas <- "lik = dnbinom_mu(cases, theta, rho * H, give_log);"

sir.step <- "
  double rate[6];
  double dN[6];
  double P;
  P = S + I + R;
  rate[0] = mu * P;       // birth
  rate[1] = beta * I / P; // transmission
  rate[2] = mu;           // death from S
  rate[3] = gamma;        // recovery
  rate[4] = mu;           // death from I
  rate[5] = mu;           // death from R
  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  H += dN[1];
"

pomp(data = data.frame(cases = NA, time = seq(0, 10, by=1/52)),
     times = "time", t0 = -1/52, dmeasure = Csnippet(dmeas),
     rmeasure = Csnippet(rmeas), rprocess = euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/52/20),
     statenames = c("S", "I", "R", "H"),
     paramnames = c("gamma", "mu", "theta", "beta", "popsize","rho", "S.0", "I.0", "R.0"), zeronames=c("H"),
     initializer=function(params, t0, ...) {
       fracs <- params[c("S.0", "I.0", "R.0")]
       setNames(c(round(params["popsize"]*fracs/sum(fracs)),0),c("S","I","R","H"))
       }, 
     params = c(popsize = 500000, beta = 400, gamma = 26,
                mu = 1/50, rho = 0.1, theta = 100, S.0 = 26/400,
                I.0 = 0.002, R.0 = 1)) -> sir1

simulate(sir1, seed = 1914679908L) -> sir1
plot(sir1)


# Complex SIR -------------------------------------------------------------

seas.sir.step <- "
  double rate[6];
  double dN[6];
  double Beta;
  double dW;
  Beta = exp(b1 + b2 * cos(M_2PI * Phi) + b3 * sin(M_2PI * Phi));
  rate[0] = births;                // birth
  rate[1] = Beta * (I + iota) / P; // infection
  rate[2] = mu;                    // death from S
  rate[3] = gamma;                 // recovery
  rate[4] = mu;                    // death from I
  rate[5] = mu;                    // death from R
  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);
  dW = rnorm(dt, sigma * sqrt(dt));
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  P = S + I + R;
  Phi += dW;
  H += dN[1];
  noise += (dW - dt) / sigma;
"

pomp(sir1, rprocess = euler.sim(
  step.fun = Csnippet(seas.sir.step), delta.t = 1/52/20),
  dmeasure = Csnippet(dmeas), rmeasure = Csnippet(rmeas),
  tcovar = "time", zeronames = c("H", "noise"),
  statenames = c("S", "I", "R", "H", "P", "Phi", "noise"),
  paramnames = c("gamma", "mu", "popsize", "rho","theta","sigma",
                 "S.0", "I.0", "R.0", "b1", "b2", "b3", "iota"),
  initializer = function(params, t0, ...) {
    fracs <- params[c("S.0", "I.0", "R.0")]
    setNames(c(round(params["popsize"]*c(fracs/sum(fracs),1)),0,0,0),
             c("S","I","R","P","H","Phi","noise"))
    },
  params = c(popsize = 500000, iota = 5, b1 = 6, b2 = 0.2,
             b3 = -0.1, gamma = 26, mu = 1/50, rho = 0.1, theta = 100,
             sigma = 0.3, S.0 = 0.055, I.0 = 0.002, R.0 = 0.94)) -> sir2

simulate(sir2, seed = 619552910L) -> sir2
plot(sir2)
