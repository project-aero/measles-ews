# make-pomp-filtering-function.R
#  R function to generate a pomp model for running a particle filter
#  to estimate states and make one-week-ahead predictions. This version
#  of the pomp model includes a random walk on the transmission rate
#  so that the "best" filtered value is chosen at each time step. Filtering
#  is done at the MLE parameters.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)

make_pomp_filter <- function(obs_data, covar_data, mles, rw_value = 0.001)
{
  library(tidyverse)
  library(pomp)
  
  # Define stochastic process (SDEs) 
  
  measles_process <- Csnippet(
    "
    // Define the variables
    int nrate = 3;          // number of rates
    double rate[nrate];	  	// transition rates
    double trans[nrate];   	// transition numbers
    double lambda;          // force of infection
    double beta;            // transmission rate
    double eta = 365/8;     // infectious rate (8 days latent)
    double gamma = 365/5;   // recovery rate (5 days infectious)
    double dW;              // white noise
    double seas;            // seasonality term
    double dN0S, dN0I, dNSE, dNEI, dNIR;  // transitions
    
    // Beta random walk
    beta_t *= rgammawn(rw_intensity, dt)/dt;

    // Calculate force of infection
    seas = (1 + exp(dot_product(K, &xi1, &b1)));
    beta = beta_t*seas;
    lambda = beta*I/N;
    
    // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
    dW = rgammawn(beta_sd, dt);
    
    // Compute the transition rates
    rate[0] = lambda*dW/dt; // force of infection
    rate[1] = eta;          // infectious rate from latent
    rate[2] = gamma;	      // recovery from infectious
    
    // Compute the state transitions
    reulermultinom(1, S, &rate[0], dt, &trans[0]);
    reulermultinom(1, E, &rate[1], dt, &trans[1]);
    reulermultinom(1, I, &rate[2], dt, &trans[2]);
    
    // Transitions
    dN0S = rpois(0.3 * mu * N * dt);
    dN0I = rpois(iota * dt);
    dNSE = trans[0];
    dNEI = trans[1];
    dNIR = trans[2];
    
    // Balance the equations
    S += dN0S - dNSE;
    E +=        dNSE - dNEI;
    I += dN0I        + dNEI - dNIR;
    
    cases += dNIR;  // cases are cumulative reports at end of infectious period (I->R)
    RE = (beta_t / gamma) * (S / N);
    RE_seas = (beta / gamma) * (S / N);
    cases_state = rnbinom_mu(1/tau, rho*cases);
    "
  )
  
  
  # Define likelihood function 
  
  measles_dmeasure <- Csnippet(
    "
    double mean;
    double f;
    mean = cases*rho;
    if (ISNA(reports)) {  // for catching missing observations
    lik = (give_log) ? 0 : 1;
    } else {
    f = dnbinom_mu(reports, 1/tau, mean, give_log);  // neg binomial likelihood
    }
    
    lik = (give_log) ? log(f) : f;
    "
  )
  
  
  # Define process simulator for observations 
  
  measles_rmeasure <- Csnippet(
    "
    reports = rnbinom_mu(1/tau, rho*cases);  // neg binomial measurement process
    
    if (reports > 0.0) {
    reports = nearbyint(reports);
    } else {
    reports = 0.0;
    }
    "
  )
  
  
  # Define parameter transformation scales 
  
  from_estimation <- Csnippet(
    "
    Tbeta_mu = exp(beta_mu);
    Tiota = exp(iota);
    Trho = expit(rho);
    Tbeta_sd = exp(beta_sd);
    TS_0 = expit(S_0);
    TE_0 = expit(E_0);
    TI_0 = expit(I_0);
    Ttau = exp(tau);
    "
  )
  
  to_estimation <- Csnippet(
    "
    Tbeta_mu = log(beta_mu);
    Tiota = log(iota);
    Trho = logit(rho);
    Tbeta_sd = log(beta_sd);
    TS_0 = logit(S_0);
    TE_0 = logit(E_0);
    TI_0 = logit(I_0);
    Ttau = log(tau);
    "
  )
  
  initial_values <- Csnippet(
    "
    S = nearbyint(N*S_0);
    E = nearbyint(N*E_0);
    I = nearbyint(N*I_0);
    cases = nearbyint(N*I_0);
    cases_state = nearbyint(N*I_0);
    RE = 0;
    RE_seas = 0;
    beta_t = beta_mu;
    "
  )
  
  params <- unlist(mles)
  
  global_str <- paste0(
    "int K = 6;", " double rw_intensity = ", rw_value, ";"
  )
  
  filtering_pomp <- pomp(
    data = obs_data,
    times = "time",
    covar = covar_data,
    tcovar = "time",
    t0 = min(obs_data$time),
    rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
    rmeasure = measles_rmeasure,
    dmeasure = measles_dmeasure,
    initializer = initial_values,
    statenames = c("S", "E", "I", "cases", "cases_state", 
                   "RE", "RE_seas", "beta_t"),
    toEstimationScale = to_estimation,
    fromEstimationScale = from_estimation,
    paramnames = names(params),
    params = params,
    globals = global_str,
    zeronames = c("cases", "cases_state")
  )
  
  return(filtering_pomp)
}
