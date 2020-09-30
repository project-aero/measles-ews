# make-pomp-simulator-function.R
#  R function to generate a pomp model for simulating measles dynamics based
#  on fitted model parameters, but with no further relation to data.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)

make_pomp_simulator <- function(do_city, mles, years_to_sim = 30, 
                                initial_population_size, susc_discount = 1,
                                exposed_discount = 1, infected_discount = 1,
                                vacc_coverage_ts = NULL)
{
  library(tidyverse)
  library(pomp)
  
  # Define stochastic process (SDEs) 
  
  measles_process <- Csnippet(
    "
    // Define the variables
    int nrate = 7;            // number of rates
    double rate[nrate];	  	  // transition rates
    double trans[nrate];   	  // transition numbers
    double lambda;            // force of infection
    double beta;              // transmission rate
    double eta = 365/8;       // infectious rate (8 days latent)
    double gamma = 365/5;     // recovery rate (5 days infectious)
    double psi;               // rate of infection importation 
    double dW;                // white noise
    double seas;              // seasonality term
    double dNS0, dN0S, dNSE;  // S transitions
    double dNE0, dNEI;        // E transitions
    double dN0I, dNI0, dNIR;  // I transitions
    double dN0R, dNR0;        // R transitions

    // Calculate force of infection
    seas = (1 + exp(dot_product(K, &xi1, &b1)));
    beta = beta_mu*seas;
    lambda = beta*I/N;
    
    // Gamma noise, mean=dt, variance=(beta_sd^2 dt)
    dW = rgammawn(beta_sd, dt);
    
    // Compute the transition rates
    rate[0] = nu;           // susceptible deaths
    rate[1] = lambda*dW/dt; // force of infection
    rate[2] = nu;           // exposed deaths
    rate[3] = eta;          // rate of exit from latent compartment
    rate[4] = nu;           // infected deaths
    rate[5] = gamma;	      // recovery from infectious
    rate[6] = nu;           // recovered deaths
    
    // Compute the state transitions
    reulermultinom(2, S, &rate[0], dt, &trans[0]);
    reulermultinom(2, E, &rate[2], dt, &trans[2]);
    reulermultinom(2, I, &rate[4], dt, &trans[4]);
    reulermultinom(1, R, &rate[6], dt, &trans[6]);
    
    // Transitions
    dNS0 = trans[0];                                // susceptible deaths
    dN0S = rpois(vacc_discount * mu * N * dt);      // births, unvaccinated
    dNSE = trans[1];                                // S -> E
    dNE0 = trans[2];                                // exposed deaths
    dNEI = trans[3];                                // E -> I
    dN0I = rpois(psi * dt);                         // imported infections
    dNI0 = trans[4];                                // infected deaths
    dNIR = trans[5];                                // I -> R
    dN0R = rpois((1-vacc_discount) * mu * N * dt);  // births, vaccinated
    dNR0 = trans[6];                                // recovered deaths
    
    // Balance the equations
    S += dN0S - dNSE               - dNS0;
    E +=        dNSE - dNEI        - dNE0;
    I += dN0I        + dNEI - dNIR - dNI0;
    R += dN0R               + dNIR - dNR0;
    
    cases += dNIR;  // cumulative reports at end of infectious period (I->R)
    RE_seas = ((eta*beta) / ((eta+nu)*(gamma+nu))) * (S / N);
    // RE_seas = (beta / gamma) * (S / N);  // simple formula for fitting model
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
    S = nearbyint(N*S_0*sfact);
    E = nearbyint(N*E_0*efact);
    I = nearbyint(N*I_0*sfact);
    R = nearbyint(N - (N*S_0*sfact + N*E_0*efact + N*I_0*ifact));
    cases = rho*N*I_0;
    RE_seas = 0;
    "
  )
  
  all_times <- seq(0, years_to_sim, by = 1/365)
  report_times <- all_times[seq(1, length(all_times), 7)]
  
  # Make data tables
  the_data <- tibble(
    time = report_times,
    reports = NA
  ) %>%
    filter(time > 0.5)  # start in trough of seasonal transmission
  
  # Generate basis functions for seasonality
  covar_table <- periodic.bspline.basis(
    all_times,
    nbasis = 6,
    degree = 3,
    period = 1,
    names = "xi%d"
  ) %>%
    as_tibble() %>%
    mutate(
      time = all_times
    )
  
  if(is.null(vacc_coverage_ts) == FALSE){
    covar_table <- covar_table %>%
      mutate(
        vacc_discount = 1 - vacc_coverage_ts
      )
  }
  
  sfact <- susc_discount
  efact <- exposed_discount
  ifact <- infected_discount
  N1 <- initial_population_size
  global_str <- paste0(
    "int K = 6; double mu = 0.05; double nu = 0.05;", 
    " double sfact = ", sfact, ";",
    " double efact = ", efact, ";",
    " double ifact = ", ifact, ";",
    " double N = ", N1, ";"
  )
  
  if(is.null(vacc_coverage_ts) == TRUE){
    vacc_discount <-  0.3
    global_str <- paste0(
      "int K = 6; double mu = 0.05; double nu = 0.05;", 
      " double sfact = ", sfact, ";",
      " double efact = ", efact, ";",
      " double ifact = ", ifact, ";",
      " double N = ", N1, ";",
      " double vacc_discount = ", vacc_discount, ";"
    )
  }
  
  simulator_pomp <- pomp(
    data = the_data,
    times = "time",
    covar = covar_table,
    tcovar = "time",
    t0 = min(the_data$time),
    rprocess = euler.sim(step.fun = measles_process, delta.t = 1/365),
    rmeasure = measles_rmeasure,
    dmeasure = measles_dmeasure,
    initializer = initial_values,
    statenames = c("S", "E", "I", "R", "cases", "RE_seas"),
    toEstimationScale = to_estimation,
    fromEstimationScale = from_estimation,
    paramnames = names(mles),
    params = unlist(mles),
    globals = global_str,
    zeronames = c("cases")
  )
  
  return(simulator_pomp)
}
