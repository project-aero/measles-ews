lna_stan_code <- "
  functions {
        // System of ODEs to be integrated to compute the LNA moments - state is an array with the current state
        real[] lna_fcn(real t, real[] state, real[] lna_params, real[] rdummy, int[] idummy) {
        
                // compute hazards and jacobian
                vector[2] haz;
                vector[2] phi;
                matrix[2,2] jac;
                real newstate[6];
                
                // extract the current drift and diffusion
                vector[2] Z      = to_vector(state[1:2]);
                vector[2] neg_Z  = -Z;
                vector[2] neg_2Z = -2*Z;
                vector[2] exp_Z  = exp(Z);
                matrix[2,2] V0   = to_matrix(state[3:6], 2, 2);
                
                // compute compartment volumes
                real S = lna_params[3] - exp_Z[1]; // S0 - N_{SI}
                real I = lna_params[4] + exp_Z[1] - exp_Z[2];
                
                // compute some quantities that will be reused
                vector[2] ito_diff = exp(neg_Z) - 0.5 * exp(neg_2Z); // exp(log(exp(-Z) - 0.5*exp(-2*Z)))
                vector[2] d_ito_diff = exp(neg_2Z) - exp(neg_Z); // exp(log(exp(-2*Z) - exp(-Z)))
                
                // numerically stable versions of the above, but again broken. not sure why.
                // real I = lna_params[4] + exp(log_diff_exp(Z[1], Z[2])); // I0 + N_{SI} - N_{IR}
                
                // vector[2] ito_diff;
                // vector[2] d_ito_diff;
                // ito_diff[1] = exp(log_diff_exp(neg_Z[1], log(0.5) + neg_2Z[1]));
                // ito_diff[2] = exp(log_diff_exp(neg_Z[2], log(0.5) + neg_2Z[2]));
                // 
                // d_ito_diff[1] = exp(log_diff_exp(neg_2Z[1], neg_Z[1]));
                // d_ito_diff[2] = exp(log_diff_exp(neg_2Z[1], neg_Z[1]));
                
                // compute the hazards
                phi[1] = lna_params[1]*S*I;
                phi[2] = lna_params[2]*I;
                haz = ito_diff .* phi;
                
                // compute the jacobian
                jac[1,1] = lna_params[1] * (d_ito_diff[1]*S*I - ito_diff[1]*exp_Z[1]*I + ito_diff[1]*S*exp_Z[1]);
                jac[1,2] = -lna_params[1] * ito_diff[1] * S * exp_Z[2];
                jac[2,1] = lna_params[2] * ito_diff[2] * exp_Z[1];
                jac[2,2] = lna_params[2] * (d_ito_diff[2]*I - ito_diff[2]*exp_Z[2]);
                
                // concatenate and return
                newstate[1:2] = to_array_1d(haz);
                newstate[3:6] = to_array_1d(V0*jac' + jac*V0);
                newstate[3] = newstate[3] + exp(neg_2Z[1])*haz[1];
                newstate[6] = newstate[6] + exp(neg_2Z[2])*haz[2];
                
                return newstate;
        }
        
        // integrate over one time interval
        real[] lna_step(real t_l, real[] t_r, real[] state, real[] lna_params, real[] rdummy, int[] idummy) {
        
                // changing the integration method to rk45 results in an error.
                return integrate_ode_bdf(lna_fcn, state, t_l, t_r, lna_params, rdummy, idummy, 1e-6, 1e-6, 100000)[1];
        }
        
        // function to map the standard normal draws to the log-LNA increments and return incidence (increments in N_SI)
        vector get_lna(vector N_raw, real[] theta, vector X0, matrix stoich, vector net_effect, real[] times, real[] rdummy, int[] idummy) {
                // number of times at which incidence is recorded
                int n_times = size(times) - 1;
                
                // containers current state
                vector[2] log_LNA;                      // log-LNA increment: (log(N_SI), log(N_IR))
                vector[2] nat_LNA;                      // LNA increment on the natural scale
                vector[2] c_incid = rep_vector(0.0, 2); // cumulative incidence
                vector[3] SIR_init = X0 - net_effect;   // initial value for SIR compartment volumes
                vector[3] SIR_cur = SIR_init;           // SIR compartment volumes
                vector[n_times] incidence;              // increments for incidence (natural scale)
                real lna_params[5];                     // (beta, mu, S_t, I_t, R_t)
                
                // containers for LNA moments
                real zero_6[6] = rep_array(0.0, 6); // for zero-ing out the state vector
                real state[6]  = zero_6;            // LNA moments - (mu, Sigma) - systems of ODEs
                
                vector[2] mu;         // LNA mean
                matrix[2,2] Sig_chol; // Cholesky decomposition of the covariance
                
                // initialize lna parameters
                lna_params[1:2] = theta; // parameters
                
                // map the standard normal draws to the log-LNA increments
                for(k in 1:n_times) {
                      
                      // set the comaprtment volumes
                      lna_params[3:5] = to_array_1d(SIR_cur);
                
                      // LNA transition density
                      state = lna_step(times[k],
                      times[(k+1):(k+1)],
                      state,
                      lna_params,
                      rdummy,
                      idummy);
                      
                      // extract moments
                      mu        = to_vector(state[1:2]);
                      Sig_chol  = cholesky_decompose(to_matrix(state[3:6], 2, 2)); // compute the cholesky decomp.
                      
                      // map N_raw values to the sampled LNA value
                      log_LNA = mu + Sig_chol * N_raw[(2*k - 1):(2*k)];
                      
                      // save new state
                      nat_LNA      = exp(log_LNA); // LNA increment on the natural scale
                      incidence[k] = nat_LNA[1];   // save incidence
                      
                      // update the cumulative incidence and set the new initial state
                      state   = zero_6;                       // reset the LNA moments
                      c_incid = c_incid + nat_LNA;            // update cumulative incidence
                      SIR_cur = SIR_init + stoich * c_incid;  // update the current compartment counts
                      lna_params[3:5] = to_array_1d(SIR_cur); // set the new initial values
                }
                
                return incidence;
        }
}
data {
        real popsize;
        int n_times;          // number of census interval endpoints
        real times[n_times];  // census interval endpoint times
        int cases[n_times-1]; // observed incidence
        matrix[3,2] stoich;   // stoichiometry matrix
        vector[3] net_effect; // net change in comparment volumes for one of each event
        vector<lower=0>[3] X0; // initial compartment volumes
}
transformed data{
        real log_popsize = log(popsize);
        real rdummy[0]; // no real data arguments to LNA ODEs
        int idummy[0];  // no integer data arguments to LNA ODEs
}
parameters {
        /*
        Raw model parameters:
        log_R0: log(beta*N/mu)
        log_mu: log recovery rate
        rho:    sampling probability
        phi:    negative binomial overdispersion
        N_raw:  multivariate standard normal draws that drive the stochasticity
        */
        real log_R0;
        real log_mu;
        real logit_rho;
        real log_phi;
        
        vector[2*(n_times-1)] N_raw; // initial state is fixed for now
}
transformed parameters {
        real<lower=0> theta[2]; // theta = (beta, mu)
        real<lower=0,upper=1> rho;
        real<lower=0> phi;
        vector<lower=0>[n_times-1] lna_incidence; // expected LNA incidence
        
        // transform parameters
        theta = {exp(log_R0 + log_mu - log_popsize), exp(log_mu)};
        rho   = inv_logit(logit_rho);
        phi   = exp(log_phi);
        lna_incidence = get_lna(N_raw, theta, X0, stoich, net_effect, times, rdummy, idummy); // map normal draws to LNA incidence path
}
model {
        // parameter sampling statements
        log_R0    ~ normal(0.45, 1);
        log_mu    ~ normal(1.2, 1);
        logit_rho ~ normal(0.1,1);
        log_phi   ~ normal(2.5, 0.5);
        
        // implies the LNA paths are multivariate normal with moments given by ODE solns
        N_raw ~ normal(0,1);
        
        // emission probabilities - parameterized by mean and overdispersion
        cases ~ neg_binomial_2(rho * lna_incidence, phi);
}
"

# Data are negative binomial samples of weekly incidence from an epidemic simulated under SIR
# dynamics with a basic reproductive number of R_0 = 1.5 and mean infectious period duration 
# of 1 week. Cases were detected with negative binomial mean of 0.2, and negative binomial 
# overdispersion coefficient of 5. The initial time is 0.
incidence_data <- data.frame(time=1:15, cases = c(3,3,19,17,39,34,23,105,63,227,227,69,92,51,39))

# list containing the arguments to the data block. In this case, we assume that the initial
# state of the population at week 0 is known. 
dat <- list(
  popsize = 10015,
  times   = as.numeric(0:15),                               # times at which to evaluate the path (t0 = 0)
  n_times = 16,                                             # number of times 
  stoich  = matrix(as.numeric(c(-1,1,0,0,-1,1)), ncol = 2), # stochiometry matrix = change in compartments from each event
  cases   = incidence_data[,2],                                 # dataset
  X0      = as.numeric(c(10000, 10, 5)),                    # initial state
  net_effect = as.numeric(c(-1,0,1))                        # net effect on compartment volumes from one of each transition
)

# initial values, centered at the truth for demonstration purposes
init_fcn <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(log_R0    = log(1.5) + rnorm(1,0,0.01),
       log_mu    = log(1) + rnorm(1, 0, 0.01),
       logit_rho = -log(1/0.2 - 1) + rnorm(1, 0, 0.1),
       log_phi   = log(5) + rnorm(1, 0, 0.1),
       N_raw     = array(rnorm(30), dim = c(30)))  
}

# run the MCMC
lna_results <- stan(model_code = lna_stan_code, 
                    data = dat, 
                    init = init_fcn,
                    iter = 20,
                    chains = 1, 
                    warmup = 10,
                    cores = 1,
                    control = list(adapt_delta = 0.95))

lna <- summary(lna_results)$summary
states <- lna[grep("lna",rownames(lna)),]
plot(states[,1], type = "l")
lines(states[,4], lty = 2)
lines(states[,8], lty = 2)
points(incidence_data)
