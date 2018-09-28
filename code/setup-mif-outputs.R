# setup-mif-outputs.R
#  This script saves empty CSVs for storing results from MIF when run
#  on the HPC.
#
# Author:
#  Andrew Tredennick


# Set up output files -----------------------------------------------------

mles <- data.frame(
  do_grid = NA,
  loglik = NA,
  loglik_se = NA,
  beta_mu = NA,
  beta_sd = NA,
  b1 = NA,
  b2 = NA,
  b3 = NA,
  b4 = NA,
  b5 = NA,
  b6 = NA,
  iota = NA,
  rho = NA,
  S_0 = NA,
  E_0 = NA,
  I_0 = NA,
  tau = NA
)

ll_file <- "initial-mif-lls.csv"
write.table(mles, ll_file, sep = ",", col.names = T, append = T, row.names = FALSE)

mf_traces <- data.frame(
  do_grid = NA,
  iteration = NA,
  loglik = NA,
  nfail = NA,
  beta_mu = NA,
  beta_sd = NA,
  b1 = NA,
  b2 = NA,
  b3 = NA,
  b4 = NA,
  b5 = NA,
  b6 = NA,
  iota = NA,
  rho = NA,
  S_0 = NA,
  E_0 = NA,
  I_0 = NA,
  tau = NA
)

trace_file <- "initial-mif-traces.csv"
write.table(mf_traces, trace_file, sep = ",", col.names = T, append = T, row.names = FALSE)