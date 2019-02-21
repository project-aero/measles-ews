# measles-ews
Data and code for project aimed at detecting critical slowing down in measles dynamics from four cities in Niger.

**NOTE**: The code in this repository relies on `pomp` version 1.18. 
You may get errors if using other versions. 
To download this specific version, use `devtools`:

```
install_packages("devtools")
devtools::install_version("pomp", version = "1.18", repos = "http://cran.us.r-project.org")
```

## Overview
The main goal of this project is to detect critical slowing in real disease dynamics.
To do so, we use a model-based approach, which allows us to simulate approaches to transcritical bifurcations.
The model, however, is empirically-based because we use SEIR models fit to time series data from four cities in Niger.
Thus, this work represents a step in the direction of an empirical test of CSD, though with still a lot of help from the model.
We use early warning signals to detect critical slowing down.

There are five main steps in the analysis.

1. Fit mechanistic SEIR models to weekly case count data from four cities in Niger (requires High Performance Computing).
2. Estimate uncertainty around maximum likelihood parameters via parametric bootstrap (requires High Performance Computing).
3. Simulate replicate emergence and elimination time series using fitted model parameters (does not require High Performance Computing, but is sped up by using multiple cores on the Drake Lab high memory machine).
4. Calculate early warning signals (EWS) in null and test intervals of the replicate simulation time series (does not require High Performance Computing, but is sped up by using multiple cores on the Drake Lab high memory machine).
5. Compute Area Under the Curve values for each EWS using the distribution of EWS calculated from the replicate simulations (no High Performance Computing!).

## Notebooks
Many of my ideas and work in progress beyond that presented in the manuscript file and the description in this ReadMe can be found in a set of R Markdown Notebooks in the `measles-ews/docs/notebooks/` subdirectory.
Likewise, progress during the model fitting stage -- different model fitting strategies, approaches for estimating parameter uncertainty -- can be found in the physical TREDENNICK Lab Notebook (#00074).

## Directory and file information
- `aux/`: Contains R scripts that are no longer in use, either because they were dead-ends or because they were revamped for production.
- `code/`: Contains R scripts necessary to reproduce the all analyses, from model fitting to simulating to calculation of early warning signals. I describe the key scripts below, other scripts are mainly for plotting/checking intermediate results.
  + `analyze-elimination-grid-sims.R`: Calculates early warning signals in two windows from the replicate elimination simulations and the area under the curve statistics.
  + `analyze-emergence-grid-sims.R`: Calculates early warning signals in two windows from the replicate emergence simulations and the area under the curve statistics. 
  + `analyze-mvw-elimination-grid-sims.R`: Calculates early warning signals in moving windows from the replicate elimination simulations and the area under the curve statistics. Only presented in supplementary material.
  + `analyze-mvw-emergence-grid-sims.R`: Calculates early warning signals in moving windows from the replicate emergence simulations and the area under the curve statistics. Only presented in supplementary material.
  + `boot-mif-job*.sh`: Bash scripts (* = 1,2,3,4,5) for running maximization by iterated filtering (MIF) on the Olympus High Performance Computing cluster. The scripts run instances of the `bootstrap-fit-mif.R` script in parallel, in 1000 core batches.
  + `bootstrap-fit-mif.R`: Fits an *SEIR* model to weekly case count data using MIF as implemented with the `pomp::mif()` function. In this case, the "data" are stochastic realizations from the fitted model -- thus, model fitting is for a parametric bootstrap. This script **designed to be run on a High Performance Computing cluster only**. See `boot-mif-job*.sh` description above.
  + `define-continuous-measles-pomp.R`: This script generates `pomp` models to be used for model fitting. One `pomp` model is generated for each city and saved as `measles-pomp-object-*.RDS`, where * is the name of the focal city. These files are stored in the code directory for easy access.
  + `estimate-transmission-state.R`: This script peforms a plain vanilla particle filter at the fitted MLE parameters to generate one-step-ahead predictions for calculating a generalized *R*<sup>2</sup> for model fit. It was also initially designed to estimate *R*<sub>E</sub> over time when transmission rate was allowed to take a random walk. However, further analyses showed that the highest likelihood models do not include a random walk tranmission rate. In other words, we found no evidence of a directional temporal trend in transmission rate.
  + `fetch-clean-data.R`: This scripts reads in the raw case count and demographic data to make clean data structures for use in the `pomp` models and elsewhere. Clean data reside in the `data/clean-data/` subdirectory.
  + `find-beta-random-walk-intensity.R`: Runs the `pomp` MIF2 algorithm starting at the MLEs but let's tranmission rate take a random walk. The idea is to test for a directional trend in tranmission. The models suggest there is no directional trend.
  + `global-search-mif.R`: This is the key script for model fitting. It implements the `pomp` MIF2 algorithm across a grid of starting values for model parameters. This script **designed to be run on a High Performance Computing cluster only**. Depends on `fetch-clean-data.R` and `define-continuous-measles-pomp.R`.
  + `initial-mif-job*.sh`: Bash scripts (* = 1,2,3,4,5) for running maximization by iterated filtering (MIF) on the Olympus High Performance Computing cluster. The scripts run instances of the `global-search-mif.R` script in parallel, in 1000 core batches.
  + `make-pomp-filtering-function.R`: Defines a function for generating a `pomp` model for particle filtering.
  + `make-pomp-simulating-function.R`: Defines a function for generating a `pomp` model for making long run simulations. The function is called to make simulators of emergence and elimination.
  + `setup-mif-outputs.R`: A utility script to generates empty CSV files for storage of model parameters from MIF on the Olympus HPC.
  + `simulate-bootstrap-series.R`: Simulates replicate stochastic realizations of 11 years of weekly case counts from the MLE parameter sets for each city. These series are used for fitting parameters again to perform a parametric bootstrap and estimate uncertainty around parameter values.
  + `simulate-elimination-grid.R`: Simulates replicate elimination scenarios at different speeds of vaccination uptake.
  + `simulate-emergence-grid.R`: Simulates replicate emergence scenarios at different levels of susceptible depletion.
- `data/`: Contains raw and cleaned versions of case count and demographic data. See ReadMes therein and the supplementary material on demographic data. Note that the raw data are provided to Project AERO by Matthew Ferrari. All subsequent use should be approved by M. Ferrari, J. Drake, and P. Rohani. There is also a subdirectory of spatial data for plotting the location of the four focal cities in this analysis.
- `docs/`: Contains the manuscript and supplementary information documents. These are `.Rmd` files that are knitted to PDF. This directory also contains subdirectories for project notebooks (`notebooks/`) and the original project protocol (`protocol/`). Other subdirectories are made on the fly when knitting the Rmd to PDF.
- `figures/`: This directory is now deprecated. It contains figures for the manuscript and SI that I am keeping here for posterity, but all manuscript and SI figures are made now made on the fly in the Rmd docs.
- `results/`: Contains CSV and RDS files of analysis results. These are either intermediary and called by scripts for further analysis or are used in manuscript/SI docs to make figures/tables.
- `simulations/`: Contains RDS files that store simulation results.

---

### Contact
Andrew Tredennick ([atredenn@gmail.com](mailto:atredenn@gmail.com))
