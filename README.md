# measles-ews
Data and code for project aimed at detecting critical slowing down in measles dynamics from four cities in Niger.

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
  +


---

### Contact
Andrew Tredennick [atredenn@gmail.com](mailto:atredenn@gmail.com]
