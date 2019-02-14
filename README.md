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
3. Simulate replicate emergence and elimination time series using fitted model parameters.
4. Calculate early warning signals (EWS) in null and test intervals of the replicate simulation time series.
5. Compute Area Under the Curve values for each EWS using the distribution of EWS calculated from the replicate simulations.

---

### Contact
Andrew Tredennick [atredenn@gmail.com](mailto:atredenn@gmail.com]
