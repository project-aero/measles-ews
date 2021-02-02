#!/usr/bin/env bash

dvc run \
    -d code/simulate-elimination-grid.R \
    -d code/make-pomp-simulator-function.R \
    -d results/initial-mif-lls-Agadez.csv \
    -d results/initial-mif-lls-Maradi.csv \
    -d results/initial-mif-lls-Niamey.csv \
    -d results/initial-mif-lls-Zinder.csv \
    -d code/measles-pomp-object-Agadez.RDS \
    -d code/measles-pomp-object-Maradi.RDS \
    -d code/measles-pomp-object-Niamey.RDS \
    -d code/measles-pomp-object-Zinder.RDS \
    -o simulations/elimination-simulations-grid-Maradi-1.5e-05.RDS \
    -o simulations/elimination-simulations-grid-Maradi-1e-04.RDS \
    -o simulations/elimination-simulations-grid-Maradi-3.2e-05.RDS \
    -o simulations/elimination-simulations-grid-Maradi-4.9e-05.RDS \
    -o simulations/elimination-simulations-grid-Maradi-6.6e-05.RDS \
    -o simulations/elimination-simulations-grid-Maradi-8.3e-05.RDS \
    -o simulations/elimination-simulations-grid-Agadez-1.5e-05.RDS \
    -o simulations/elimination-simulations-grid-Agadez-1e-04.RDS \
    -o simulations/elimination-simulations-grid-Agadez-3.2e-05.RDS \
    -o simulations/elimination-simulations-grid-Agadez-4.9e-05.RDS \
    -o simulations/elimination-simulations-grid-Agadez-6.6e-05.RDS \
    -o simulations/elimination-simulations-grid-Agadez-8.3e-05.RDS \
    -o simulations/elimination-simulations-grid-Niamey-1.5e-05.RDS \
    -o simulations/elimination-simulations-grid-Niamey-1e-04.RDS \
    -o simulations/elimination-simulations-grid-Niamey-3.2e-05.RDS \
    -o simulations/elimination-simulations-grid-Niamey-4.9e-05.RDS \
    -o simulations/elimination-simulations-grid-Niamey-6.6e-05.RDS \
    -o simulations/elimination-simulations-grid-Niamey-8.3e-05.RDS \
    -o simulations/elimination-simulations-grid-Zinder-1.5e-05.RDS \
    -o simulations/elimination-simulations-grid-Zinder-1e-04.RDS \
    -o simulations/elimination-simulations-grid-Zinder-3.2e-05.RDS \
    -o simulations/elimination-simulations-grid-Zinder-4.9e-05.RDS \
    -o simulations/elimination-simulations-grid-Zinder-6.6e-05.RDS \
    -o simulations/elimination-simulations-grid-Zinder-8.3e-05.RDS \
    --force \
    -n simulate-for-elimination-ews \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript simulate-elimination-grid.R

dvc run \
    -d code/analyze-elimination-grid-sims.R \
    -d simulations/elimination-simulations-grid-Maradi-1.5e-05.RDS \
    -d simulations/elimination-simulations-grid-Maradi-1e-04.RDS \
    -d simulations/elimination-simulations-grid-Maradi-3.2e-05.RDS \
    -d simulations/elimination-simulations-grid-Maradi-4.9e-05.RDS \
    -d simulations/elimination-simulations-grid-Maradi-6.6e-05.RDS \
    -d simulations/elimination-simulations-grid-Maradi-8.3e-05.RDS \
    -d simulations/elimination-simulations-grid-Agadez-1.5e-05.RDS \
    -d simulations/elimination-simulations-grid-Agadez-1e-04.RDS \
    -d simulations/elimination-simulations-grid-Agadez-3.2e-05.RDS \
    -d simulations/elimination-simulations-grid-Agadez-4.9e-05.RDS \
    -d simulations/elimination-simulations-grid-Agadez-6.6e-05.RDS \
    -d simulations/elimination-simulations-grid-Agadez-8.3e-05.RDS \
    -d simulations/elimination-simulations-grid-Niamey-1.5e-05.RDS \
    -d simulations/elimination-simulations-grid-Niamey-1e-04.RDS \
    -d simulations/elimination-simulations-grid-Niamey-3.2e-05.RDS \
    -d simulations/elimination-simulations-grid-Niamey-4.9e-05.RDS \
    -d simulations/elimination-simulations-grid-Niamey-6.6e-05.RDS \
    -d simulations/elimination-simulations-grid-Niamey-8.3e-05.RDS \
    -d simulations/elimination-simulations-grid-Zinder-1.5e-05.RDS \
    -d simulations/elimination-simulations-grid-Zinder-1e-04.RDS \
    -d simulations/elimination-simulations-grid-Zinder-3.2e-05.RDS \
    -d simulations/elimination-simulations-grid-Zinder-4.9e-05.RDS \
    -d simulations/elimination-simulations-grid-Zinder-6.6e-05.RDS \
    -d simulations/elimination-simulations-grid-Zinder-8.3e-05.RDS \
    -o results/elimination-bandwidths.csv \
    -o results/ews-elimination.csv \
    -o results/elimination-grid-aucs.csv \
    -o simulations/single-elimination-example.RDS \
    --force \
    -n analyze-elimination-ews \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript analyze-elimination-grid-sims.R