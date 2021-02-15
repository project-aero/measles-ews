#!/usr/bin/env bash

dvc run \
    -d code/simulate-emergence-grid.R \
    -d code/make-pomp-simulator-function.R \
    -d results/initial-mif-lls-Agadez.csv \
    -d results/initial-mif-lls-Maradi.csv \
    -d results/initial-mif-lls-Niamey.csv \
    -d results/initial-mif-lls-Zinder.csv \
    -d code/measles-pomp-object-Agadez.RDS \
    -d code/measles-pomp-object-Maradi.RDS \
    -d code/measles-pomp-object-Niamey.RDS \
    -d code/measles-pomp-object-Zinder.RDS \
    -o simulations/emergence-simulations-grid-Agadez-1e-04.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.1.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.2.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.3.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.4.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.5.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.6.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.7.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.8.RDS \
    -o simulations/emergence-simulations-grid-Agadez-0.9.RDS \
    -o simulations/emergence-simulations-grid-Agadez-1.RDS \
    -o simulations/emergence-simulations-grid-Maradi-1e-04.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.1.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.2.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.3.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.4.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.5.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.6.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.7.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.8.RDS \
    -o simulations/emergence-simulations-grid-Maradi-0.9.RDS \
    -o simulations/emergence-simulations-grid-Maradi-1.RDS \
    -o simulations/emergence-simulations-grid-Zinder-1e-04.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.1.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.2.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.3.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.4.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.5.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.6.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.7.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.8.RDS \
    -o simulations/emergence-simulations-grid-Zinder-0.9.RDS \
    -o simulations/emergence-simulations-grid-Zinder-1.RDS \
    -o simulations/emergence-simulations-grid-Niamey-1e-04.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.1.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.2.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.3.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.4.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.5.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.6.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.7.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.8.RDS \
    -o simulations/emergence-simulations-grid-Niamey-0.9.RDS \
    -o simulations/emergence-simulations-grid-Niamey-1.RDS \
    --force \
    -n simulate-for-emergence-ews \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript simulate-emergence-grid.R
    
dvc run \
    -d code/analyze-emergence-grid-sims.R \
    -d simulations/emergence-simulations-grid-Agadez-1e-04.RDS \
    -d simulations/emergence-simulations-grid-Agadez-0.1.RDS \
    -d simulations/emergence-simulations-grid-Agadez-0.2.RDS \
    -d simulations/emergence-simulations-grid-Agadez-0.3.RDS \
    -d simulations/emergence-simulations-grid-Agadez-0.4.RDS \
    -d simulations/emergence-simulations-grid-Agadez-0.5.RDS \
    -d simulations/emergence-simulations-grid-Maradi-1e-04.RDS \
    -d simulations/emergence-simulations-grid-Maradi-0.1.RDS \
    -d simulations/emergence-simulations-grid-Maradi-0.2.RDS \
    -d simulations/emergence-simulations-grid-Maradi-0.3.RDS \
    -d simulations/emergence-simulations-grid-Maradi-0.4.RDS \
    -d simulations/emergence-simulations-grid-Maradi-0.5.RDS \
    -d simulations/emergence-simulations-grid-Zinder-1e-04.RDS \
    -d simulations/emergence-simulations-grid-Zinder-0.1.RDS \
    -d simulations/emergence-simulations-grid-Zinder-0.2.RDS \
    -d simulations/emergence-simulations-grid-Zinder-0.3.RDS \
    -d simulations/emergence-simulations-grid-Zinder-0.4.RDS \
    -d simulations/emergence-simulations-grid-Zinder-0.5.RDS \
    -d simulations/emergence-simulations-grid-Niamey-1e-04.RDS \
    -d simulations/emergence-simulations-grid-Niamey-0.1.RDS \
    -d simulations/emergence-simulations-grid-Niamey-0.2.RDS \
    -d simulations/emergence-simulations-grid-Niamey-0.3.RDS \
    -d simulations/emergence-simulations-grid-Niamey-0.4.RDS \
    -d simulations/emergence-simulations-grid-Niamey-0.5.RDS \
    -o results/emergence-bandwidths.csv \
    -o results/ews-emergence.csv \
    -o results/emergence-grid-aucs.csv \
    --force \
    -n analyze-emergence-ews \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript analyze-emergence-grid-sims.R