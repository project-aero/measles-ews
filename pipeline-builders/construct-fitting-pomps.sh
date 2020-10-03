#!/usr/bin/env bash

dvc run \
    -d code/define-continuous-measles-pomp.R \
    -d data/clean-data \
    -o code/measles-pomp-object-Agadez.RDS \
    -o code/measles-pomp-object-Maradi.RDS \
    -o code/measles-pomp-object-Zinder.RDS \
    -o code/measles-pomp-object-Niamey.RDS \
    --force \
    -n construct-fitting-pomps \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript define-continuous-measles-pomp.R 
