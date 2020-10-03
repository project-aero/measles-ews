#!/usr/bin/env bash

dvc run \
    -d data/raw-data \
    -d code/fetch-clean-data.R \
    -o data/clean-data \
    --force \
    -n clean-data \
    podman run -w /root/code --rm \
    --mount type=bind,src='$PWD',dst=/root \
    docker.io/eamon/2019measles:v20200928 \
    Rscript fetch-clean-data.R 
