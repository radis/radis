#!/bin/bash

conda_env=$(conda env list | grep radis-env)

echo $conda_env

if [[ $conda_env == *"radis-env"* ]]; then
    exit 0
fi

conda env create --file environment.yml --solver=libmamba
conda init 
