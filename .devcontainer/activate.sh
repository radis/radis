#!/bin/bash

conda_env=$(conda env list | grep radis-env)

echo $conda_env

if [[ $conda_env == *"radis-env"* ]]; then
    pip install -e .[dev] -v
    exit 0
fi

conda activate radis-env
pip install -e .[dev] -v

