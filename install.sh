#!/bin/bash

# Install this repository to a new Python environment

# Install pip first so we don't install the package to base
mamba install pip || conda install pip

# Install CLIMADA
curl -o env_climada.yml https://raw.githubusercontent.com/CLIMADA-project/climada_python/main/requirements/env_climada.yml
conda env update -n my_env --file env_climada.yaml
pip install climada

# Install this package's routines
pip install .

# Create input/output directories
mkdir input
mkdir output
mkdir plots