#!/bin/bash

# Install this repository to a new Python environment

# Install pip first so we don't install the package to base
mamba install pip || conda install pip

# Install this package's routines
pip install .

# Create input/output directories
mkdir input
mkdir output

# Install CLIMADA