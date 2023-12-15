# Coastal risk assessment using ROMS

This directory contains scripts and executables for running the Regional Ocean Model System (ROMS) to assess coastal risk due to tropical cyclone storm surges.

## Setup

1. Use env_surge.yml and conda/mamba to create a new Python environment with the necessary dependencies:
``` conda create -n roms_env -f env_surge.yml ```
2. Install the modules for `crisk_surge`:
``` pip install -e . ```
The `-e` flag means that the modules are editable. The functions are found in `src/crisk_surge`.
3. To generate ROMs grids, the PyROMS library is used, which requires a separate installation. It is best to clone the PyROMS repository and install it using `pip`. See here for instructions. Sometimes PyROMS will have trouble finding <library> on the system. In this case, go into `pyroms/pyroms/src/...` and change line X to reflect its path on your system. This is probably somewhere like `usr/lib`.

## How to use

## How accurate is this approach?
