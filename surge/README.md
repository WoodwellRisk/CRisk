# Coastal risk assessment using ROMS

This directory contains scripts and functions for running the Regional Ocean Model System (ROMS) to assess coastal risk due to tropical cyclone storm surges. For a specified location, thousands of years of synthetically (statistically) generated tropical cyclones are used to force the ROMS ocean surface. Tracks from the STORM dataset (Bloemendaal et al., 2020) are expanded into 2D wind and pressure fields using the parametric model of Holland, (1980). Wind fields are then converted to wind stress at the ocean surface using a quadratic function of windspeed and a drag coefficient according to (Peng et al., 2020). The [ParaTC](https://pypi.org/project/paratc/) python package is used to generate forcing. A validation against observed storm surges at tide gauges shows a mean absolute error in maximum surge of ~14cm and a correlation of 84%. For more information on the methodology and accuracy assessment, see (methodologydoc).

## Setup

1. Use env_surge.yml and conda/mamba to create a new Python environment with the necessary dependencies:
``` conda create -n roms_env -f env_surge.yml ```
2. Install the modules for `crisk_surge` by running the install inside `CRISK/surge`:
``` pip install -e . ```
The `-e` flag means that the modules are editable. The functions are found in `src/crisk_surge` and be imported into Python using `import crisk_surge`, or by importing modules directory e.g. `from crisk_surge import validation postprocessing forcing` 
3. To generate ROMs grids, the PyROMS library is used, which requires a separate installation. You can follow their instructions to do this.

## Structure

A new 'project' directory is created for each area of study and model domain. Each project directory is stored and controlled from the `/projects/` directory in this repository. Each project directory eventually contains:

1. `roms_grd.nc` : The ROMS grid file. This is generated using the `make_grid.py` script in `projects/`. This script uses PyROMS (and by extension `GridGen`) to generate a rectangular ROMS grid file from the command line. 
2. `roms.in.template` : The ROMS control input file. This is a template file, meaning there are missing lines waiting to be filled by our python scripts. The up to date version of this files is in `surge/ROMS`.
3. `romsM`: The ROMS executable. You can generate this by following the ROMS compilation instructions and using `ROMS/stormsurge.h` as a header file. Once compiled, the executable should be kept in `surge/ROMS` if you will be using `make_project.sh`.
4. `maxima/`: This is the directory where the storm surge maxima for each year of simulation will be stored.
5. `plots/`: If you are doing validation, initial plots will be stored here.
6. `analysis/`: Final analysis will go in here.
7. `roms_frc.nc`: ROMS wind and pressure forcing file. This is generated automatically when calling one of `run_synthetic_ensemble.py` or `run_ibtracs_validation.py`.

## How to use
An overview is provided below. However, to get more detailed information on each script you can use the `-h` flag, e.g. `python make_grid.py -h`.

**0. Creating a new project directory**
You can use the `make_project.sh` script in `projects/` to quickly generate a new project directory, automatically copying over the required files. This will create a new directory, the necessary subdirectories and copy the roms template file into the directory. `roms.in.template` and `romsM` will be copied from `CRISK/surge/ROMS` (assuming you have compiled the model). For example:

```
cd projects
make_project.sh woods_hole
```

**1. Making a grid**

`make_grid.py`: Creates a rectangular ROMS grid file from the top corner points. Top corners are defined using -c1 and -c2 (c1 = left, c2 = right) and the vertical length of the domain (-cdist) in degrees. Grid is generated using PyROMS package and Gridgen. Resolution will be approximated as best as possible, but will not be constant throughout the domain. You can specify the approximate resolution in km using the `-r` argument. `-b` is the path of the bathymetry geotiff file to use (from ETOPO or GEBCO) and is specified relative to the project directory (by default it will just look for `projects/<name>/bathy.tiff`. Make sure the data in the bathymetry file completely contains the entire grid you are making. You can download appropriate bathymetry from (GEBCO2023)[https://download.gebco.net/].

For example:

```
python make_grid.py woods_hole  -res 10 -c1 -80 37 -c2 -70.8 43.8 -b bathy.tiff -cdist 10 -hmax 5
```

This will have created a new grid file `roms_grd.nc` in `projects/woods_hole` that has top left corner at (80W, 37N), top right corner at (70.8W, 43.8N), extends ~10 degrees perpendicularly from those corners, has resolution of ~10km and can flood land up to 5m elevation. You can quickly plot this file using:

```
python plot_grid.py woods_hole
```

Which creates `roms_grd.png` in `woods_hole/`. 

**2. Running a synthetic ensemble**

Once you have generated a grid file, you can run an ensemble of simulations with synthetic tropical cyclone using `run_synthetic_ensemble.py`. This will extract all storm that pass within some distance of a point (default 2 degrees), run an individual simulation for each (across multiple cores) and save the maximum surge envelope of point time series for each. To use this script, you will need a STORM track file, which is just a text file containing 1000 years. This should be linked or copied to the project directory with a name like `tracks_IBTRACS.txt`. You can then specify which tracks file to use by including the `-tracks` flag. To run more than 1000 years, concatenate the STORM files into a single file.

For our woods hole example:

```
python woods_hole -ni 4 -nj 4 -slon -74.03 -slat 40.58 -nyears 1000 -tracks IBTRACS -basin NA
```

This command will run 1000 years of STORM synthetic tracks through the model. It will only run tracks that approach within 2 degrees of (74.03W, 40.58N). The simulations will be split across 16 cores (4x4). The script will open a file called tracks_IBTRACS.txt in the project directory.

**3. Running a validation with real storms**

