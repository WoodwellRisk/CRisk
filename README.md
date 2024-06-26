# Coastal risk assessment using synthetic tropical cyclone tracks

<p align="center">
  <img src="https://github.com/WoodwellRisk/CRisk/blob/main/assets/crisk_ex.png" />
</p>

Bespoke coastal risk assessments play a crucial role in empowering coastal communities to build resilience to the impacts of tropical cyclones. This framework evaluates the magnitude of these risks, going from cyclone tracks to very high resolution projections of flood extent. Here, we have functions and code for:

* **Storm Surge**: Storm surge return level estimation using ROMS.
* **Wind Speed**: Tropical cyclone wind return level estimation using CLIMADA.

## Setup
If you are going to run the ROMS analysis (storm surge), you will need to follow all of the steps below. For just wind speed analysis, steps 1 and 2 are sufficient.

1. **Install Dependencies:** Use env.yml and conda/mamba to create a new Python environment with the necessary dependencies: `conda create -n roms_env -f env.yml`
2. **Install CRISK Python Functions**: The functions are found in `src/crisk`. From the base directory of this repository: `pip install -e .`
3. **Install PyROMS**: To generate ROMs grids, the PyROMS library is used, which requires a separate installation. You can follow [their instructions](https://github.com/ESMG/pyroms) to do this.
4. **Compile ROMS**: Follow the instructions on the [ROMS website](https://www.myroms.org/) to setup your environment to compile and run ROMS. You will need to compile ROMS v3.6, and you can use `ROMS/build_roms.sh` and `ROMS/stormsurge.h`.
5. **Get Data**: You will need STORM track files (concatenated into 10000 year files), a ROMS executable `romsM` (from the step above) and a bathymetry file from [GEBCO2023](https://download.gebco.net/) which covers your whole analysis area.

## Structure

Scripts can be run from anywhere. Prior to running a script, you should create a 'project' directory, which will contain inputs for ROMs and outputs. You can automatically generate this directory using `make_project.sh`:

```
cd scripts
make_project.sh woods_hole
```

This will generate the required subdirectories and copy `roms.in.template` and `romsM` if they exist. Each project directory should contain:

1. `roms_grd.nc` : The ROMS grid file. This is generated using the `make_grid.py` script in `projects/` (see below).  
2. `roms.in.template` : The ROMS template input file. The up to date version of this files is in `surge/ROMS`.
3. `romsM`: The ROMS executable. 
4. `maxima/`: This is the directory where the storm surge maxima for each year of simulation will be stored.
5. `plots/`: If you are doing validation, initial plots will be stored here.
6. `analysis/`: Final analysis will go in here.

## Example Usage
For help on each script, you can use the `-h` flag, e.g. `python make_grid.py -h`.

**1. Making a ROMS grid**

`make_grid.py`: Generates a rectangular ROMS grid file from the top corner points. Top corners are defined using -c1 and -c2 (c1 = left, c2 = right) and the vertical length of the domain (-cdist) in degrees. You can specify the approximate resolution in km using the `-r` argument. To generate a grid for woods hole at a resolution (`-r`) of 10km, length of 10 degrees, a maximum inundation depth of 5m and using `woods_hole/bathy.tiff`:

```
python make_grid.py woods_hole  -r 10 -c1 -80 37 -c2 -70.8 43.8 -b bathy.tiff -cdist 10 -hmax 5
```

This will have created a new grid file `roms_grd.nc` in `projects/woods_hole`. You can quickly plot this file using:

```
python plot_grid.py woods_hole
```

Which creates `roms_grd.png` in `woods_hole/`. 

**2. Running a synthetic ensemble**

You can run an ensemble of simulations to estimate return periods using `run_synthetic_ensemble.py`. This will extract all storm that pass within some distance of a point (default 2 degrees), run an individual simulation for each (across multiple cores) and save the maximum surge envelope of point time series for each. Tracks file should be be linked or copied to the project directory with a name like `tracks_IBTRACS.txt`.

For our woods hole example:

```
cp ../data/STORM/STORM_10000years_CMCC_NA.txt woods_hole/tracks_IBTRACS.txt
python woods_hole -ni 4 -nj 4 -slon -74.03 -slat 40.58 -nyears 1000 -tracks IBTRACS -basin NA
```

This command will run 1000 years of STORM synthetic tracks through the model. It will only run tracks that approach within 2 degrees of (74.03W, 40.58N). The simulations will be split across 16 cores (4x4). The script will open a file called tracks_IBTRACS.txt in the project directory. To run all tracks in the domain, do not specify any of the `-s` variables.

**3. Running a validation with real storms**

Validation is important, and it shoudl be done anywhere that observations are available. You can run an ensemble of real storms through the model using `run_ibtracs_ensemble.py`. To perform the validation, download the relevant tide gauge from the research quality [University of Hawaii Sea Level Center database](https://uhslc.soest.hawaii.edu/network/). This netCDF file should be placed in the project directory and names `obs.nc`. The result will be a file called `validation.csv`, containing a row comparing modelled and observed values for each simulated storm.

```
python run_ibtracs_validation.py woods_hole -ni 4 -nj 4
```

## Methodology Overview

Thousands of years of synthetically (statistically) generated tropical cyclones are expanded into 2D wind and pressure fields using the parametric model of Holland, (1980). At the moment, we use tracks from the STORM dataset (Bloemendaal et al., 2020). Wind fields are then converted to wind stress at the ocean surface using a quadratic function of windspeed and a drag coefficient according to (Peng et al., 2020). A validation against observed storm surges at tide gauges shows a mean absolute error in maximum surge of ~14cm and a correlation of 84%. For more information on the methodology and accuracy assessment, see (methodologydoc).

<p align="center">
  <img src="https://github.com/WoodwellRisk/CRisk/blob/main/assets/git_val.png" />
</p>

## References
Bloemendaal, N., Haigh, I.D., de Moel, H. et al. Generation of a global synthetic tropical cyclone hazard dataset using STORM. Sci Data 7, 40 (2020). https://doi.org/10.1038/s41597-020-0381-2

Bloemendaal et al., A globally consistent local-scale assessment of future tropical cyclone risk.Sci. Adv.8,eabm8438(2022). 

Holland, G. J., 1980: An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Mon. Wea. Rev., 108, 1212–1218


