# Coastal Ocean simulation using ROMS

## Overview

## Setup

## Directory Structure

## Functions

**make_grid.py**

Create a rectangular ROMS grid file. Creation of a new grid is done by defining the top corners (c1 = left, c2 = right) and the vertical length of the domain (cdist) in degrees. Call this from your projects directory and specify the project using -proj.

```
optional arguments:
  -h, --help        show this help message and exit
  -proj PROJ        Project name / directory
  -res RES          Approximate resolution of domain in km
  -c1 C1 [C1 ...]   Coordinates of domain top left (lon, lat)
  -c2 C2 [C2 ...]   Coordinates of domain top right (lon, lat)
  -cdist CDIST      How far to project grid into the ocean (deg)
  -b B              Path to bathymetry tiff file to interpolate to model grid
  -rx0_max RX0_MAX  Maximum rx0 to use for bathy smoothing (larger = more smoothing)
  -hmin HMIN        Minimum bathymetric depth. All depths will be clipped to this value.
  -hmax HMAX        Maximum topo value. Elevations greater than this will be masked as a solid wall.
                    If not using wetting and drying, set this to 0.
  -o O              Name of output grid netCDF file
```

**plot_grid.py**

Makes an exploratory plot of a ROMS 2D grid file. Will generate a .png image files containing bathymetry, solid mask and initial wet-dry mask. 

```
optional arguments:
  -h, --help  show this help message and exit
  -proj PROJ  Project name/directory
  -f F        Path to grid file. [Default = roms_grd.nc]
  -o O        Output file name. By default, uses same basename as input file (+.png)
```

**run_ibtracs_validation.py**

 Run ensemble of IBTRACS simulations for a project or list of projects. You can specify a number of run parameters. If you want to compare with observations, you must make sure obs.nc is in the project directory.

```
optional arguments:
  -h, --help            show this help message and exit
  -proj PROJ [PROJ ...]
                        Project name(s) / directory(s)
  -basin BASIN          Name of IBTRACS basin to search for storms
  -fp_tg FP_TG          Name of tidegauge file for validation in project directory
  -fp_val FP_VAL        Name of output validation file
  -o O                  Name of combined validation file when multiple projects are
                        defined
  -fp_grd FP_GRD        Name of grid file in project directory
  -year_start YEAR_START
                        First year to use storms
  -year_end YEAR_END    Last year to use storms
  -ntilei NTILEI        Number of parallel I tiles
  -ntilej NTILEJ        Number of parallel J tiles
  -tau TAU              Wind stress parameterization
  -cdmax CDMAX          Maximum value of stress Cd
  -wind WIND            Wind model
  -inflow_angle INFLOW_ANGLE
                        Constant inflow angle
  -bstress BSTRESS      Bottom stress
  -scale SCALE          Wind Scaling
  -dt DT                Coarse grid delta time
  -nest NEST            1 = nest, 0 = no nest
```
