## Tropical Cyclone Analysis
Routines and scripts for doing analysis of synthetic tropical cyclone tracks.
At the moment, this repository is geared around doing an analysis of STORM synthetic tracks [1,2].
A key library that needs to be installed for many of these scripts is CLIMADA (see here).

Scripts:

1. `analyse_tc_windspeed.py`. Calculate return levels for tropical cyclone windspeeds using the Holland parametric tropical cyclone model. This can be done on a grid or at point locations.
2. `analyse_tc_category.py`. Calculate return periods for each tropical cyclone category passing within 100km. This can be done on a grid or at point locations.

These scripts should be used from the command line using `python analyse_tc_windspeed.py [args]`. You can use the `-h` flag to get more information on how to use each script, arguments and options. 

For example:

`python analyse_tc_windspeed.py -h`


## Setup

### 1. Install Script

The easiest way to set everything up for this package is to use the `install.sh` script.
This script will install pip, install dependencies for the `tc_analysis` scripts and install
climada to your active python environment.

Steps:

* Clone this repository somewhere onto your system using `git clone`.
* Create a new Python environment using `conda` or `mamba` and activate it.
* Change directory into the repository directory (where you cloned it to).
* Run `./install.sh`. You may need to change permissions on the file using `chmod` first.

This scripts will also create `input` and `output` directories. You can keep data wherever you want
by using the arguments in the scripts, however they are setup to look for the above directories by
default.

### 2. Get track data from bucket

Input data required to run these scripts lives at `gs://cmip5_data/tropical_cyclones/`. 
Track data lives at `gs://cmip5_data/tropical_cyclones/STORM/tracks`. These files are the
main input to the scripts in this repository. You can quickly download all of the required
files by calling `./get_data.sh`. This will create a folder called `./input/STORM`.

As mentioned above, by default the scripts in this repository will search for `./input/STORM`.
