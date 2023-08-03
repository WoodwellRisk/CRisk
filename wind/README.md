# WINDSPEED ANALYSIS

This folder contains scripts for the analysis Tropical Cyclone windspeed return periods and category return periods.
These analyses directly analyse synthetic cyclone tracks. By default, these scripts are set up to work with tracks form
the STORM datset [1,2].

The scripts in this folder can be run using a command like:

```
python analyse_tc_windspeed.py [args]
```

There are various inputs and arguments that can be provided to each. To see useage examples and input arguments, use the `-h` flag:

```
python analyse_tc_windspeed.py -h
```

### analyse_tc_windspeed.py

Uses STORM synthetic tracks to calculate wind speed return periods associated with Tropical Cyclones. 
You can generate these return periods on a regular grid or at set point locations.
Analysis will be saved to a netCDF file with name template tc_wpsd_returnlevels_grid_{name}_{nyears}years.nc.
The script is run at the command line using something like:

This script uses a parametric approach to determine wind speed return periods due to synthetically generated tropical cyclones. The following analysis is performed:

1. Each track (of 10000 years of tracks) is interpolated in time to a 1 hour frequency to ensure fine resolution track characteristics.
2. For every time step in each track, a 2 dimensional wind field is approximated using a Holland wind profile (Holland et al., 1980).
3. For each track, an envelope of maximum winds is calculated.
4. The envelope from each track is stored for the domain of study and used to determine wind speed return periods.

This script makes use of the CLIMADA library for wind field expansion and return period fitting.

### analyse_tc_category.py

Uses STORM synthetic tracks to determine the annual probability (and return period) of each tropical cyclone category passing with 100km of a set of points or
collection of grid points. The analysis is as follows:

1. Each track (of 10000 years of tracks) is interpolated in time to a 1/2 hour frequency to ensure fine resolution track characteristics.
2. For each time step in each track, the distance from storm center is calculated using a haversine function.
3. Where the distance is less than 100km (or `radius` if specified), the category of the storm is stored for all of those locations.
4. For each point in the study set, the maximum category of the passing storm is stored into a 'maximum category envelope'.
5. Envelopes are used to determine annual passing probabilities for each category.