## Functions

**make_grid.py**
```
 Create a rectangular ROMS grid file. Creation of a new grid is done by defining the top corners (c1 = left, c2 = right) and the vertical length of the domain (cdist) in degrees. Call this from your projects directory and specify the project using -proj.

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
  -hmax HMAX        Maximum topo value. Elevations greater than this will be masked as a solid wall. If not using wetting and drying,
                    set this to 0.
  -o O              Name of output grid netCDF file
```
