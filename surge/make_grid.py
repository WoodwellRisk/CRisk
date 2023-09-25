import os
from pyroms import _iso
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import pygridtools as pgt
import matplotlib.pyplot as plt
import rasterio, rioxarray
import pyroms
import xarray as xr
from bathy_smoother import *
import warnings
import argparse
import sys
import scipy.ndimage as ndi

# Filter warnings
warnings.filterwarnings("ignore")

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

description = (' Create a rectangular ROMS grid file. ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_grid',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-res', type=float, help='Reciprocal of horizontal resolution (30 = 1/30 degree)', default=None)
parser.add_argument('-c1', type=float, help='Coordinates of coastline left (lon, lat)', default=[-93, 32], nargs='+')
parser.add_argument('-c2', type=float, help='Coordinates of coastline right (lon, lat)', default=[-86, 32], nargs='+')
parser.add_argument('-cdist', type=float, help='How far to project grid into the ocean (deg)', default=2)

parser.add_argument('-proj', type=str, help='Projection', default='merc')
parser.add_argument('-b', type=str, help='Path to bathymetry tiff file to interpolate to model grid', default='bathy.tiff')
parser.add_argument('-rx0_max', type=float, help='Maximum rx0 to use for bathy smoothing (larger = more smoothing)', default=0.35)
parser.add_argument('-hmin', type=float, help='Minimum bathymetric depth. All depths will be clipped to this value', default=0)
parser.add_argument('-hmax', type=float, help='Maximum topo value. Elevations greater than this will be masked. If not using wetting and drying, set this to 0.', default=10)
parser.add_argument('-o', type=str, help='Name of output grid netCDF file', default='roms_grd.nc')
parser.add_argument('-grd_name', type=str, help='Name of grid (attributed of output file)', default='TEST')

args = parser.parse_args()

# Get unit normal for extending grid away from coastline
res = 1 / args.res
vec_dx = args.c1[0] - args.c2[0]
vec_dy = args.c1[1] - args.c2[1]
vec_norm = np.array( [-vec_dy, vec_dx] )
vec_norm_unit = vec_norm / np.sqrt( vec_norm[0]**2 + vec_norm[1]**2 )
Mm = int(args.cdist / res)

# Grid corners -- top left, top right, bottom left, bottom right
p_tl = args.c1
p_tr = args.c2
p_bl = p_tl + args.cdist*vec_norm_unit
p_br = p_tr + args.cdist*vec_norm_unit

# If resolution provided, we need to edit the grid bounds to reflect this -- first X direction (Mm)
vec_dx = p_tr[0] - p_tl[0]
vec_dy = p_tr[1] - p_tl[1]
vec_dist = np.sqrt(vec_dx**2 + vec_dy**2)
Lm = int(vec_dist / res )                                  # Number of x points in ROMS speak

# Now Y direction (Lm)
# vec_dx = p_tl[0] - p_bl[0]
# vec_dy = p_tl[1] - p_bl[1]
# vec_dist = np.sqrt(vec_dx**2 + vec_dy**2)
# vec_unit = np.array( [vec_dx, vec_dy] ) / vec_dist
# Mm = int(vec_dist / res )                                  # Number of y points in ROMS speak

print(f'Number of coordinates: X (Xi) = {Lm}, Y (Eta) = {Mm}')

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Get the bounds of the grid and central point
lon_min = min(p_tl[0], p_tr[0], p_bl[0], p_br[0])
lon_max = max(p_tl[0], p_tr[0], p_bl[0], p_br[0])
lon_0 = (lon_min + lon_max) / 2.
lat_min = min(p_tl[1], p_tr[1], p_bl[1], p_br[1])
lat_max = max(p_tl[1], p_tr[1], p_bl[1], p_br[1])
lat_0 = (lat_min + lat_max) / 2.

# Basemap projection
map = Basemap(projection=args.proj, llcrnrlon=lon_min, llcrnrlat=lat_min,
              urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
              resolution='h')

# Inputs to gridgen
lonp = np.array([p_tl[0], p_bl[0], p_br[0], p_tr[0]])
latp = np.array([p_tl[1], p_bl[1], p_br[1], p_tr[1]])
beta = np.array([1, 1, 1, 1])
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

# Interpolate bathymetry
ds_bathy = xr.open_dataset(args.b)
topo = ds_bathy.band_data[0].values
lons = ds_bathy.x.values
lats = ds_bathy.y.values

# depth positive
topo = -topo

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
my_interpolating_function = RegularGridInterpolator((lats, lons), topo, method='linear')
h = my_interpolating_function((hgrd.lat_rho, hgrd.lon_rho))

# Make solid mask and remove lakes
mask_rho = h>-args.hmax
mask_labels = ndi.label(mask_rho)
for ii in range(2, mask_labels[1]):
    mask_rho[ np.where(mask_labels[0] == ii) ] = 0
hgrd.mask_rho = mask_rho

# Minimum depth if > 0
if args.hmin > 0:
    idx = np.where(hgrd.mask_rho <= hmin)
    h[idx] = hmin

# save raw bathymetry
hraw = h.copy()

# Smooth the bathy
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
h = bathy_smoothing.smoothing_Positive_rx0(h>0, h, args.rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# Vertical coordinates, even though not required
theta_b = 2
theta_s = 7.0
Tcline = 50
N = 30
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# Create the pyroms grid object and write to file
grd = pyroms.grid.ROMS_Grid(args.grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename=args.o)
