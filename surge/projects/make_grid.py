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
from roms_func import postprocessing

# Filter warnings
warnings.filterwarnings("ignore")

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

description = (' Create a rectangular ROMS grid file. Creation of a new grid is done by defining the top corners (c1 = left, c2 = right) and the vertical length of the domain (cdist) in degrees. Call this from your projects directory and specify the project using -proj.')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_grid',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-proj', type=str, help='Project name / directory', default = '.')
parser.add_argument('-res', type=float, help='Approximate resolution of domain in km', default=10)
parser.add_argument('-c1', type=float, help='Coordinates of domain top left (lon, lat)', default=[-93, 32], nargs='+')
parser.add_argument('-c2', type=float, help='Coordinates of domain top right (lon, lat)', default=[-86, 32], nargs='+')
parser.add_argument('-cdist', type=float, help='How far to project grid into the ocean (deg)', default=2)

parser.add_argument('-b', type=str, help='Path to bathymetry tiff file to interpolate to model grid', default='bathy.tiff')
parser.add_argument('-rx0_max', type=float, help='Maximum rx0 to use for bathy smoothing (larger = more smoothing)', default=0.35)
parser.add_argument('-hmin', type=float, help='Minimum bathymetric depth. All depths will be clipped to this value.', default=0)
parser.add_argument('-hmax', type=float, help='Maximum topo value. Elevations greater than this will be masked as a solid wall. If not using wetting and drying, set this to 0.', default=10)
parser.add_argument('-o', type=str, help='Name of output grid netCDF file', default='roms_grd.nc')

args = parser.parse_args()

# Get unit normal for extending grid away from coastline
res = 1 / args.res
vec_dx = args.c1[0] - args.c2[0]
vec_dy = args.c1[1] - args.c2[1]
vec_norm = np.array( [-vec_dy, vec_dx] )
vec_norm_unit = vec_norm / np.sqrt( vec_norm[0]**2 + vec_norm[1]**2 )

# Grid corners -- top left, top right, bottom left, bottom right
p_tl = args.c1
p_tr = args.c2
p_bl = p_tl + args.cdist*vec_norm_unit
p_br = p_tr + args.cdist*vec_norm_unit

# # If resolution provided, we need to edit the grid bounds to reflect this -- first X direction (Mm)
# vec_dx = p_tr[0] - p_tl[0]
# vec_dy = p_tr[1] - p_tl[1]
# vec_dist = np.sqrt(vec_dx**2 + vec_dy**2)
# Lm = int(vec_dist / res )                                  # Number of x points in ROMS speak
# Mm = int(args.cdist / res)

#********
dist_x = postprocessing.haversine( p_tr[0], p_tr[1], p_tl[0], p_tl[1],
                                   radians = False )
dist_y = postprocessing.haversine( p_tl[0], p_tl[1], p_bl[0], p_bl[1],
                                   radians = False )

# Get grid dimensions
Lm = int( dist_x / (args.res) ) + 4
Mm = int( dist_y / (args.res) ) + 4
#*********

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
map = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
              urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
              resolution='h')

# Inputs to gridgen
lonp = np.array([p_tl[0], p_bl[0], p_br[0], p_tr[0]])
latp = np.array([p_tl[1], p_bl[1], p_br[1], p_tr[1]])
beta = np.array([1, 1, 1, 1])
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

for xx,yy in map.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy,np.float32)
    vv = np.zeros((xa.shape[0],2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv,mask_value=0)
tmp_mask = hgrd.mask_rho.copy().astype(bool)

# Get resolution and print
x_res, y_res = postprocessing.get_average_grid_resolution( hgrd.lon_rho, hgrd.lat_rho )
x_res_mean = np.round( np.nanmean(x_res), 2 )
y_res_mean = np.round( np.nanmean(y_res), 2 )
print(f' Mean X resolution: {x_res_mean}km')
print(f' Mean Y resolution: {y_res_mean}km')

# Interpolate bathymetry - Variable name changes accoridng to ETOPO or GEBCO
fp_bathy = os.path.join( args.proj, args.b )
ds_bathy = xr.open_dataset(fp_bathy)
if 'band_data' in ds_bathy:
    topo = ds_bathy.band_data[0].values
elif 'elevation' in ds_bathy:
    topo = ds_bathy.elevation[0].values
lons = ds_bathy.x.values
lats = ds_bathy.y.values

# depth positive
topo = -topo

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
my_interpolating_function = RegularGridInterpolator((lats, lons), topo, method='linear')
h = my_interpolating_function((hgrd.lat_rho, hgrd.lon_rho))

# Make solid mask and remove lakes in solid mask
mask_rho = h>-args.hmax
mask_labels = ndi.label(mask_rho)
for ii in range(2, mask_labels[1]):
    mask_rho[ np.where(mask_labels[0] == ii) ] = 0
hgrd.mask_rho = mask_rho

# Minimum depth if > 0
if args.hmin > 0:
    idx = np.where(hgrd.mask_rho <= hmin)
    h[idx] = hmin

# Clip values on land to be at least 1m
clip_idx = np.logical_and( ~tmp_mask, h>=0 )  
h[clip_idx] = -.25

import matplotlib.pyplot as plt
f,a = plt.subplots(1,1)
a.imshow( clip_idx )
f.savefig('test.png')

# save raw bathymetry
hraw = h.copy()

# Create the pyroms grid object and write to file
grdname = args.proj
if args.proj == '.':
    grdname = 'NONAME'
grd = pyroms.grid.ROMS_Grid(grdname, hgrd)
grd.vgrid.h = h
grd.vgrid.hraw=hraw
fp_out = os.path.join(args.proj, args.o)
pyroms.grid.write_ROMS_grid(grd, fp_out)

# Save some data to grid_info file
fp_info = os.path.join( args.proj, 'grid_info.txt' )
fid = open(fp_info, 'w')
fid.write('ROMS grid created using make_grid.py \n')
fid.write('')
fid.write(f' c1      :: {args.c1}\n')
fid.write(f' c2      :: {args.c2}\n')
fid.write(f' cdist   :: {args.cdist}\n')
fid.write(f' res     :: {args.res}\n')
fid.write(f' rx0_max :: {args.rx0_max}\n')
fid.close()
