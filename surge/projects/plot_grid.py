import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import argparse
import numpy as np
import xarray as xr
import warnings
import os
from pathlib import Path

# Filter warnings
warnings.filterwarnings("ignore")

description = (' Makes an exploratory plot of a ROMS 2D grid file. Will generate a .png image files containing bathymetry, solid mask and initial wet-dry mask. ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_grid',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('proj', type=str, help='Project name/directory')
parser.add_argument('-f', type=str, help='Path to grid file. [Default = roms_grd.nc]', default='roms_grd.nc')
parser.add_argument('-o', type=str, help='Output file name. By default, uses same basename as input file (+.png)', default=None)
args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# MAKE PLOTS
#- - - - - - - - - - - - - - - -#

ds = xr.open_dataset( os.path.join( args.proj, args.f))
minlon = np.min(ds.lon_rho) - 0.05
maxlon = np.max(ds.lon_rho) + 0.05
minlat = np.min(ds.lat_rho) - 0.05
maxlat = np.max(ds.lat_rho) + 0.05

f,a = plt.subplots(1,3, subplot_kw={'projection': ccrs.PlateCarree()},
                       sharey=True, sharex=True, figsize=((10,5)))
a = a.ravel()

for ii in range(3):
    gl = a[ii].gridlines(draw_labels=True, linewidth=1, 
                         color='gray', alpha=0.5, linestyle='--')
    # manipulate `gridliner` object
    gl.xlabels_top = False
    gl.ylabels_right = False
    if ii >0:
        gl.ylabels_left = False
    a[ii].set_extent([minlon, maxlon, minlat, maxlat], 
                     crs=ccrs.PlateCarree())
    a[ii].coastlines(resolution='10m', color=[1, 102/255, 1], linewidth=.5)

a[0].pcolormesh(ds.lon_rho, ds.lat_rho, ds.h, cmap = plt.get_cmap('BrBG',12), vmin=-20, vmax=20)
a[0].set_title('Clipped bathymetry')

a[1].pcolormesh(ds.lon_rho, ds.lat_rho, ds.h>0, cmap='Blues')
a[1].set_title('Positive/Negative bathy')

a[2].pcolormesh(ds.lon_rho, ds.lat_rho, ds.mask_rho, cmap='Reds')
a[2].set_title('Solid Mask')
    
f.tight_layout()

# Figure out output name
if args.o is None:
    fp_out = Path(args.f).with_suffix('.png')
    fp_out = os.path.join( args.proj, fp_out )

f.savefig(fp_out)