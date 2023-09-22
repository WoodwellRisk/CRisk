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

description = (' Makes an exploratory plot of a ROMS 2D grid file ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_grid',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', type=str, help='Path to grid file', default='./roms_grd.nc')
parser.add_argument('-o', type=str, help='Path to output .png file. By default, uses same basename as input file (+.png)', default=None)
parser.add_argument('-p', type=float, help='Point location to add to plot', default=None, nargs='+')
args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# MAKE PLOTS
#- - - - - - - - - - - - - - - -#

ds = xr.open_dataset(args.f)
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

if args.p is not None:
    a[0].scatter(args.p[0], args.p[1],  marker='o', s=25, facecolor='none', edgecolor='r', linewidth=1.5)
    a[1].scatter(args.p[0], args.p[1],  marker='o', s=25, facecolor='none', edgecolor='r', linewidth=1.5)
    a[2].scatter(args.p[0], args.p[1],  marker='o', s=25, facecolor='none', edgecolor='r', linewidth=1.5)
    
f.tight_layout()

# Figure out output name
if args.o is None:
    fp_out = Path(args.f).with_suffix('.png')

f.savefig(fp_out)