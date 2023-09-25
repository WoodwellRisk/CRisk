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

description = (' Plots surge envelope ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_grid',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', type=str, help='Path to z envelope file', default='./roms_zenv.nc')
parser.add_argument('-o', type=str, help='Path to output .png file. By default, uses same basename as input file (+.png)', default='surge_envelope.png')
parser.add_argument('-p', type=float, help='Point location to add to plot', default=None, nargs='+')
parser.add_argument('-frc', type=str, help='Path to forcing file if track wanted for plotting', default=None)
parser.add_argument('-cmax', type=float, help=' Maximum for color bar, default auto', default=None )
args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# MAKE PLOTS
#- - - - - - - - - - - - - - - -#

ds = xr.open_dataset(args.f)
minlon = np.min(ds.lon_rho) - 0.05
maxlon = np.max(ds.lon_rho) + 0.05
minlat = np.min(ds.lat_rho) - 0.05
maxlat = np.max(ds.lat_rho) + 0.05

f,a = plt.subplots(1,1, subplot_kw={'projection': ccrs.PlateCarree()},
                       sharey=True, sharex=True, figsize=((6,5)))


gl = a.gridlines(draw_labels=True, linewidth=1, 
                     color='gray', alpha=0.5, linestyle='--')
# manipulate `gridliner` object
gl.xlabels_top = False
gl.ylabels_right = False
a.set_extent([minlon, maxlon, minlat, maxlat], 
                 crs=ccrs.PlateCarree())

a.pcolormesh( ds.lon_rho[1:-1, 1:-1], 
              ds.lat_rho[1:-1, 1:-1], 
              (ds.h < 0)[1:-1,1:-1], cmap=plt.get_cmap('spring'))

# Adjust z so that we mask out all land areas
z = ds.zenv
z = xr.where(ds.h<0, np.nan, z)
if args.cmax is None:
    pc = a.pcolormesh(ds.lon_rho[1:-1, 1:-1], ds.lat_rho[1:-1, 1:-1], 
                      z[1:-1, 1:-1], cmap = plt.get_cmap('Blues',8), vmin=0,)
else:
    pc = a.pcolormesh(ds.lon_rho[1:-1, 1:-1], ds.lat_rho[1:-1, 1:-1], 
                      z[1:-1, 1:-1], cmap = plt.get_cmap('Blues',8), vmin=0, vmax=args.cmax)
a.set_title('Surge Envelope')

if args.p is not None:
    a.scatter(args.p, args.p[1],  marker='o', s=25, facecolors='none', edgecolor='r', linewidth=1.5)

# Plot tracks if present
if args.frc is not None:
    ds_frc = xr.open_dataset(args.frc)
    a.scatter(ds_frc.track_lon, ds_frc.track_lat, c='r', facecolors='none', s=5, linewidths=0.5)
    a.scatter(ds_frc.track_lon, ds_frc.track_lat, c='r', facecolors='none', s=5, linewidths=0.5)

f.subplots_adjust(right=0.9)
cbar_ax = f.add_axes([0.88, 0.15, 0.05, 0.7])
f.colorbar(pc, cax=cbar_ax, extend='max')

f.tight_layout()

# Figure out output name
f.savefig(args.o)