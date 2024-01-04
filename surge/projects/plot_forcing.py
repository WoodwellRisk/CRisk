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
                prog='plot_forcing',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('proj', type=str, help='Project name/directory', default='.')
parser.add_argument('-f', type=str, help='Name of forcing file', default='roms_frc.nc')
parser.add_argument('-o', type=str, help='Name of output file. By default, uses same basename as input file (+.png)', default='forcing.png')
args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# MAKE PLOTS
#- - - - - - - - - - - - - - - -#

ds = xr.open_dataset( os.path.join( args.proj, args.f) )
minlon = np.min(ds.lon_rho) - 0.05
maxlon = np.max(ds.lon_rho) + 0.05
minlat = np.min(ds.lat_rho) - 0.05
maxlat = np.max(ds.lat_rho) + 0.05

f,a = plt.subplots(1,2, subplot_kw={'projection': ccrs.PlateCarree()},
                       sharey=True, sharex=True, figsize=((8,5)))
a = a.ravel()

for ii in range(2):
    gl = a[ii].gridlines(draw_labels=True, linewidth=1, 
                         color='gray', alpha=0.5, linestyle='--')
    # manipulate `gridliner` object
    gl.xlabels_top = False
    gl.ylabels_right = False
    if ii >0:
        gl.ylabels_left = False
    a[ii].set_extent([minlon, maxlon, minlat, maxlat], 
                     crs=ccrs.PlateCarree())
    a[ii].coastlines(resolution='10m', color='b', linewidth=.5)

s = np.sqrt( ds.sustr.values[:,:-1,]**2 + ds.svstr.values[:,:,:-1]**2)
s_env = np.max(s, axis=0)
pc_s = a[0].pcolormesh(ds.lon_rho, ds.lat_rho, s_env, cmap = plt.get_cmap('Reds',8), vmin=0,)
a[0].set_title('Wind Stress Envelope')
f.colorbar(pc_s, ax=a[0], fraction=0.05)

#a[1].pcolormesh(ds.lon_rho, ds.lat_rho, ds.h>0, cmap='Blues')
a[1].set_title('Pressure Envelope')
pmin = np.min( ds.Pair.values, axis=0 )
pc_p = a[1].pcolormesh(ds.lon_rho, ds.lat_rho, pmin, cmap = plt.get_cmap('Greens_r',8))
f.colorbar(pc_p, ax=a[1], fraction=0.05)

# Plot tracks if present
if 'track_lon' in ds and 'track_lat' in ds:
    a[0].scatter(ds.track_lon, ds.track_lat, c='k', facecolors='none', s=1.5, linewidths=0.5)
    a[1].scatter(ds.track_lon, ds.track_lat, c='k', facecolors='none', s=1.5, linewidths=0.5)
    
f.tight_layout()

# Figure out output name
fp_out = os.path.join( args.proj, args.o )
f.savefig(fp_out)