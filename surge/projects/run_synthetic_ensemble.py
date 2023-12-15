import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from climada.hazard import TCTracks
import pandas as pd
from shapely.geometry import Point
from crisk_surge import forcing, input_control, validation, nesting
from crisk_surge import plot as rplot
from paratc import track_tools, tctools
import subprocess
import os
import argparse
import warnings
import sys
warnings.filterwarnings('ignore')

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

n_years = 250
basin = 'NA'
fp_grd = './roms_grd.nc'
proj = 'new_york'
ctr_lon = -74.04
ctr_lat = 40.58
stormmodel = 'IBTRACS'
storm_dir = '../data/STORM/'
fp_zenv = None

bstress = 1.5e-3
dt = 10
ntilei = 4
ntilej = 4
ncpu = ntilei * ntilej

fp_tracks = os.path.join( storm_dir, f'STORM_10000years_{stormmodel}_{basin}.txt' )

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Change into project directory
os.chdir(proj)
print(f' *****************-> PROJECT: {proj}', flush=True)

# Read tracks between start and end years
print(' Opening best tracks file..', flush=True)
tracks = TCTracks.from_simulations_storm(fp_storm)
tracks.data = track_tools.shift_tracks_lon(tracks.data) #TODO
sid = [tr.sid for tr in tracks.data]
tracks = [track_tools.climada_to_dataframe(ds) for ds in tracks.data]

# Open grid file and get grid polygon
ds_grd = xr.open_dataset( fp_grd )
grid_poly = forcing.get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho)

# Subset tracks in grid domain
keep_idx = track_tools.subset_tracks_in_poly(tracks, grid_poly )
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Subsetting to grid domain: {len(tracks)}', flush=True)

# Filter out storms that don't pass close to central tide gauge
min_dist = 1.5
ctr_poly = Point( [ctr_lon, ctr_lat] ).buffer(min_dist)
keep_idx = track_tools.subset_tracks_in_poly(tracks, tg_poly)
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Keeping only passing storms: {len(tracks)}', flush=True)

# Filter out weak storms
keep_idx = track_tools.filter_tracks_by_column( tracks, col_min = 33 / 0.9 )
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Removed very weak storms: {len(tracks)}', flush=True)

track_years = []
for tr in tracks.data:
    track_years.append( float(sid.split('-')[-2]) )
track_years = np.array(track_years)

n_storms = len(tracks.data)
n_storms


year_list = np.arange(0,n_years)
year_start = 0
year_end = np.max(year_list)
subprocess.run('rm -rf maxima/*', shell=True)

for yy in year_list:

    yystr = str(yy).zfill(4)
    
    # Extract tracks for this year
    year_idx = np.where(track_years == yy)[0]

    # If empty, save empty envelope and continue onto the next year
    n_storms_year = len(year_idx)
    if n_storms_year == 0:
        ds_zmax = postprocessing.make_zero_surge_envelope( )
        ds_zmax['year'] = yy
        ds_zmax.to_netcdf( f'./maxima/zmax_y{yystr}.nc') 
        continue
    
    tracks_yy = TCTracks()
    tracks_yy.data = [tracks.data[ii] for ii in year_idx]

    z_envelopes = []
    
    for ii in range(n_storms_year):

        if os.path.exists('roms_his.nc'):
            os.remove('roms_his.nc')
        
        storm = tracks_yy.data[ii]
        print(f'Year {yy} / {year_end} --> {ii+1} / {n_storms_year} ---> {storm.sid}')

        # Make forcing for this storm
        forcing.make_forcing( ds_grd, storm, args.cdmodel, 1, args.cdmax,
                                './roms_frc.nc', args.scale )
        
        # Make input control file
        input_control.make_infile_from_files( ntilei = ntilei, 
                                              ntilej = ntilej, 
                                              dt= dt, 
                                              bstress = bstress,
                                              fp_frc = 'roms_frc.nc',
                                              fp_grd = 'roms_grd.nc')
    
        # Run model
        subprocess.run( f'mpirun -np {ncpu} --oversubscribe romsM roms.in > log.txt', shell=True,)
    
        # Make annual max
        z_envelopes.append(postprocessing.calculate_surge_envelope())

    ds_zmax = xr.concat(z_envelopes, dim='storm').max(dim='storm')
    ds_zmax['year'] = yy
    ds_zmax.to_netcdf( f'./maxima/zmax_y{yystr}.nc') 

# Aggregate
fp_list = glob('./maxima/zmax*')
fp_list.sort()
ds_list = [xr.open_dataset(fp, chunks={}) for fp in fp_list]

ds = xr.concat(ds_list, dim='storm', coords='different').drop('h')

ds['lon_rho'] = ds_grd.lon_rho
ds['lat_rho'] = ds_grd.lat_rho

if fp_zenv is None:
    ds.to_netcdf(f'./zenv_{stormmodel}.nc')
else:
    ds.to_netcdf( fp_zenv )