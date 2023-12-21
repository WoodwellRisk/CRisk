import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from climada.hazard import TCTracks
import pandas as pd
from shapely.geometry import Point
from crisk_surge import forcing, input_control, validation, nesting, postprocessing
from crisk_surge import plot as rplot
from paratc import track_tools, tctools
import subprocess
import os
import argparse
import warnings
import sys
from glob import glob
warnings.filterwarnings('ignore')

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

n_years = 999
basin = 'NA'
fp_grd = './roms_grd.nc'
proj = 'charleston'
ctr_lon = -79.9
ctr_lat = 32.78
stormmodel = 'IBTRACS'
storm_dir = '../../data/STORM/'

bstress = 1e-3
cdmax = 2.1e-3
cdmodel = 'peng_li15'
scale = 0.89
dt = 6
ntilei = 5
ntilej = 5
ncpu = ntilei * ntilej
min_dist = 3000

fp_tracks = os.path.join( storm_dir, f'STORM_10000years_{stormmodel}_{basin}.txt' )
fp_tracks = os.path.join( storm_dir, 'STORM_DATA_IBTRACS_NA_1000_YEARS_4.txt')

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Change into project directory
os.chdir(proj)
print(f' *****************-> PROJECT: {proj}', flush=True)

# Read tracks between start and end years
print(' Opening synthetic tracks file..', flush=True)
tracks = TCTracks.from_simulations_storm(fp_tracks)
tracks.data = track_tools.shift_tracks_lon(tracks.data) #TODO
sid = [tr.sid for tr in tracks.data]
tracks = [forcing.climada_to_dataframe(ds) for ds in tracks.data]

# Open grid file and get grid polygon
ds_grd = xr.open_dataset( fp_grd )
grid_poly = forcing.get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho)

# Subset tracks in grid domain
keep_idx = track_tools.subset_tracks_in_poly(tracks, grid_poly )
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Subsetting to grid domain: {len(tracks)}', flush=True)

# Filter out storms that don't pass close to central tide gauge
ctr_poly = Point( [ctr_lon, ctr_lat] ).buffer(min_dist)
keep_idx = track_tools.subset_tracks_in_poly(tracks, ctr_poly)
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Keeping only passing storms: {len(tracks)}', flush=True)

# Filter out weak storms
keep_idx = track_tools.filter_tracks_by_column( tracks, col_min = 33 / 0.9 )
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Removed very weak storms: {len(tracks)}', flush=True)

track_years = []
for ii, tr in enumerate(tracks):
    track_years.append( float(sid[ii].split('-')[-2]) )
track_years = np.array(track_years)

n_storms = len(tracks)
print(f' Running an average of {n_storms/1000} per year')

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
        # ds_zmax = postprocessing.make_zero_surge_envelope( )
        # ds_zmax['year'] = yy
        # ds_zmax.to_netcdf( f'./maxima/zmax_y{yystr}.nc') 
        continue
    
    tracks_yy = [tracks[ii] for ii in year_idx]
    sid_yy = [sid[ii] for ii in year_idx]

    z_envelopes = []
    timeseries = []
    
    for ii in range(n_storms_year):

        if os.path.exists('roms_his.nc'):
            os.remove('roms_his.nc')
        
        storm = tracks_yy[ii]
        print(f'Year {yy} / {year_end} --> {ii+1} / {n_storms_year} ---> {sid_yy[ii]}')

        # Make forcing for this storm
        forcing.make_forcing( ds_grd, storm, cdmodel, 1.5, cdmax,
                              './roms_frc.nc', scale )
        
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
        zenvii = postprocessing.calculate_surge_envelope(mask_land=False)
        z_envelopes.append(zenvii)

        # Pull out time series at nearest point
        ds_ts = postprocessing.get_nearest_timeseries_from_file( ctr_lon, ctr_lat, window=48)
        timeseries.append(ds_ts)
        

    ds_zmax = xr.concat(z_envelopes, dim='storm').max(dim='storm')
    ds_zmax['year'] = yy
    ds_zmax.to_netcdf( f'./maxima/zmax_y{yystr}.nc') 

    ds_ts = xr.concat(timeseries, dim='storm')
    ds_ts['year'] = yy
    ds_ts.to_netcdf( f'./maxima/timeseries_y{yystr}.nc') 

# Aggregate
fp_list = glob('./maxima/zmax*')
fp_list.sort()
ds_list = [xr.open_dataset(fp, chunks={}) for fp in fp_list]
ds = xr.concat(ds_list, dim='storm', coords='different').drop('h')
ds['lon_rho'] = ds_grd.lon_rho
ds['lat_rho'] = ds_grd.lat_rho
ds['n_years'] = n_years
fp_zenv = f'./analysis/zenv_{stormmodel}.nc'
if os.path.exists(fp_zenv):
    os.remove(fp_zenv)
ds.to_netcdf(fp_zenv)

fp_list = glob('./maxima/timeseries*')
fp_list.sort()
ds_list = [xr.open_dataset(fp) for fp in fp_list]
ds = xr.concat(ds_list, dim='storm', coords='different').drop('h')
ds['n_years'] = n_years
fp_ts = f'./analysis/timeseries_{stormmodel}.nc'
if os.path.exists(fp_ts):
    os.remove(fp_ts)
ds.to_netcdf(fp_ts)
