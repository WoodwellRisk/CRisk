import xarray as xr
import rioxarray
import matplotlib.pyplot as plt
import numpy as np
from climada.hazard import TCTracks
import pandas as pd
from shapely.geometry import Point
from crisk.surge import forcing, input_control, validation, nesting, postprocessing
from crisk.surge import plot as rplot
import crisk.utils as utils
from paratc import track_tools, tctools
import subprocess
import os
import argparse
import warnings
import sys
import geopandas as gpd
from glob import glob
warnings.filterwarnings('ignore')

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

description = (' Run ensemble of synthetic simulations for a project. You can specify a number of run parameters. If you want to compare with observations, you must make sure obs.nc is in the project directory.')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='run_synthetic_ensemble.py',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('proj', type=str, help='Project name(s) / directory(s)')
parser.add_argument('-nyears', type=int, help='Last year to use storms. default=2022', default=100)
parser.add_argument('-ni', type=int, help='Number of parallel I tiles. default=1', default = 1)
parser.add_argument('-nj', type=int, help='Number of parallel J tiles. Default=1.', default = 1)
parser.add_argument('-slon', type=float, help='Longitude of points around which to contruct track neighbourhood.', default=None, nargs='+')
parser.add_argument('-slat', type=float, help='Latitude of points around which to contruct track neighbourhood.', default=None, nargs='+')
parser.add_argument('-sfile', type=str, help='Name of shapefile containing track neighbourhood polygon', default=None)
parser.add_argument('-sbuff', type=float, help='Buffer width (degrees) around track neighbourhood to use tracks', default=2)
parser.add_argument('-tracks', type=str, help='Name of STORM synthetic track model to use in filename construction. Will also be used in names of analysis files.', default='IBTRACS')
parser.add_argument('-cdmodel', type=str, help='Wind stress parameterization. Default=peng_li15', default='peng_li15')
parser.add_argument('-cdmax', type=float, help='Maximum value of stress Cd. Default=2.1e-3', default=2.1e-3)
parser.add_argument('-bstress', type=float, help='Bottom stress. Default=2e-3', default = 1e-3)
parser.add_argument('-scale', type=float, help='Wind Scaling. Default=0.89', default = .88)
parser.add_argument('-dt', type=float, help='Coarse grid delta time. Default=10', default = 10)

args = parser.parse_args()

fp_grd = './roms_grd.nc'
ncpu = args.ni * args.nj

fp_tracks = f'./tracks_{args.tracks}.txt'

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Change into project directory
os.chdir(args.proj)
print(f' **********-> PROJECT: {args.proj} <-**********', flush=True)
print(f' **********-> Running {args.nyears} years of STORM_{args.tracks} simulations on {ncpu} CPUs <-**********')

# Read tracks between start and end years
print(' Processing synthetic tracks..', flush=True)
tracks = TCTracks.from_simulations_storm( os.path.abspath(fp_tracks) )
tracks.data = track_tools.shift_tracks_lon(tracks.data)
sid = [tr.sid for tr in tracks.data]
tracks = [utils.climada_to_dataframe(ds) for ds in tracks.data]
print(f'    Total number of tracks in file: {len(tracks)}', flush=True)

# Remove tracks with 1 point
track_lens = np.array([len(tr) for tr in tracks])
keep_idx = np.where( track_lens > 1 )[0]
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Removing single point tracks: {len(tracks)}', flush=True)

# Open grid file and get grid polygon
ds_grd = xr.open_dataset( fp_grd )
grid_poly = forcing.get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho)

# Subset tracks in grid domain
keep_idx = track_tools.subset_tracks_in_poly(tracks, grid_poly, buffer=2 )
tracks = [tracks[ii] for ii in keep_idx]
sid = [sid[ii] for ii in keep_idx]
print(f'    Subsetting to grid domain: {len(tracks)}', flush=True)

# Filter out storms that don't pass close to central tide gauge
spoly = None
if args.slon is not None and args.slat is not None:
    points = [Point(args.slon[ii], args.slat[ii]).buffer(args.sbuff) for ii in range(len(args.slon))]
    gpd_poly = gpd.GeoDataFrame( geometry = points ).dissolve()
    spoly = gpd_poly.geometry[0]
elif args.sfile is not None:
    spoly = gpd.read_file(args.sfile).buffer(args.sbuff).geometry[0]
    
if spoly is not None:
    keep_idx = track_tools.subset_tracks_in_poly(tracks, spoly)
    tracks = [tracks[ii] for ii in keep_idx]
    sid = [sid[ii] for ii in keep_idx]
    print(f'    Keeping only storms in study area / near points: {len(tracks)}', flush=True)

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
if n_storms == 0:
    raise Exception('No storms found that meet criteria')
print(f'    Running an average of {n_storms/args.nyears} per year')

year_list = np.arange(0,args.nyears)
year_start = 0
year_end = np.max(year_list)
subprocess.run('rm -rf maxima/*', shell=True)
print(' Running ROMS simulations...')

for yy in year_list:

    yystr = str(yy).zfill(4)
    
    # Extract tracks for this year
    year_idx = np.where(track_years == yy)[0]

    # If empty, save empty envelope and continue onto the next year
    n_storms_year = len(year_idx)
    if n_storms_year == 0:
        continue
    
    tracks_yy = [tracks[ii] for ii in year_idx]
    sid_yy = [sid[ii] for ii in year_idx]

    z_envelopes = []
    timeseries = []
    
    for ii in range(n_storms_year):

        if os.path.exists('roms_his.nc'):
            os.remove('roms_his.nc')
        
        storm = tracks_yy[ii]
        print(f'    Year {yy} / {year_end} --> {ii+1} / {n_storms_year} ---> {sid_yy[ii]}')

        # Make forcing for this storm
        try:
            forcing.make_forcing( ds_grd, storm, args.cdmodel, 0.1, args.cdmax,
                                  './roms_frc.nc', args.scale )
        except:
            print('Failed to make forcing')
            continue
        
        # Make input control file
        input_control.make_infile_from_files( ntilei = args.ni, 
                                              ntilej = args.nj, 
                                              dt= args.dt, 
                                              bstress = args.bstress,
                                              fp_frc = 'roms_frc.nc',
                                              fp_grd = 'roms_grd.nc')
    
        # Run model
        subprocess.run( f'/usr/bin/mpirun -np {ncpu} --oversubscribe romsM roms.in > log.txt', shell=True,)

        # Check that run was a success (ERROR will be in log file if not)
        if open('log.txt', 'r').read().find('ERROR') >=0:
            raise Exception(f' ERROR: Model crashed at storm {sid}. See {args.proj}/log.txt')
    
        # Make annual max
        zenvii = postprocessing.calculate_surge_envelope(mask_land=False)
        zenvii['year'] = yy
        z_envelopes.append(zenvii)

        # Pull out time series at nearest point
        ds_his = xr.open_dataset('roms_his.nc')
        ts_yy= postprocessing.get_coastal_data( ds_his.h, ds_his[['zeta']] )
        ts_yy = postprocessing.align_timeseries_by_max( ts_yy.zeta )
        ts_yy['year'] = yy
        timeseries.append(ts_yy)

    if len(z_envelopes) > 0:
        ds_zmax = xr.concat(z_envelopes, dim='storm').max(dim='storm')
        ds_zmax.to_netcdf( f'./maxima/zmax_y{yystr}.nc') 
        ds_ts = xr.concat(timeseries, dim='storm')
        ds_ts.to_netcdf( f'./maxima/timeseries_y{yystr}.nc') 
        

# Aggregate zmax
fp_list = glob('./maxima/zmax*')
fp_list.sort()
ds_list = [xr.open_dataset(fp, chunks={}) for fp in fp_list]
ds = xr.concat(ds_list, dim='storm', coords='different').drop('h')
ds['lon_rho'] = ds_grd.lon_rho
ds['lat_rho'] = ds_grd.lat_rho
ds['nyears'] = args.nyears
fp_zenv = f'./analysis/zenv_{args.tracks}.nc'
if os.path.exists(fp_zenv):
    os.remove(fp_zenv)
ds.to_netcdf(fp_zenv)

# Aggregate time series
fp_list = glob('./maxima/timeseries*')
fp_list.sort()
ds_list = [xr.open_dataset(fp) for fp in fp_list]
ds = xr.concat(ds_list, dim='storm', coords='different')
ds['nyears'] = args.nyears
fp_ts = f'./analysis/timeseries_{args.tracks}.nc'
if os.path.exists(fp_ts):
    os.remove(fp_ts)
ds.to_netcdf(fp_ts)

# POST ANALYSIS: Mean surge profiles
print('Calculating mean surge profiles.')
ts_mean = postprocessing.analyse_mean_surge( ds )
ts_mean.to_netcdf(f'./analysis/mean_surge_{args.tracks}.nc')

# POST ANALYSIS: Return periods
print('Calculating Return Periods.')
fp_rp = f'./analysis/return_periods_grid_{args.tracks}.nc'
fp_coast = f'./analysis/return_periods_coastal_{args.tracks}.gpkg'
postprocessing.analyse_return_periods(fp_zenv = fp_zenv, 
                                      fp_grd = 'roms_grd.nc',
                                      fp_out_rp = fp_rp,
                                      fp_out_coast = fp_coast)
### Make tifs
# ds = xr.open_dataset( f'./analysis/return_periods_grid_{args.tracks}.nc' ) 
# grd = xr.open_dataset('roms_grd.nc')
# rp = ds.rp.values

# # Interpolate to regular grid
# lonbnds = ( ds.lon_rho.min().values, ds.lon_rho.max().values )
# latbnds = ( ds.lat_rho.min().values, ds.lat_rho.max().values )

# # Create regular grid
# grid_lon = np.arange( lonbnds[0], lonbnds[1], r/222 )
# grid_lat = np.arange( latbnds[0], latbnds[1], r/222 )
# grid = xr.Dataset( coords = dict( lon = grid_lon, lat=grid_lat ) )

# # Interpolate
# ds = ds.rename({'lon_rho':'lon','lat_rho':'lat'})
# rg = xe.Regridder( ds, grid, method = 'nearest_s2d')
# h_rg = rg(grd.h)
# grid['mask'] = h_rg > 0
# ds['mask'] = grd.h > 0
# rg = xe.Regridder( ds, grid, method = 'bilinear', 
#                    extrap_method='nearest_s2d')
# ds_rg = rg( ds )
# ds_rg = ds_rg.rename({'lat':'y','lon':'x'})

# for ii, rpii in enumerate(rp):
#     dsii = ds_rg.sel(rp = rpii).surge
#     dsii = dsii.rio.write_crs(4326)
#     dsii.rio.to_raster(f'surge_rp_{rpii}yr.tiff')
