import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from climada.hazard import TCTracks
import pandas as pd
from shapely.geometry import Point
from crisk_surge import forcing, input_control, validation
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

description = (' Run ensemble of IBTRACS simulations for a project or list of projects. You can specify a number of run parameters, use the -h flag for more information. If you want to compare with observations, you must make sure obs.nc is in the project directory. This should be a tide gauge file from the UHSLC database, containing sea_level and time variables.')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='rub_ibtracs_validation.py',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('proj', type=str, help='Project name(s) / directory(s)', nargs='+')
parser.add_argument('-basin', type=str, help='Name of IBTRACS basin to search for storms, default=NA', default='NA')
parser.add_argument('-fp_tg', type=str, help='Name of tidegauge file for validation in project directory, Default=obs.nc', 
                    default = 'obs.nc')
parser.add_argument('-o', type=str, help = 'Name of combined validation file when multiple projects are defined, default None', default = None)
parser.add_argument('-fp_grd', type=str, help='Name of grid file in project directory. default=roms_grd.nc', default='roms_grd.nc')
parser.add_argument('-yr_start', type=int, help='First year to use storms. default=2000', default=2000)
parser.add_argument('-yr_end', type=int, help='Last year to use storms. default=2022', default=2022)
parser.add_argument('-ni', type=int, help='Number of parallel I tiles. default=1', default = 1)
parser.add_argument('-nj', type=int, help='Number of parallel J tiles. Default=1.', default = 1)
parser.add_argument('-cdmodel', type=str, help='Wind stress parameterization. Default=garratt77', default='garratt77')
parser.add_argument('-cdmax', type=float, help='Maximum value of stress Cd. Default=0.0035', default=2e-3)
parser.add_argument('-bstress', type=float, help='Bottom stress. Default=2e-3', default = 1e-3)
parser.add_argument('-scale', type=float, help='Wind Scaling. Default=1.', default = .83)
parser.add_argument('-dt', type=float, help='Coarse grid delta time. Default=10', default = 10)
parser.add_argument('-nest', type=int, help='1 = nest, 0 = no nest. Default=0', default = 0)

args = parser.parse_args()

yr_start = args.yr_start
yr_end = args.yr_end
basin = args.basin
fp_tg = args.fp_tg
fp_grd = args.fp_grd
ni = args.ni
nj = args.nj
ncpu = ni * nj

# How many projects are being run
n_proj = len(args.proj)
proj_list = args.proj

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#
completed_proj = []

# Loop over projects
for pp, proj in enumerate( proj_list ):

    # Change into project directory
    os.chdir(proj)
    print(f' *****************-> PROJECT: {proj}', flush=True)
    
    # Open tide gauge file and get time bounds
    ds_tg = xr.open_dataset(fp_tg)
    yr_start = max( pd.to_datetime( ds_tg.time.min().values ).year + 1, yr_start )
    yr_end = min( pd.to_datetime( ds_tg.time.max().values ).year -1, yr_end )
    
    # Read tracks between start and end years
    print(' Opening best tracks file..', flush=True)
    tracks = TCTracks.from_ibtracs_netcdf(basin=basin, year_range=[yr_start,yr_end])
    sid = [tr.sid for tr in tracks.data]
    tracks = [forcing.climada_to_dataframe(ds) for ds in tracks.data]
    
    # Open grid file and get grid polygon
    ds_grd = xr.open_dataset( fp_grd )
    grid_poly = tctools.get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho)
    
    # Subset tracks in grid domain
    print(' Filtering storms...')
    keep_idx = track_tools.subset_tracks_in_poly(tracks, grid_poly )
    tracks = [tracks[ii] for ii in keep_idx]
    sid = [sid[ii] for ii in keep_idx]
    print(f'    Subsetting to grid domain: {len(tracks)}', flush=True)
    
    # Filter out storms that don't pass close to central tide gauge
    min_dist = 1.5
    tg_lon = ds_tg.lon.values[0]
    if tg_lon > 180:
        tg_lon = tg_lon-360
    tg_lat = ds_tg.lat.values[0]
    tg_poly = Point( [tg_lon, tg_lat] ).buffer(min_dist)
    keep_idx = track_tools.subset_tracks_in_poly(tracks, tg_poly)
    tracks = [tracks[ii] for ii in keep_idx]
    sid = [sid[ii] for ii in keep_idx]
    print(f'    Keeping only passing storms: {len(tracks)}', flush=True)
    
    # Filter out weak storms
    keep_idx = track_tools.filter_tracks_by_column( tracks, col_min = 33 / 0.9 )
    tracks = [tracks[ii] for ii in keep_idx]
    sid = [sid[ii] for ii in keep_idx]
    print(f'    Removed very weak storms: {len(tracks)}', flush=True)

    # # Filter bad rmw
    if proj != 'galveston':
        keep_idx = track_tools.filter_bad_rmw( tracks )
        tracks = [tracks[ii] for ii in keep_idx]
        sid = [sid[ii] for ii in keep_idx]
        print(f'    Removing very bad RMW observations: {len(tracks)}', flush=True)

    # Get final number of storms and print, check larger than 0
    n_storms = len(tracks)
    print(f' Running {n_storms} storms.', flush=True)

    if n_storms == 0:
        os.chdir('..')
        continue
    completed_proj.append(proj)
    
    val = []
    reject = []
    year_list = np.arange(yr_start, yr_end+1)
    subprocess.run('rm -rf maxima/*', shell=True)
    subprocess.run('rm -rf plots/*', shell=True)
    
    for yy in year_list:
    
        # Extract tracks for this year
        keep_idx = track_tools.subset_tracks_in_year( tracks, yy )
        tracks_yy = [tracks[ii] for ii in keep_idx]
        sid_yy = [sid[ii] for ii in keep_idx]
    
        # If empty, save empty envelope and continue onto the next year
        n_storms_year = len(tracks_yy)
        if n_storms_year == 0:
            print(f'Year {yy} / {yr_end} --> No tracks found', flush=True)
            continue
        
        for ii in range(n_storms_year):
    
            if os.path.exists('roms_his.nc'):
                os.remove('roms_his.nc')
            
            storm = tracks_yy[ii]
            print(f'Year {yy} / {yr_end} --> {ii+1} / {n_storms_year} ---> {sid_yy[ii]}', flush=True)
            
            forcing.make_forcing( ds_grd, storm, args.cdmodel, 1, args.cdmax,
                                    './roms_frc.nc', args.scale )
            os.chdir('..')
            subprocess.run( f'python ./plot_forcing.py {proj} -o plots/frc_{sid_yy[ii]}.png', shell=True,)
            os.chdir(proj)

            # Check there are enough obs to warrant running the model
            ndata_ha, ndata_frac = validation.check_obs_from_file( fp_tg )
        
            # Don't reject but do flag data
            rejii = False
            if ndata_ha < 3*30*24:
                rejii = True
            if ndata_frac < 0.75:
                rejii = True
            
            # Make input control file
            input_control.make_infile_from_files( ntilei = ni, 
                                                  ntilej = nj, 
                                                  dt=args.dt, 
                                                  bstress=args.bstress,
                                                  fp_frc = 'roms_frc.nc',
                                                  fp_grd = 'roms_grd.nc')
        
            # Run model
            subprocess.run( f'mpirun -np {ncpu} --oversubscribe romsM roms.in > log.txt', shell=True,)

            #################
            # if args.nest > 0:
            #     # Make forcing for this storm
            #     subprocess.run( f'cp roms_his.nc roms_his_c.nc', shell=True)
            #     print('FINE RUN START', flush=True)
            #     ds_bry = nesting.make_boundaries_from_parent( 'roms_grd.nc', 'roms_grd_fine.nc', 
            #                                                   'roms_his.nc', 'roms_bry.nc' )
            #     ds_grd_f = xr.open_dataset('roms_grd_fine.nc')
            #     nesting.interp_to_forcing_time( 'roms_frc_f.nc', 'roms_bry.nc' )
                
            #     # Make input control file
            #     input_control.make_infile_from_files( ntilei = ni, 
            #                                           ntilej = nj, 
            #                                           dt=4, 
            #                                           bstress=args.bstress,
            #                                           fp_out = 'roms_f.in',
            #                                           fp_frc = 'roms_frc_f.nc',
            #                                           fp_grd = 'roms_grd_fine.nc')
        
            #     # Run model
            #     subprocess.run( f'mpirun -np {ncpu} --oversubscribe romsM_fine roms_f.in > log_f.txt', shell=True,)
            #################
        
            # Get validation stats
            ds_v, obs, mod = validation.validate_from_file(fp_tg = fp_tg, sid=sid_yy[ii])
            val.append(ds_v)
            reject.append(rejii)
        
            # Make plots
            f,a = rplot.plot_validation_timseries( mod, obs )
            f.savefig( f'./plots/valts_{sid_yy[ii]}' )
            plt.close('all')

    df_val = pd.concat(val)
    df_val['reject'] = reject
    df_val.to_csv('./validation.csv')
    os.chdir('..')

# If -o is specified, combine validation files
if args.o is not None:
    val_list = [pd.read_csv(f'./{prj}/validation.csv') for prj in completed_proj]
    for ii, vv in enumerate( val_list ):
        vv['proj'] = proj_list[ii]
    df_val = pd.concat(val_list)
    df_val.to_csv(args.o)