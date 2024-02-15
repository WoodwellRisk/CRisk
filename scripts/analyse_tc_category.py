import pandas as pd
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import dask
import dask.delayed as delayed
from dask.distributed import Client
from glob import glob
from crisk.wind import radial_analysis
import multiprocessing
import rioxarray
from climada.hazard import TCTracks
from crisk import utils
from paratc import track_tools, tctools
import warnings
import argparse
from argparse import RawDescriptionHelpFormatter

if __name__ == '__main__':

    # Make client
    client = Client()

    # Silence warnings
    warnings.filterwarnings("ignore")
    
    description = ('')
                     
    # Parse input arguments
    parser = argparse.ArgumentParser(
                    prog='analyse_tc_category',
                    description= description,
                    epilog='Author: David Byrne, Woodwell Climate Research Center',
                    formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument('-basin', type=str, help='Name of STORM basin to use tracks from.', required=True,
                        choices = ['EP','NA','NI','SI','SP','WP'])
    parser.add_argument('-lon', type=float, help='Grid longitude (min, max)', nargs='+')
    parser.add_argument('-lat', type=float, help='Grid latitude (min, max)', nargs='+')
    parser.add_argument('-res', type=float, help='Grid resolution in degrees [Default=0.05]', default=0.05)
    parser.add_argument('-name', type=str, default='TEST',
                        help='Name of analysis run. Used in output filename. [Default=TEST]')
    parser.add_argument('-timestep', type=int, help='Timestep of tracks to use. Interpolates to finer track. [Default=1].',
                        default=.25)
    parser.add_argument('-nyears', type=int, help='Number of synthetic years to use. [Default=10000].', default=5000)
    parser.add_argument('-ncpus', type=int, default=1, 
                        help='Number of cpus to use in Pathos parallel processing [Default=1].')
    parser.add_argument('-tracks_dir', type=str, default='../data/STORM',
                        help='Directory containing STORM concatenated files (10000 years) [Default=./input/STORM].')
    parser.add_argument('-radius', type=str, help='Storm Radius', default=100)
    parser.add_argument('-model', type=str, help='Track file names, or list of', nargs='+', default=None)
    outdir = '.'
    args = parser.parse_args()
    margin = 3
    
    # Set default model names if not provided
    if args.model is None:
        model_list = ['IBTRACS','HADGEM','ECEARTH','CMCC','CNRM']
    else:
        model_list = args.model
        
    if type(model_list) is not list:
        model_list = list(model_list)
    
    # Make the grid
    grid_lon = np.arange( args.lon[0], args.lon[1], args.res )
    grid_lat = np.arange( args.lat[0], args.lat[1], args.res )
    lon2, lat2 = np.meshgrid( grid_lon, grid_lat )
    grid_poly = utils.get_grid_poly( lon2, lat2 )

    # Loop over models in sepcified list
    ds_list = []
    for model_name in model_list:
    
        # Get list of inputs and read into dataframe
        fp_in = f'{args.tracks_dir}/STORM_10000years_{model_name}_{args.basin}.txt'
        print(f' Opening... {fp_in}')
        tracks = TCTracks.from_simulations_storm( fp_in )
        tracks.data = tracks.data
        tracks.data = track_tools.shift_tracks_lon(tracks.data)
        sid = [tr.sid for tr in tracks.data]
        tracks = [utils.climada_to_dataframe(ds) for ds in tracks.data]

        # Get years of STORM tracks
        track_years = []
        for ii, tr in enumerate(tracks):
            track_years.append( float(sid[ii].split('-')[-2]) )
        track_years = np.array(track_years)

        # Take first n years of tracks
        keep_idx = np.where(track_years < args.nyears)[0]
        tracks = [tracks[ii] for ii in keep_idx]
        sid = [sid[ii] for ii in keep_idx]
        track_years = track_years[keep_idx]

        # Remove 1 point tracks
        tlens = np.array([len(tr) for tr in tracks])
        keep_idx = np.where( tlens > 1 )[0]
        tracks = [tracks[ii] for ii in keep_idx]
        sid = [sid[ii] for ii in keep_idx]
        track_years = track_years[keep_idx]
    
        # Remove tracks outside of the grid domain
        keep_idx = track_tools.subset_tracks_in_poly(tracks, grid_poly, buffer=2 )
        tracks = [tracks[ii] for ii in keep_idx]
        sid = [sid[ii] for ii in keep_idx]
        track_years = track_years[keep_idx]
        print(f'    Subsetting to grid domain: {len(tracks)}', flush=True)

        # Do the analysis for this model
        ds_pasana = radial_analysis.passing_storm_stats( tracks, lon2, lat2, track_years,
                                                         client=client, radius=args.radius)

        # Format output and save to file
        ds_pasana['annual_count'] = args.nyears / ds_pasana['annual_count']
        ds_pasana['event_count'] = ds_pasana['event_count'] / args.nyears
        ds_pasana = ds_pasana.rename({'annual_count':'return_period',
                        'event_count':'events_per_year'})
        ds_pasana.attrs = {'name':f'Hurricane frequencies within {args.radius}km according to the STORM dataset',
                    'nyears':f'{args.nyears}',
                    'model':f'{model_name}'}
        ds_pasana.return_period.attrs = {'dscription':f'Annual return period of storm exceeding bin intensity'}
        ds_pasana.events_per_year.attrs = {'description':f'Average number of events per year exceeding bin intensity'}
        ds_pasana.bin.attrs = {'description':'Hurricane intensity bins',
                        'units':'ms-1'}
        ds_pasana['x'] = (['x'], grid_lon)
        ds_pasana['y'] = (['y'], grid_lat)
        ds_pasana.to_netcdf(f'passing_hurricanes_{args.name}_{model_name}_{args.radius}km.nc')