from climada.hazard import TCTracks
from climada.hazard import Centroids, TropCyclone
from pathos.pools import ProcessPool as Pool
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from tc_analysis import windspeed_analysis
import os
import argparse
import warnings

def main(min_lat, max_lat, min_lon, max_lon, resolution, 
         basin, timestep, n_years,
         ncpus, rperiods, out_dir, tracks_dir,
         points_lon, points_lat):
    '''
    MAIN function for windspeed_return_periods.py. Makes use of windspeed_analysis.tracks_to_wspd_rp()
    and cliamda TCTracks() and Tropical_Cyclone() classes to analyse windspeeds for each of the 5
    STORM models. Arguments are comparable to those in script header.
    '''

    # Define model names for analysis and file finding
    model_names = ['IBTRACS','HADGEM','ECEARTH','CNRM','CMCC']

    # Create a parallel pool to pass to climada functions
    pool=Pool(ncpus=ncpus)

    # Do return period analysis for each model individually
    ds_out_list = []
    for model in model_names:
        print(f'Analyzing model: {model}')
        
        # Create tracks file path for this loop
        fp_tracks = os.path.join( args.tracks_dir, f'STORM_10000years_{model}_{basin}.txt')
        
        # Read tracks from files and interpolate
        print('   > Reading Tracks')
        tracks = TCTracks.from_simulations_storm( fp_tracks, years=np.arange(0, n_years))
        tracks.equal_timestep(time_step_h = timestep, pool=pool)
    
        # Construct centroids (I.E. grid)
        if points_lon is None or points_lat is None:
            print('   > Constructing Centroids for GRID:')
            print(f'Latitude bounds: [{min_lat}, {max_lat}]')
            print(f'Longitude bounds: [{min_lon}, {max_lon}]')
            cent = Centroids.from_pnt_bounds((min_lon, min_lat, max_lon, max_lat), res=resolution)
            cent.check()
        
            # Do windspeed analysis
            print('   > Expanding tracks into wind envelopes.')
            ds_rp = windspeed_analysis.tracks_to_wspd_rp( tracks, cent, pool, n_years,
                                                          rperiods = rperiods)
            # Make output file
            fp_out = os.path.join(args.out_dir, f'./tc_wpsd_returnlevels_grid_{args.name}_{args.nyears}years.nc')
        
        else:
            print('   > Constructing Centroids for POINTS.')
            cent = Centroids( np.array(points_lat), np.array(points_lon) )
            cent.check()
    
            # Do windspeed analysis
            print('   > Expanding tracks into wind envelopes.')
            ds_rp = windspeed_analysis.tracks_to_wspd_rp( tracks, cent, pool, n_years,
                                                          rperiods = rperiods, reshape_2d=False)
            fp_out = os.path.join(args.out_dir, f'./tc_wpsd_returnlevels_points_{args.name}_{args.nyears}years.nc')

        ds_out_list.append(ds_rp)

    # Concatenate output datasets
    ds_out = xr.concat(ds_out_list, dim='model')
    ds_out['model'] = model_names

    # Save output dataset to file
    print(f'Saving to file: {fp_out}')
    if os.path.exists(fp_out):
        os.remove(fp_out)
    ds_out.to_netcdf(fp_out)

if __name__ == "__main__":
    '''
    Calls main() after parsing input arguments.
    '''

    # Silence warnings
    warnings.filterwarnings("ignore")

    description = '''Uses STORM synthetic tracks [1] to calculate wind speed return periods associated with Tropical Cyclones.
                     You can generate these return periods on a regular grid or at set point locations.
                     To generate on a grid, you must specify at least: min_lat, max_lat, min_lon, max_lon
                     To generate at set points, you must at least provide points_lon and points_lat.
                     See arguments and options below for more information.'''
                     
    # Parse input arguments
    parser = argparse.ArgumentParser(
                    prog='calculate_windspeed_return_periods',
                    description= description,
                    epilog='Version 2023.07.19')
    parser.add_argument('-basin', type=str, help='Name of STORM basin to use tracks from', required=True)
    parser.add_argument('-lat1', type=float, help='Minimum grid latitude', default=24)
    parser.add_argument('-lat2', type=float, help='Maximum grid latitude', default=31)
    parser.add_argument('-lon1', type=float, help='Minimum grid longitude', default=275)
    parser.add_argument('-lon2', type=float, help='Maximum grid longitude', default=285)
    parser.add_argument('-res', type=float, help='Grid resolution in degrees [Default=0.05]', default=0.05)
    parser.add_argument('-timestep', type=int, help='Timestep of tracks to use. Interpolates to finer track. [Default=1].',
                        default=1)
    parser.add_argument('-rperiods', nargs = '+', type=int,
                        help='Return periods to analyse at [Default=5 10 20 50 100 200 500].',
                        default=[5, 10, 20, 50, 100, 200, 500])
    parser.add_argument('-nyears', type=int, help='Number of synthetic years to use. [Default=10000].', default=10000)
    parser.add_argument('-ncpus', type=int, default=1, 
                        help='Number of cpus to use in Pathos parallel processing [Default=1].')
    parser.add_argument('-tracks_dir', type=str, default='./input/STORM',
                        help='Directory containing STORM concatenated files (10000 years) [Default=../input/STORM].')
    parser.add_argument('-name', type=str, default='TEST',
                        help='Name of analysis run. Used in output filename: tc_wpsd_returnlevels_{name}_{n_years}years_.nc [Default=TEST]')
    parser.add_argument('-out_dir', type=str, help='Path to output directory [Default=./output]', default='./output')
    parser.add_argument('-points_lon', type=float, help='Longitude point locations [Default = None]', 
                        default=None, nargs='+')
    parser.add_argument('-points_lat', type=float, help='Latitude point locations [Default = None]', 
                        default=None, nargs='+')
    args = parser.parse_args()
    
    # Call main function
    main(args.lat1, args.lat2, args.lon1, args.lon2, args.res, 
         args.basin, args.timestep, args.nyears, args.ncpus, 
         args.rperiods, args.out_dir, args.tracks_dir, 
         args.points_lon, args.points_lat)