from climada.hazard import TCTracks
from climada.hazard import Centroids, TropCyclone
import climada.util.coordinates as u_coord
from pathos.pools import ProcessPool as Pool
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from tc_analysis import windspeed_analysis
import os
import argparse
from argparse import RawDescriptionHelpFormatter
import warnings
from time import time
from tc_analysis import plot as tcplot


def main(min_lat, max_lat, min_lon, max_lon, resolution, 
         basin, timestep, n_years,
         ncpus, rperiods, out_dir, tracks_dir,
         pt_lon, pt_lat, pt_name, plot_dir):
    '''
    MAIN function for windspeed_return_periods.py. Makes use of windspeed_analysis.tracks_to_wspd_rp()
    and cliamda TCTracks() and Tropical_Cyclone() classes to analyse windspeeds for each of the 5
    STORM models. Arguments are comparable to those in script header.
    '''

    # Define model names for analysis and file finding
    model_names = ['IBTRACS','HADGEM','ECEARTH','CNRM','CMCC']

    # Create a parallel pool to pass to climada functions
    pool=Pool(ncpus=ncpus)

    # Make centroids and calculate distance to coast
    print('Creating centroids and related parameters..')
    if pt_lon is None or pt_lat is None:
        gridded = True
        cent = Centroids.from_pnt_bounds((min_lon, min_lat, max_lon, max_lat), res=resolution)
    else:
        gridded = False
        cent = Centroids( np.array(pt_lat), np.array(pt_lon) )

    # Calculate distance to land
    cent.set_dist_coast( precomputed=True, signed = False)
    #cent.set_on_land()
    n_cent = len(cent.lon)
    print(f'Total number of centroids: {n_cent} \n')

    # Do return period analysis for each model individually
    ds_out_list = []
    for model in model_names:
        print(f'Analyzing model: {model}')
        
        # Create tracks file path for this loop
        fp_tracks = os.path.join( args.tracks_dir, f'STORM_10000years_{model}_{basin}.txt')
        
        # Read tracks from files and interpolate
        print('   > Reading Tracks')
        tracks = TCTracks.from_simulations_storm( fp_tracks, years=np.arange(0, n_years))
        print('   > Interpolating Tracks')
        tracks.equal_timestep(time_step_h = timestep, pool=pool)

        # Do windspeed analysis
        print('   > Expanding tracks into wind envelopes.')
        ds_rp = windspeed_analysis.tracks_to_wspd_rp( tracks, cent, pool, n_years,
                                                      rperiods = rperiods, reshape_2d = gridded)
        ds_out_list.append(ds_rp)
        del tracks

    # Make output file name
    fp_out = os.path.join(args.out_dir, f'tc_wpsd_returnlevels_{args.name}_{args.nyears}years.nc')

    # Concatenate output datasets
    ds_out = xr.concat(ds_out_list, dim='model')
    ds_future = ds_out.return_level.isel(model=[1,2,3,4])
    ds_out['wspd_ibtracs'] = ds_out.return_level.isel(model=0)
    ds_out['wspd_mean'] = ds_future.median(dim='model')
    ds_out['wspd_median'] = ds_future.mean(dim='model')
    ds_out['wspd_std'] = ds_future.std(dim='model')
    ds_out = ds_out.drop('return_level')

    # Set point names if provided
    if pt_name is not None:
        ds_rp['loc'] = (['loc'], pt_name)

    # Save output dataset to file
    print(f'Saving to file: {fp_out}')
    if os.path.exists(fp_out):
        os.remove(fp_out)
    ds_out.to_netcdf(fp_out)

    # Make plots
    print(f'Saving figures to {plot_dir}')
    if gridded:
        for rp in rperiods:
            f,a = tcplot.compare_windspeed_grid( ds_out.wspd_ibtracs.sel(return_period = rp),
                                                 ds_out.wspd_median.sel(return_period=rp) )
            fp_fig = os.path.join( plot_dir, f'windspeed_grid_{rp}year_{args.name}.png')
            plt.savefig( fp_fig, bbox_inches='tight', dpi=300 )

def check_args(args):
    if args.lat2 <= args.lat1:
        raise ValueError(' lat2 must be greater than or equal to lat1.')
    if args.lon2 <= args.lon1:
        raise ValueError(' lon2 must be greater than or equal to lon1.')
    

if __name__ == "__main__":
    '''
    Calls main() after parsing input arguments.
    '''

    # Silence warnings
    warnings.filterwarnings("ignore")

    description = ('Use STORM synthetic tracks to calculate wind speed return periods associated with Tropical Cyclones.\n'
                  'You can generate these return periods on a regular grid or at set point locations.\n'
                  'Analysis will be saved to a netCDF file with name template tc_wpsd_returnlevels_grid_{name}_{nyears}years.nc. \n'
                  'To generate on a grid, you must specify at least: lat1, lat2, lon1, lon2.\n'
                  'To generate at a set of locations, you must at least provide pt_lon and pt_lat.\n'
                  'If you provide both grid and points, only the points will be analysed.\n'
                  'See arguments and options below for more information.\n\n'
                  'Example Usage: \n \n'
                  '# Generate windspeeds for pt at New York and Boston \n'
                  'python ./analyse_tc_windspeed.py -pt_lon -74 -71.1 -pt_lat 40.71 42.4'
                       ' -pt_name new_york boston -ncpus 12 -name NYBO -basin NA\n\n'
                  '# Generate windspeeds for a grid around New York \n'
                  'python ./analyse_tc_windspeed.py -lat1 40 -lat2 41.5 -lon1 -75 -lon2 -71'
                       ' -ncpus 12 -name NYBO -basin NA\n' )
                     
    # Parse input arguments
    parser = argparse.ArgumentParser(
                    prog='analyse_tc_windspeed',
                    description= description,
                    epilog='Author: David Byrne, Woodwell Climate Research Center',
                    formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument('-basin', type=str, help='Name of STORM basin to use tracks from.', required=True,
                        choices = ['EP','NA','NI','SI','SP','EP', 'WP'])
    parser.add_argument('-lat1', type=float, help='Minimum grid latitude')
    parser.add_argument('-lat2', type=float, help='Maximum grid latitude')
    parser.add_argument('-lon1', type=float, help='Minimum grid longitude')
    parser.add_argument('-lon2', type=float, help='Maximum grid longitude')
    parser.add_argument('-res', type=float, help='Grid resolution in degrees [Default=0.05]', default=0.05)
    parser.add_argument('-pt_lon', type=float, help='Longitude point locations', 
                        default=None, nargs='+')
    parser.add_argument('-pt_lat', type=float, help='Latitude point locations', 
                        default=None, nargs='+')
    parser.add_argument('-pt_name', type=str, help='Names of pt for output file [Default=None]', 
                        default=None, nargs='+')
    parser.add_argument('-name', type=str, default='TEST',
                        help='Name of analysis run. Used in output filename. [Default=TEST]')
    parser.add_argument('-timestep', type=int, help='Timestep of tracks to use. Interpolates to finer track. [Default=1].',
                        default=1)
    parser.add_argument('-rperiods', nargs = '+', type=int,
                        help='Return periods to analyse at [Default=5 10 100 500].',
                        default=[5, 10, 100, 500])
    parser.add_argument('-nyears', type=int, help='Number of synthetic years to use. [Default=10000].', default=10000)
    parser.add_argument('-ncpus', type=int, default=1, 
                        help='Number of cpus to use in Pathos parallel processing [Default=1].')
    parser.add_argument('-tracks_dir', type=str, default='./input/STORM',
                        help='Directory containing STORM concatenated files (10000 years) [Default=./input/STORM].')
    parser.add_argument('-out_dir', type=str, help='Path to output directory [Default=./output]', default='./output')
    parser.add_argument('-plot_dir', type=str, help='Path to plots directory [Default=./plots]', default='./plots')
    args = parser.parse_args()

    
    # Print some things 
    print(' ')
    print(f' *** Starting analysis with {args.nyears} years of simulated storms ***')
    if args.pt_lon is None or args.pt_lat is None:
        print(f' Analyzing on grid with resolution {args.res} degrees:')
        print(f'    Latitude bounds: [{args.lat1}, {args.lat2}]')
        print(f'    Longitude bounds: [{args.lon1}, {args.lon2}]')
    else:
        n_pt = len( args.pt_lon )
        print(f' Analyzing for {n_pt} pt.')
    print(' ******')
    print(' ')

    # Check inputs
    check_args(args)
        
    # Call main function
    tick = time()
    main(args.lat1, args.lat2, args.lon1, args.lon2, args.res, 
         args.basin, args.timestep, args.nyears, args.ncpus, 
         args.rperiods, args.out_dir, args.tracks_dir, 
         args.pt_lon, args.pt_lat, args.pt_name, args.plot_dir)
    tock = time()
    time_elapsed = np.round((tock-tick)/60,2) #Minutes
    print(' ')
    print(f' Done. Total Walltime: {time_elapsed} minutes')
