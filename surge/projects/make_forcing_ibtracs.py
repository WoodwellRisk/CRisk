import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from tc_analysis import windspeed_analysis
from tc_analysis import plot as tcplot
from climada.hazard import TCTracks
from climada.hazard import Centroids, TropCyclone
import tc_analysis
from roms_func import forcing
import pandas as pd
import xesmf as xe
import argparse
import os

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

description = (' Make forcing for a single storm in the IBTRAC database ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_forcing_ibtracs',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-proj', type=str, help='Project name / directory', default = '.')
parser.add_argument('-g', type=str, help='Name of grid file', default='roms_grd.nc')
parser.add_argument('-o', type=str, help='Name of output file', default='roms_frc.nc')
parser.add_argument('-basin', type=str, help='Hurricane basin to search (required if no SID)', default=None)
parser.add_argument('-sid', type=str, help='Hurricane ibtracs ID to use ', default=None)
parser.add_argument('-year', type=int, help='Hurricane season to search for name (required if no SID)')
parser.add_argument('-freq', type=float, help='Timestep frequency of forcing (hours)', default=0.5)
parser.add_argument('-name', type=str, help='Name of hurricane to search for (required if no SID)', default=None)
parser.add_argument('-buffer', type=float, help='Buffer around grid to clip track (degrees)', default=3)
parser.add_argument('-model', type=str, help='Parametric wind model', default='H1980')

args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Read in storm
if args.sid is None:
    track = tc_analysis.read_one_from_ibtracs( year = args.year, name = args.name,
                                               basin = args.basin)
else:
    track = tc_analysis.read_one_from_ibtracs( sid = args.sid )
track.equal_timestep( time_step_h = args.freq )

# Open grid dataset
ds_grd = xr.open_dataset( os.path.join( args.proj, args.g) ).load()

ds_frc = forcing.make_forcing_from_track( ds_grd, track.data[0], rotate_vectors = True,
                                          tau_model = 'peng', wind_model=args.model,
                                          domain_buffer = args.buffer )

fp_out = os.path.join( args.proj, args.o)
if os.path.exists(fp_out):
    os.remove(fp_out)
ds_frc.to_netcdf(fp_out)