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

# Extract the track that intersects with model domain
pol = forcing.get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho )
track = tc_analysis.clip_track_to_poly( track, pol, max_dist=args.buffer )

# Generate pressure field
press = forcing.make_pressure_field( ds_grd.lon_rho.values, 
                                     ds_grd.lat_rho.values,
                                     track.data[0]) 
    
# Generate wind vectors
wind_u, wind_v = forcing.make_windfield( ds_grd.lon_rho.values, 
                                         ds_grd.lat_rho.values, 
                                         track, model=args.model, 
                                         angle=ds_grd.angle.values )

# Convert wind vectors to stresses
tau, sustr_rho, svstr_rho = forcing.tau_andreas(wind_u, wind_v)
#tau2, sustr_rho2, svstr_rho2 = forcing.tau_large_pond(wind_u, wind_v)

# Make dataset from stresses on rho grid
ds_r = xr.Dataset()
ds_r['lon'] = (['eta','xi'], ds_grd.lon_rho.values)
ds_r['lat'] = (['eta','xi'], ds_grd.lat_rho.values)
ds_r['sustr'] = (['sms_time','eta','xi'], sustr_rho)
ds_r['svstr'] = (['sms_time','eta','xi'], svstr_rho)

# Regrid stresses onto U and V grids
ds_u = xr.Dataset()
ds_u['lon'] = (['eta_u','xi_u'], ds_grd.lon_u.values)
ds_u['lat'] = (['eta_u','xi_u'], ds_grd.lat_u.values)
regridder = xe.Regridder(ds_r, ds_u, "bilinear")
ds_u = regridder( ds_r[['sustr']] )

ds_v = xr.Dataset()
ds_v['lon'] = (['eta_v','xi_v'], ds_grd.lon_v.values)
ds_v['lat'] = (['eta_v','xi_v'], ds_grd.lat_v.values)
regridder = xe.Regridder(ds_r, ds_v, "bilinear")
ds_v = regridder( ds_r[['svstr']] )

# Create output dataset and save to file
ds_out = forcing.make_forcing_dataset( ds_grd, track.data[0].time.values,
                                       sustr = ds_u.sustr.values, 
                                       svstr = ds_v.svstr.values, 
                                       press = press,
                                       track_lon = track.data[0].lon.values, 
                                       track_lat = track.data[0].lat.values )

fp_out = os.path.join( args.proj, args.o)
if os.path.exists(fp_out):
    os.remove(fp_out)
ds_out.to_netcdf(fp_out)