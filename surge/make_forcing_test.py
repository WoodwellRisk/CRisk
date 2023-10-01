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
from datetime import datetime
import xesmf as xe
import argparse
import os

#- - - - - - - - - - - - - - - -#
# HANDLE INPUTS
#- - - - - - - - - - - - - - - -#

description = (' Make test forcing for testing domains ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_test_forcing',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-proj', type=str, help='Project name/directory', default='.')
parser.add_argument('-g', type=str, help='Name of grid file', default='roms_grd.nc')
parser.add_argument('-o', type=str, help='Name of output file', default='roms_frc.nc')
parser.add_argument('-u', type=float, help='Direction of u winds (relative to grid)', default=0)
parser.add_argument('-v', type=float, help='Direction of v winds (relative to grid)', default=0)
parser.add_argument('--rotate', type=bool, help='Rotate stress vectors to be relative to N/E', action=argparse.BooleanOptionalAction,
                    default=False)

args = parser.parse_args()
print(args)

#- - - - - - - - - - - - - - - -#
# SCRIPT
#- - - - - - - - - - - - - - - -#

# Make dates and open grid file
dates = pd.date_range( datetime(2010,1,1), datetime(2010,1,31), freq='1H' )
ds_grd = xr.open_dataset( os.path.join( args.proj, args.g) )

# Generate winds and stresses
if args.rotate:
    angle = ds_grd.angle.values
else:
    angle = None
wind_u, wind_v =  forcing.make_uniform_windfield( ds_grd.lon_rho, ds_grd.lat_rho, 
                                                  dates, args.u, args.v, 
                                                  angle=angle)
tau, sustr_rho, svstr_rho = forcing.tau_andreas( wind_u, wind_v)

# Generate pressure
xy_shape = ds_grd.lon_rho.shape
press = np.zeros( (len(dates), xy_shape[0], xy_shape[1]) )
for ii in range(xy_shape[0]):
    press[:,ii] = 1000-3*ii

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
ds_out = forcing.make_forcing_dataset( ds_grd, time = dates,
                                       sustr = ds_u.sustr.values, 
                                       svstr = ds_v.svstr.values, 
                                       press = press)

fp_out = os.path.join( args.proj, args.o)
if os.path.exists(fp_out):
    os.remove(fp_out)
ds_out.to_netcdf(fp_out)