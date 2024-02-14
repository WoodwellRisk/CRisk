import numpy as np
from climada.hazard.tc_tracks import TCTracks, estimate_rmw
from climada.hazard import Centroids, TropCyclone
from climada.hazard import trop_cyclone
import climada.util.coordinates as u_coord
from shapely.geometry import Polygon, Point, LineString
import xarray as xr
import pandas as pd
from . import postprocessing
from paratc import track_tools
import os

def get_forcing_info( fp_frc = './roms_frc.nc' ):

    ds_frc = xr.open_dataset(fp_frc)
    min_time = pd.to_datetime( ds_frc.sms_time.min().values )
    max_time = pd.to_datetime( ds_frc.sms_time.max().values )
    n_seconds = (max_time - min_time).total_seconds()
    return min_time, max_time, n_seconds

def get_grid_poly( lon, lat ):

    c1 = (lon[0,0], lat[0,0])      # Top left
    c2 = (lon[0,-1],lat[0,-1])     # Top right
    p1 = (lon[-1,0], lat[-1,0])    # Bottom left
    p2 = (lon[-1,-1], lat[-1,-1])  # Bottom right
    pol = Polygon( (c1, c2, p2, p1) )
    return pol

def rotate_winds( x,y, angle ):
    ''' Rotate wind vectors counter clockwise '''

    cosang = np.cos(angle)
    sinang = np.sin(angle)
    new_x = x * cosang - y * sinang
    new_y = x * sinang + y * cosang

    return new_x, new_y

def make_uniform_windfield( lon, lat, time, u, v, angle=None ):

    if hasattr( u, '__len__'):
        wind_u = np.ones( (len(time), lon.shape[0], lon.shape[1]) )
        wind_v = np.ones( (len(time), lon.shape[0], lon.shape[1]) )
        for ii in range(len(u)):
            wind_u[ii] = u[ii]
            wind_v[ii] = v[ii]
    else:
        wind_u = np.ones( (len(time), lon.shape[0], lon.shape[1]) )*u
        wind_v = np.ones( (len(time), lon.shape[0], lon.shape[1]) )*v

    if angle is not None:
        wind_u, wind_v = rotate_winds(wind_u, wind_v, -np.radians(angle))
    
    return wind_u, wind_v

def make_forcing(   ds_grd, track, 
                    cd_model = 'peng_li15', domain_buffer = 1, 
                    cdmax = 2.5e-3, fp_out = 'roms_grd.nc', 
                    scale_winds = 0.81 ):

    from paratc.tc_models import Holland1980 as h80

    # Extract the track that intersects with model domain
    pol = get_grid_poly( ds_grd.lon_rho, ds_grd.lat_rho )
    track = track_tools.clip_track_to_poly( track, pol, max_dist=domain_buffer)

    # Create storm instance and look at dataset
    storm = h80( track, ds_grd.lon_rho.values,
                 ds_grd.lat_rho.values, B_model='powell05', 
                 interp_timestep=.1, rmw_model = 'VW08' )
    storm.scale_winds( scale_winds )
    storm.apply_inflow_angle( inflow_model = 'nws' )
    storm.add_background_winds( bg_model = 'constant', bg_alpha = .55, bg_beta = 20 )
    storm.make_wind_stress( cd_model, cd_max=cdmax )
    storm.to_ROMS(ds_grd)
    storm.to_netcdf( fp_out )