import xarray as xr
from . import forcing
from glob import glob
import os 
import numpy as np
from shapely.geometry import Point, LineString
from utide import solve, reconstruct
from datetime import datetime, timedelta
import pandas as pd

def nearest_ocean_points( lon2, lat2, mask, pt_lon, pt_lat ):
    c_index = np.ones_like(lon2)
    r_index = np.ones_like(lon2)
    for ii in range( lon2.shape[0] ):
        r_index[ii] = ii
    for ii in range( lon2.shape[1] ):
        c_index[:,ii] = ii 

    c_index = c_index[ mask==1 ]
    r_index = r_index[ mask == 1]
    mask_noland = mask[ mask == 1 ]
    lon2_noland = np.radians( lon2[ mask == 1 ] )
    lat2_noland = np.radians( lat2[ mask == 1 ] )
    n_pts = len(pt_lon)
    pt_lon_rad = np.radians(pt_lon)
    pt_lat_rad = np.radians(pt_lat)
    min_ind = []
    for ii in range(n_pts):
        dist_pt = haversine( lon2_noland, lat2_noland, 
                                 pt_lon_rad[ii], pt_lat_rad[ii] )
        min_ind.append( np.argmin(dist_pt) )
    return r_index[min_ind].astype(int), c_index[min_ind].astype(int)
    

def haversine(lon1, lat1, lon2, lat2, radians=True):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    if not radians:
        lon1 = np.radians(lon1)
        lon2 = np.radians(lon2)
        lat1 = np.radians(lat1)
        lat2 = np.radians(lat2)
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

def get_tidegauges_in_grid( ds_grd, tg_list, buffer = 0.5 ):

    # Find and open all files
    fp_list = tg_list
    ds_list = [xr.open_dataset(fp) for fp in fp_list]

    # Make lon/lat vector
    lons = np.array([float(ds.lon[0].values) for ds in ds_list])
    lats = np.array([float(ds.lat[0].values) for ds in ds_list])

    # Adjust lons
    lons[lons>180] = lons[lons>180] - 360

    # Make polygon from grid
    poly = forcing.get_grid_poly( ds_grd.lon_rho.values, ds_grd.lat_rho.values )
    exterior_line = LineString(poly.exterior.coords)

    # Get indices of point that are inside the grid poly
    keep_idx = []
    for ii, coord in enumerate(zip(lons, lats)):
        point = Point(coord[0], coord[1])
        poly_dist = point.distance(poly)
        exte_dist = point.distance(exterior_line)

        idx_inside = poly_dist == 0
        idx_enough = exte_dist >= buffer
        if idx_inside and idx_enough:
            keep_idx.append( ii )

    return [ds_list[ii] for ii in keep_idx], [fp_list[ii] for ii in keep_idx]

def surge_envelope_from_file( fp_his = 'roms_his.nc', fp_out = 'roms_zenv.nc'):
    ''' Calculate surge envelope from file, and save to file '''

    # Calculate envelope
    ds_his = xr.open_dataset(fp_his, chunks={'ocean_time':100})
    zenv = calculate_surge_envelope( ds_his )

    # Write to file
    if os.path.exists( fp_out ):
        os.remove(fp_out)
    zenv.to_netcdf(fp_out)

def calculate_surge_envelope( ds_his = None, fp_his = 'roms_his.nc', 
                              fp_out = 'roms_zenv.nc', save_to_file = False,
                              mask_land = True  ):

    if ds_his is None:
        ds_his = xr.open_dataset(fp_his, chunks={'ocean_time':100})
    
    zenv = ds_his.zeta.max(dim='ocean_time')
    zenv = zenv.to_dataset(name='zenv')
    zenv['h'] = ds_his.h

    if mask_land:
        zenv = zenv.where( ds_his.h > 0 )

    if save_to_file:
        if os.path.exists( fp_out ):
            os.remove(fp_out)
        zenv.to_netcdf(fp_out)
    else:
        return zenv

def make_zero_surge_envelope( fp_grd = './roms_grd.nc', ds_grd = None, 
                              fp_out = 'roms_zenv.nc', save_to_file = False,
                              mask_land = True  ):

    if ds_grd is None:
        ds_grd = xr.open_dataset(fp_grd)
    
    zenv = ds_grd[['lon_rho','lat_rho']]
    zenv['zenv'] = (['eta_rho','xi_rho'], np.zeros_like( zenv.lon_rho ) )
    zenv.attrs = {}

    if mask_land:
        zenv = zenv.where( ds_grd.h > 0 )

    if save_to_file:
        if os.path.exists( fp_out ):
            os.remove(fp_out)
        zenv.to_netcdf(fp_out)
    else:
        return zenv
