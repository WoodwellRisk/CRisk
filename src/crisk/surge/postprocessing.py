import xarray as xr
from . import forcing
from glob import glob
import os 
import numpy as np
from shapely.geometry import Point, LineString
from utide import solve, reconstruct
from datetime import datetime, timedelta
import pandas as pd
import scipy.ndimage as ndi
import scipy.stats as stats
from scipy.interpolate import interp1d
import geopandas as gpd

def analyse_mean_surge( ts ):
    ''''''
    n_storm, n_time, n_index = ts.zeta.shape
    # Get percentiles of maxima
    pc_bins = np.arange(80, 105, 5)/100
    n_bins = len(pc_bins) - 1
    mean_surge = np.zeros(( n_index, n_bins, n_time ))
    tsmax = ts.max(dim='time')
    tsmax_pc = tsmax.quantile(dim='storm', q=pc_bins)
    
    for ii in range( n_index ):
        tsmax_ii = tsmax.zeta.isel(index = ii)
        z_ii = ts.zeta.isel( index=ii ) / tsmax_ii
        pc_ii = tsmax_pc.zeta.isel(index=ii).values
        for bb in range(n_bins):
            bb_idx = np.logical_and( tsmax_ii >= pc_ii[bb], tsmax_ii < pc_ii[bb+1] )
            bb_idx = np.where( bb_idx.values )[0]
            z_ii_bin = z_ii.isel( storm = bb_idx ).mean(dim='storm')
            mean_surge[ii, bb] = z_ii_bin.values
    
    ds_out = xr.Dataset()
    ds_out['bin'] = pc_bins[:-1]
    ds_out['time'] = np.arange( n_time ) + 1
    ds_out['surge_norm'] = (['index','bin','time'], mean_surge)
    ds_out['pc_val'] = (['index','bin'], tsmax_pc.zeta.values[:-1,:].transpose())
    return ds_out

def analyse_return_periods( fp_zenv = './analysis/zenv_IBTRACS.nc', 
                            fp_grd = 'roms_grd.nc',
                            fp_out_rp = './analysis/return_periods.nc',
                            fp_out_coast = './analysis/coastal_return_periods.gpkg'):

    grd = xr.open_dataset(fp_grd)
    zenv = xr.open_dataset( fp_zenv )
    n_years = zenv.year.max().values + 1
    zenv = adjust_inundation( zenv, grd.h ).zenv

    ds_rp = return_period_emp_grid( zenv, n_time=n_years )
    rp_ocean, rp_land = separate_land_ocean( ds_rp.surge, grd.h )
    rp_land.name='flood'
    ds_rp = rp_ocean.to_dataset()
    ds_rp['flood'] = rp_land
    rp_coast= get_coastal_data( grd.h, ds_rp, kr=1 )
    gdf_coast = coastal_rp_to_gdf( rp_coast )

    gdf_coast.to_file(fp_out_coast)
    ds_rp.to_netcdf(fp_out_rp)


def align_timeseries_by_max( ds_ts, window=48 ):

    n_ts, n_pts = ds_ts.shape
    # Index window around maximum
    window_z = np.zeros( (2*window + 1, n_pts) ) * np.nan
    max_idx = np.argmax( ds_ts.values, axis=0 ).astype(int)
    for ii in range(n_pts):
        idx0 = max_idx[ii] - window
        idx1 = max_idx[ii] + window
        zeta = ds_ts.isel(index=ii).values

        if idx0 < 0:
            wind_idx0 = np.abs(idx0)
            window_z[wind_idx0:window+1, ii] = zeta[:max_idx[ii]+1]
        else:
            window_z[:window+1, ii] = zeta[idx0:max_idx[ii]+1]

        if idx1 >= len(zeta):
            wind_idx1 = len(zeta) - idx1
            window_z[window+1:wind_idx1-1, ii] = zeta[max_idx[ii]+1:]
        else:
            window_z[window+1:, ii] = zeta[ max_idx[ii]+1:idx1+1 ]
        
    ds_aligned = xr.Dataset()
    ds_aligned['zeta'] = (['time','index'], window_z)
    ds_aligned['lon'] = (['index'], ds_ts.lon_rho.values)
    ds_aligned['lat'] = (['index'], ds_ts.lat_rho.values)
    ds_aligned = ds_aligned.set_coords(['lon','lat'])
    return ds_aligned

def adjust_inundation( z, h, dcrit = 0.2):
    z = z.where( h>0, z + h - dcrit )
    return z

def separate_land_ocean( z, h ):
    # Separate into ocean and land arrays
    z_ocean = z.where(h>0)
    z_land = z.where(h<=0)
    return z_ocean, z_land

def coastal_rp_to_gdf( rp_coast ):
    n_rp = len(rp_coast.rp)
    lon = rp_coast.lon_rho.values
    lat = rp_coast.lat_rho.values
    points = [Point( lon[ii], lat[ii] ) for ii in range(len(lon))]
    gdf = gpd.GeoDataFrame( geometry = points, crs=4326 )
    for ii, rp in enumerate( rp_coast.rp.values ):
        rlii = rp_coast.surge[ii].values
        rlii[ rlii == 0 ] = np.nan
        gdf[f'surge_{rp}yr'] = rlii
    gdf = gdf.dropna()
    return gdf

def return_period_emp_grid( bmax_grid, n_time, rp=[50, 100, 200],):
    ds_out = bmax_grid.to_dataset()[['lon_rho','lat_rho']]
    bmax_grid = bmax_grid.values
    bshp = bmax_grid.shape
    n_pts = np.multiply(*bshp[1:])
    n_rp = len(rp)
    bmax_grid = bmax_grid.reshape( (bshp[0], n_pts ))
    rl = np.zeros( (n_rp,n_pts) ) * np.nan
    for ii in range(n_pts):
        bii = bmax_grid[:,ii]
        bii = bii[ bii > 1e-5 ]
        if len(bii) > 0:
            rl[:, ii] = return_period_emp( bii, n_time, rp)
    ds_out['rp'] = rp
    ds_out['surge'] = (['rp','eta_rho','xi_rho'], rl.reshape((n_rp, *bshp[1:])))
    return ds_out

def return_period_emp( bmax, n_time, rp = [50, 100, 200],
                       interp_rp = True):

    if np.sum( ~np.isnan(bmax) ) == 0:
        return np.zeros_like(rp)*np.nan

    # Remove any NaN values
    bmax = bmax[ ~np.isnan(bmax) ]
    
    bmax = np.sort(bmax)[::-1]
    rk = len(bmax) - stats.rankdata(bmax, method='min') + 1
    emp_rp = n_time/rk

    if interp_rp:
        interpf = interp1d( emp_rp, bmax, bounds_error=False )
        return_levels = interpf( rp )
        return return_levels
    else:
        return emp_rp, bmax
    
def get_coastal_data( h, data, kr = 1, return_indices = False ):
    ''' Identify indices of ocean points neighbouring land points '''
    hmask = h.values <= 0
    n_r, n_c = hmask.shape
    hmask = fill_isolated_regions(hmask).astype(int)
    
    r_ind = []
    c_ind = []
    for rr in range(kr, n_r-1):
        for cc in range(kr, n_c-1):
            ctr = hmask[rr, cc]
            kernel = hmask[ rr-kr: rr+kr+1, cc-kr:cc+kr+1 ]
            #kernel[0,0] = 0
            #kernel[-1,-1] = 0
            #kernel[-1,0] = 0
            #kernel[0,-1] = 0
            if np.sum(kernel) > 0 and ctr == 0:
                r_ind.append(rr)
                c_ind.append(cc)

    data_coast = data.isel( eta_rho = xr.DataArray( r_ind, dims=['index'] ), 
                            xi_rho = xr.DataArray( c_ind, dims=['index'] ) )

    if return_indices:
        return data_coast, r_ind, c_ind
    else:
        return data_coast
            

def get_average_grid_resolution( lon, lat, utm=False ):

    n_r, n_c = lon.shape

    x_res = np.zeros((n_r, n_c))*np.nan
    y_res = np.zeros((n_r,n_c))*np.nan

    for rr in range( 1, n_r-1 ):
        for cc in range( 1, n_c-1 ):

            if utm:
                x_diff = np.sqrt( (lon[rr,cc-1] - lon[rr,cc+1])**2 +
                                  (lat[rr,cc-1] - lat[rr,cc+1])**2 ) / 1000
                y_diff = np.sqrt( (lon[rr-1,cc] - lon[rr+1,cc])**2 +
                                  (lat[rr-1,cc] - lat[rr+1,cc])**2 ) / 1000
            else:
                x_diff = haversine( lon[rr, cc-1], lat[rr, cc-1], 
                                                   lon[rr, cc+1], lat[rr, cc+1],
                                                   radians=False )
                y_diff = haversine( lon[rr-1, cc], lat[rr-1, cc], 
                                                   lon[rr+1, cc], lat[rr+1, cc],
                                                   radians=False )
            x_res[rr, cc] = x_diff / 2
            y_res[rr, cc] = y_diff / 2

    return x_res, y_res

def get_nearest_timeseries_from_file( pt_lon, pt_lat, 
                                      fp_grd = 'roms_grd.nc', 
                                      fp_his = 'roms_his.nc',
                                      window = None):

    ds_grd = xr.open_dataset(fp_grd)
    ds_his = xr.open_dataset(fp_his)
    ds_ts = get_timeseries_nearest_point( ds_grd, ds_his, pt_lon, pt_lat, window=window )
    return ds_ts

def get_timeseries_nearest_point( ds_grd, ds, 
                                  pt_lon, pt_lat,
                                  window = None):

    lon2 = ds_grd.lon_rho.values
    lat2 = ds_grd.lat_rho.values
    mask = (ds_grd.h.values > 0)
    mask = fill_isolated_regions(mask)

    if type(pt_lon) is not list:
        pt_lon = [pt_lon]
    if type(pt_lat) is not list:
        pt_lat = [pt_lat]

    r_ind, c_ind = nearest_ocean_points( lon2, lat2, mask, pt_lon, pt_lat )

    ds_ts = ds.isel( eta_rho = xr.DataArray(r_ind), 
                     xi_rho = xr.DataArray(c_ind) )
    ds_ts = ds_ts.rename({'dim_0':'loc'})[['zeta','h']]

    # Index window around maximum
    if window is not None:
        window_z = np.zeros( 2*window + 1 ) * np.nan
        max_idx = np.argmax( ds_ts.zeta.values )
        idx0 = max_idx - window
        idx1 = max_idx + window
        zeta = ds_ts.zeta.values.squeeze()

        if idx0 < 0:
            wind_idx0 = np.abs(idx0)
            window_z[wind_idx0:window+1] = zeta[:max_idx+1]
        else:
            window_z[:window+1] = zeta[idx0:max_idx+1]

        if idx1 >= len(zeta):
            wind_idx1 = len(zeta) - idx1
            window_z[window+1:wind_idx1-1] = zeta[max_idx+1:]
        else:
            window_z[window+1:] = zeta[ max_idx+1:idx1+1 ]
        
        ds_ts2 = xr.Dataset()
        ds_ts2['zeta'] = (['time'], window_z)
        ds_ts2['h'] = ds_ts.h
        ds_ts = ds_ts2
    
    return ds_ts

def fill_isolated_regions( mask, area = 20 ):
    
    # Fill isolated areas in mask
    mask_labels = ndi.label(mask)
    for ii in range(mask_labels[1]):
        if np.sum(mask_labels[0] == ii) <= area:
            mask[ np.where(mask_labels[0] == ii) ] = 0
        else:
            pass

    return mask

def nearest_ocean_points( lon2, lat2, mask, pt_lon, pt_lat):
    c_index = np.ones_like(lon2)
    r_index = np.ones_like(lon2)
    for ii in range( lon2.shape[0] ):
        r_index[ii] = ii
    for ii in range( lon2.shape[1] ):
        c_index[:,ii] = ii 

    c_index = c_index[ mask == 1 ]
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
                             pt_lon_rad[ii], pt_lat_rad[ii])
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

    zenv = zenv.set_coords(['lon_rho','lat_rho'])

    if save_to_file:
        if os.path.exists( fp_out ):
            os.remove(fp_out)
        zenv.to_netcdf(fp_out)
    else:
        return zenv

def make_zero_surge_envelope( fp_grd = './roms_grd.nc', ds_grd = None, 
                              fp_out = 'roms_zenv.nc', save_to_file = False,
                              mask_land = True,
                              timeseries_window = None):

    if ds_grd is None:
        ds_grd = xr.open_dataset(fp_grd)
    
    zenv = ds_grd[['lon_rho','lat_rho']]
    zenv['zenv'] = (['eta_rho','xi_rho'], np.zeros_like( zenv.lon_rho ) )
    zenv.attrs = {}

    if mask_land:
        zenv = zenv.where( ds_grd.h > 0 )

    zenv = zenv.set_coords(['lon_rho','lat_rho'])

    if timeseries_window is not None:
        zenv['timeseries'] = (['time'], np.zeros( 2*window + 1))

    if save_to_file:
        if os.path.exists( fp_out ):
            os.remove(fp_out)
        zenv.to_netcdf(fp_out)
    else:
        return zenv
