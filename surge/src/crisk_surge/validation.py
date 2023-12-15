import xarray as xr
from . import forcing, postprocessing
from glob import glob
import os 
import numpy as np
from shapely.geometry import Point, LineString
from utide import solve, reconstruct
from datetime import datetime, timedelta
import pandas as pd

def check_obs_from_file( fp_tg, fp_frc = 'roms_frc.nc'):

    # Open input files
    ds_frc = xr.open_dataset(fp_frc)
    tg = xr.open_dataset(fp_tg)
    time_hour = ds_frc.sms_time.resample(sms_time='1H').mean()
    
    # Figure out time period of model run
    roms_date0 = pd.to_datetime( time_hour[0].values )
    roms_date1 = pd.to_datetime( time_hour[-1].values )
    
    # Time period to analyse tides (2 years either side of event )
    ha_date0 = roms_date0 - timedelta(days=1*365)
    ha_date1 = roms_date1 + timedelta(days=1*365)

    # Open datasets and extract time periods
    tg = tg.sel(time=slice(ha_date0, ha_date1) )

    # Return empties if there is no data
    if len(tg.time) == 0:
        return 0, 0
    
    # Check there is enough data for tidal analysis
    numdata_ha = np.sum( ~np.isnan( tg.sea_level[0].values ) )

    # Match timings
    tg = tg.interp( time = time_hour, method='linear' ).compute()
    
    numdata_event = np.sum( ~np.isnan( tg.sea_level[0].values ) )
    numdata_frac = numdata_event / len(tg.sea_level[0].values)

    return numdata_ha, numdata_frac
                         

def validate_from_file( fp_his = 'roms_his.nc', 
                        fp_tg = None, sid=None,
                        search_hours = 24, ):

    # Open input files
    ds_his = xr.open_dataset(fp_his)
    
    # Find tidegauges that lie within the model grid
    tg = xr.open_dataset(fp_tg)
    
    # Figure out time period of model run
    roms_date0 = pd.to_datetime( ds_his.ocean_time[0].values )
    roms_date1 = pd.to_datetime( ds_his.ocean_time[-1].values )

    # Figure out frequency of model time
    roms_dt = pd.to_datetime( ds_his.ocean_time[1].values ) - pd.to_datetime( ds_his.ocean_time[0].values )
    roms_dt = roms_dt.total_seconds()/(60**2)
    
    # Time period to analyse tides (2 years either side of event )
    ha_date0 = roms_date0 - timedelta(days=2*365)
    ha_date1 = roms_date1 + timedelta(days=2*365)
    
    # Open datasets and extract time periods
    tg = tg.sel(time=slice(ha_date0, ha_date1) )

    # Check there is enough data for tidal analysis
    numdata_ha = np.sum( ~np.isnan( tg.sea_level[0].values ) )
    
    # Do a tidal analysis
    tg = analyse_tg( tg )

    # Get nearest 'ocean' points
    tg_lon = tg.lon.values
    tg_lat = tg.lat.values
    if tg_lon > 180:
        tg_lon = tg_lon - 360
    mask = (ds_his.h>=0).values
    mask = postprocessing.fill_isolated_regions(mask)
    r_ind, c_ind = postprocessing.nearest_ocean_points( ds_his.lon_rho.values, 
                                         ds_his.lat_rho.values, 
                                         mask,
                                         tg_lon, tg_lat )
    # Extract model points
    ds_his_tg = ds_his.isel( eta_rho = xr.DataArray(r_ind), 
                             xi_rho = xr.DataArray(c_ind) ).compute()
    ds_his_tg = ds_his_tg[['zeta','h']].rename({'dim_0':'loc'})

    # Match timings
    tg = tg.interp( time = ds_his.ocean_time.values, method='linear' ).compute()
    numdata_event = np.sum( ~np.isnan( tg.sea_level[0].values ) )
    numdata_frac = numdata_event / len(tg.sea_level[0].values) 

    # Put out some metrics - Maxima first
    zeta = ds_his_tg.zeta.squeeze().values
    ntr = tg.ntr.squeeze().values
    nsearch = int(search_hours / roms_dt)
    model_max = np.nanmax(zeta)
    model_argmax = np.nanargmax(zeta)
    search_idx0 = max(model_argmax - nsearch, 0)
    search_idx1 = model_argmax + nsearch
    obs_max = np.nanmax( ntr[ search_idx0 : search_idx1 ] )
    diff_max = model_max - obs_max
    
    # 90th percentile
    model_q = np.percentile( zeta[ search_idx0 : search_idx1 ], 90 )
    obs_q = np.percentile( ntr[ search_idx0 : search_idx1 ], 90 )
    diff_q = model_q - obs_q

    df_out = pd.DataFrame()
    df_out['tg_lon'] =  tg_lon
    df_out['tg_lat'] =  tg_lat
    df_out['obs_std'] = tg['std'].values
    df_out['model_lon'] =  ds_his_tg.lon_rho.values
    df_out['model_lat'] =  ds_his_tg.lat_rho.values
    df_out['max_model'] = model_max
    df_out['max_obs'] = obs_max
    df_out['max_diff'] = diff_max
    df_out['q95_model'] =  model_q
    df_out['q95_obs'] = obs_q
    df_out['q_diff'] = diff_q
    df_out['numdata_frac'] = numdata_frac
    df_out['numdata_ha'] = numdata_ha
    if sid is not None:
        df_out['sid'] = sid

    return df_out, tg, ds_his_tg

def analyse_tg( tg ):

    tg = tg.copy()
    tg['sea_level'] = tg['sea_level'] / 1000
    ha = solve(tg.time.values, tg.sea_level[0].values, lat=tg.lat.values,
               nodal=True, trend=False, method='ols',
               conf_int='linear', Rayleigh_min=0.95)
    tide = reconstruct( tg.time, ha )
    tg['tide'] = (['time'], tide.h)
    tg['ntr'] = (['time'], tg.sea_level[0].values - tide.h )
    tg['std'] = np.nanstd(tg.ntr.values)
    return tg