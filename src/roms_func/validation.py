import xarray as xr
from . import forcing, postprocessing
from glob import glob
import os 
import numpy as np
from shapely.geometry import Point, LineString
from utide import solve, reconstruct
from datetime import datetime, timedelta
import pandas as pd

def validate_from_file( fp_his = 'roms_his.nc', dir_tg = '../../data/uhslc' ):

    # Open input files
    ds_his = xr.open_dataset(fp_his)
    tg_list = glob( os.path.join(dir_tg, '*') ) 
    
    # Find tidegauges that lie within the model grid
    tg_list, fp_list = postprocessing.get_tidegauges_in_grid( ds_his, tg_list )

    # Figure out time period of model run
    roms_date0 = pd.to_datetime( ds_his.ocean_time[0].values )
    roms_date1 = pd.to_datetime( ds_his.ocean_time[-1].values )
    
    # Time period to analyse tides (2 years either side of event )
    ha_date0 = roms_date0 - timedelta(days=1*365)
    ha_date1 = roms_date1 + timedelta(days=1*365)
    
    # Open datasets and extract time periods
    for ii, tg in enumerate( tg_list ):
        tg_list[ii] = tg.sel(time=slice(ha_date0, ha_date1) )
    
    # Do a tidal analysis
    for ii, tg in enumerate( tg_list ):
        tg['sea_level'] = tg['sea_level'] / 1000
        ha = solve(tg.time.values, tg.sea_level[0].values, lat=tg.lat.values,
                   nodal=False, trend=False, method='ols',
                   conf_int='linear', Rayleigh_min=0.95)
        tide = reconstruct( tg.time, ha )
        tg['tide'] = (['time'], tide.h)
        tg['ntr'] = (['time'], tg.sea_level[0].values - tide.h )
        tg['std'] = np.nanstd(tg.ntr.values)

    # Get nearest 'ocean' points
    tg_lon = np.array([tg.lon.values for tg in tg_list])
    tg_lat = [tg.lat.values for tg in tg_list]
    tg_lon[tg_lon > 180] = tg_lon[tg_lon>180] - 360
    r_ind, c_ind = postprocessing.nearest_ocean_points( ds_his.lon_rho.values, 
                                         ds_his.lat_rho.values, 
                                         (ds_his.h>=0).values,
                                         tg_lon, tg_lat )

    # Extract model points
    ds_his_tg = ds_his.isel( eta_rho = xr.DataArray(r_ind), 
                             xi_rho = xr.DataArray(c_ind) ).compute()
    ds_his_tg = ds_his_tg[['zeta','h']].rename({'dim_0':'loc'})

    # Match timings
    ds_tg = [tg.interp( time = ds_his.ocean_time.values ) for tg in tg_list] 
    ds_tg = xr.concat(ds_tg, dim='loc')

    # Put out some metrics
    # Maxima
    model_max = ds_his_tg.zeta.max(dim='ocean_time').values
    obs_max = ds_tg.ntr.max(dim='time').values
    diff_max = model_max - obs_max
    
    # 90th percentile
    model_q = ds_his_tg.zeta.quantile(.95, dim='ocean_time').values
    obs_q = ds_tg.ntr.quantile(.95, dim='time').values
    diff_q = model_q - obs_q

    ds_out = xr.Dataset()
    ds_out['tg_lon'] = (['loc'], tg_lon[0])
    ds_out['tg_lat'] = (['loc'], tg_lat[0])
    ds_out['model_lon'] = (['loc'], ds_his_tg.lon_rho.values)
    ds_out['model_lat'] = (['loc'], ds_his_tg.lat_rho.values)
    ds_out['max_model'] = (['loc'], model_max)
    ds_out['max_obs'] = (['loc'],obs_max)
    ds_out['max_diff'] = (['loc'],diff_max)
    ds_out['q95_model'] = (['loc'], model_q)
    ds_out['q95_obs'] = (['loc'],obs_q)
    ds_out['q_diff'] = (['loc'],diff_q)

    return ds_out