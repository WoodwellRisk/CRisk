import xarray as xr
from climada.hazard import TCTracks
from climada.hazard import Centroids, TropCyclone
import numpy as np
import datetime

def tracks_to_wspd_rp(tracks, cent, pool, 
                      n_years = 10000, rperiods=[10,50,100,200,500],
                      reshape_2d = True):
    ''' Wrapper function to convert climada tracks to return periods. 
        The tracks input argument is a climada TCTracks() instance.
        cent is a climada Centroids instance. pool is a pathos 
        Pool instance, already initiated.
        The function returns an xarray dataset containing return levels.'''

    # Convert tracks into TropCyclone object at centroids
    tropcyc = TropCyclone.from_tracks(tracks, centroids=cent, 
                                      pool=pool, max_dist_inland_km=50,
                                      max_dist_eye_km = 500)
    tropcyc.set_frequency([1,n_years])
    
    # Calculate return levels for specified periods
    n_rp = len(rperiods)
    return_level = tropcyc.local_exceedance_inten(rperiods)

    # Mask 0s
    return_level[return_level == 0] = np.nan

    if reshape_2d:
        # Reshape return level array to be 2 dimensional
        cent_shape = cent.meta['height'], cent.meta['width']
        return_level = return_level.reshape((n_rp, cent_shape[0], cent_shape[1]))
        
        # Create output dataset and return
        ds_rp = xr.Dataset()
        lon2 = np.reshape(cent.lon, cent_shape)
        lat2 = np.reshape(cent.lat, cent_shape)
        ds_rp['lon'] = (['lon'], lon2[0,:] )
        ds_rp['lat'] = (['lat'], lat2[:,0] )
        ds_rp['return_period'] = (['return_period'], rperiods)
        ds_rp['return_level'] = (['return_period','lat','lon'], return_level)
    else:
        # Create output dataset and return
        ds_rp = xr.Dataset()
        ds_rp['lon'] = (['loc'], np.unique(cent.lon) )
        ds_rp['lat'] = (['loc'], np.unique(cent.lat) )
        ds_rp['return_period'] = (['return_period'], rperiods)
        ds_rp['return_level'] = (['return_period','loc'], return_level)

    # Set attributes
    ds_rp.attrs = {'title': 'Windspeed return levels calculated at specified return periods',
                   'n_years':n_years,
                   'history':f'Generated: {datetime.datetime.now()}'}
    ds_rp.return_period.attrs = {'units':'years'}
    ds_rp.return_level.attrs = {'units':'ms^-1'}
              
    # Delete the expanded wind fields to avoid memory issues
    del tropcyc

    return ds_rp
