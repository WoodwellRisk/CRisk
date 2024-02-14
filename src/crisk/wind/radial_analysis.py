import numpy as np
import xarray as xr
from crisk import utils
import pandas as pd
import dask
import dask.delayed as delayed
from dask.distributed import Client
from glob import glob
import os

def passing_storm_stats( track_list, grid_lon, grid_lat, year,
                         margin = 3,
                         radius=100,
                         vmax_bins = [18, 33, 43, 50, 58, 70],
                         n_years = None,
                         outdir = '.',
                         client = None,
                         min_events=2):
    '''

    '''
    
    n_tracks = len(track_list)
    n_bins = len(vmax_bins)

    @dask.delayed
    def par_loop(sidx, tracks, grid_lon, grid_lat, year, margin, radius, vmax_bins):
        analysis_loop = []
        for ii, tr in enumerate(tracks):
            tr = tr[['time','vmax','lat','lon']].set_index('time')
            tr = tr.resample('1H').interpolate()
            analysis_loop.append( passing_storm_grid( tr, 
                                        grid_lon, grid_lat, 
                                        margin, radius = radius,
                                        vmax_bins = vmax_bins).data_bin )
        analysis_loop = xr.concat( analysis_loop, dim='storm')
        analysis_loop['year'] = (['storm'], year)
        iistr = str(ii).zfill(6)
        analysis_loop.to_netcdf(f'{outdir}/tmp_ra_{sidx}.nc')
        del analysis_loop
        
    nthreads = [client.cluster.workers[ii].nthreads for ii in client.cluster.workers.keys() ]
    nthreads = np.sum(nthreads)
    batch_size = int( np.min( [n_tracks / nthreads + 1 , 500] ))
    to_compute = []
    for ii in range(0,n_tracks,batch_size):
        tracks_ii = track_list[ii:ii+batch_size]
        year_ii = year[ii:ii+batch_size]
        to_compute.append( par_loop( ii, tracks_ii, grid_lon, grid_lat, year_ii, 
                                     margin, radius, vmax_bins ) )
    computed = dask.compute(to_compute)[0]
    tmp_files = glob(f'{outdir}/*tmp*.nc')
    ds_tmp = [xr.open_dataset(fp, chunks={}) for fp in tmp_files]
    all_events = xr.concat(ds_tmp, dim='storm')
    annual_max = all_events.groupby('year').max()

    annual_count = []
    event_count = []
    for ii in range(n_bins):
        annual_count.append( (annual_max >= ii).sum(dim='year') )
        event_count.append( (all_events >= ii).sum(dim='storm') )
        
    annual_count = xr.concat(annual_count, dim='bin')
    annual_count['bin'] = vmax_bins
    event_count = xr.concat(event_count, dim='bin')
    event_count['bin'] = vmax_bins
    annual_count = annual_count.where( event_count >= min_events, 0)
    event_count = event_count.where( event_count >= min_events, 0)

    annual_count = annual_count.data_bin.to_dataset(name='annual_count')
    event_count = event_count.data_bin.to_dataset(name='event_count')
    ds_out = xr.merge([annual_count, event_count]).compute()

    # Delete tmp files
    for fp in tmp_files:
        os.remove(fp)
    return ds_out

def calculate_grid_radial_trackvar( df_track, 
                                    grid_lon, grid_lat, 
                                    radius=100, trackvarname = 'vmax'):
    ''' Get distances between all points in track dataframe and a grid '''

    # Convert all to radians
    grid_lon = np.radians(grid_lon)
    grid_lat = np.radians(grid_lat)
    track_lon = np.radians(df_track.lon.values)
    track_lat = np.radians(df_track.lat.values)
    trackvar = df_track[trackvarname].values

    n_pts = len(track_lon)
    n_grid = len(grid_lon)
    grid_out = np.zeros_like(grid_lon)*np.nan
    grid_count = np.zeros_like(grid_lon)

    for ii in range(n_pts):

        dist = utils.haversine( grid_lon, grid_lat, 
                                    track_lon[ii], 
                                    track_lat[ii], radians=True )
        dist_var = (dist < radius)
        grid_count[dist_var] = grid_count[dist_var] + 1
        dist_var = dist_var.astype(float)
        dist_var[ dist_var == 0 ] = np.nan
        dist_var = dist_var * trackvar[ii]
        grid_out = np.fmax( grid_out, dist_var )

    return grid_out, grid_count

def passing_storm_grid(df_track, grid_lon, grid_lat,
                        margin = 3,  
                        vmax_bins = [18, 33, 43, 50, 58, 70],
                        radius=100):
    '''
    
    '''

    # Number of intensity bins
    n_bins = len(vmax_bins)

    # Make full grid structure
    grid_poly = utils.get_grid_poly( grid_lon, grid_lat )

    # For this track, get distance to grid
    prdist = utils.distance_track_to_poly( df_track, grid_poly )
    df_track = df_track[['lon','lat','vmax']]

    # If track is too far away, don't analyse. Just return zeros
    if prdist > margin:
        grid_inten = np.zeros_like(grid_lon)
        grid_count = np.zeros_like(grid_lon)
    else:

        df_track = utils.clip_track_to_poly( df_track, grid_poly, 
                                             max_dist = margin, 
                                             round_days=False )
        grid_inten, grid_count = calculate_grid_radial_trackvar( df_track,
                                                     grid_lon, grid_lat,
                                                     radius=radius )
    bin_inten = np.zeros_like(grid_inten).astype(int) - 1
    for ii in range(n_bins):
        if ii < n_bins - 1:
            bin_idx = np.logical_and( grid_inten > vmax_bins[ii],
                                      grid_inten <= vmax_bins[ii+1] )
        else:
            bin_idx = grid_inten > vmax_bins[-1]
        bin_inten[bin_idx] = ii

    ds_out = xr.Dataset()
    ds_out['lon'] = (['y','x'], grid_lon)
    ds_out['lat'] = (['y','x'], grid_lat)
    ds_out['bin'] = (['bin'], vmax_bins)
    ds_out['data_bin'] = (['y','x'], bin_inten)
    #ds_out['data_bin'] = ds_out['data_bin'].astype('int32')
    #ds_out['data_count'] = (['y','x'], grid_count)
    return ds_out