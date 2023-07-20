import numpy as np
import track_analysis
import xarray as xr

def multiple_analysis( df_idx_list, margin, lonmin, 
                       lonmax, latmin, latmax, resolution, 
                       radius, delta):

    n_idx = len(df_idx_list)
    output_data = []

    for kk in range(n_idx):
        out_single = single_analysis( df_idx_list[kk], margin, 
                                      lonmin, lonmax, latmin, latmax, resolution,
                                      radius = radius, delta=delta)
        if out_single is None:
            continue
        else:
            ds = out_single
            output_data.append(ds)

    if len(output_data) == 0:
        return None
    
    ds_concat = xr.concat( output_data, dim='storm' )
    ds_out = xr.Dataset()
    for cat in [0, 1, 2, 3, 4, 5]:
        ds_tmp = ds_concat.data == cat
        ds_tmp = ds_tmp.sum( dim='storm' )
        ds_out[f'category_{cat}'] = ds_tmp

    return ds_out

def single_analysis(df_ii, margin, 
                    lonmin, lonmax, latmin, latmax, 
                    resolution, dir_tmp = './tmp', 
                    radius=200, delta = 1/3):
    '''

    '''

    grid_lon = np.arange(lonmin, lonmax, resolution)
    grid_lat = np.arange(latmin, latmax, resolution)
    lon2, lat2 = np.meshgrid(grid_lon, grid_lat)
    lonF = lon2.flatten()
    latF = lat2.flatten()
    
    prdist = track_analysis.track_distance_to_box( df_ii, lonmin, lonmax, latmin, latmax )
    n_r, n_c = lon2.shape
    df_ii = df_ii[['longitude','latitude','category']]

    if prdist > margin:
        t = np.zeros((6, n_r, n_c))
        return
    else:
        df_interp = track_analysis.interpolate_track( df_ii, delta=delta )
        df_interp['category'] = df_interp['category'].round(0).astype(int)
        
        t = track_analysis.track_distance_to_grid( df_interp,
                             lonF, latF,
                             radius=radius )
        t = t.reshape((6, n_r, n_c))
        
    for ii in range(6):
        t[ii] = t[ii]*(ii+1)
    t = np.max(t, axis = 0) - 1
    
    ds_out = xr.Dataset()
    ds_out['x'] = (['x'], grid_lon)
    ds_out['y'] = (['y'], grid_lat)
    ds_out['data'] = (['y','x'], t)
    return ds_out