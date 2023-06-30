'''
This module contains routines for manipulating a track dataframe.
Some assumptions are made as to the structure of this dataframe.

Expected variable names:

> year : year of storm
> ind  : Storm start indicator. This should = 1 for the first timestep of
         a track in a dataframe.
> longitude : Longitude of storm track at timestep
> latitude  : Latitude of storm track at timestep

'''

import pandas as pd
from datetime import datetime, timedelta
import numpy as np
from scipy.interpolate import interp1d
from shapely import Point, Polygon, LineString
import xarray as xr
from datetime import datetime, timedelta
from . import _utils

def separate_years( df ):
    ''' Separate a track dataframe (df) into years. This routine will
    return a new list of dataframes, each corresponding to the values
    of 'year' in the input dataset '''

    # Get years as integer array and define starting indices
    year = df.year.values.astype(int)
    start_indices = np.where(year[:-1] != year[1:])[0] + 1

    # Initialise output list
    df_list = []
    n_years = len(start_indices)

    # Loop over years and append dataframes to list
    for ii in range(n_years):
        if ii < n_years-1:
            df_ii = df[start_indices[ii]:start_indices[ii+1]] 
        else:
            df_ii = df[start_indices[-1]:]

        df_list.append(df_ii.reset_index(drop=True))
    return df_list

def separate_events( df ):
    ''' Separates discrete events from a track dataframe
        Returns a list of event dataframes. Events are identified
        by the 'ind' column, which are 1 or True for the first
        timestep of a storm. '''

    df_list = []
    start_indices = np.where(df.ind == 1)[0]
    n_events = len(start_indices)
    for ii in range(n_events):
        if ii < n_events-1:
            df_ii = df.loc[start_indices[ii]:start_indices[ii+1]-1] 
        else:
            df_ii = df.loc[start_indices[-1]:]

        df_list.append(df_ii.reset_index(drop=True))
    
    return df_list

def interpolate_track( df, delta = 1/3 ):
    ''' Interpolates storm track dataframe to a finer timestep. 
        New timestep is defined by delta, which should be a fraction.
        For example delta = 1/2 would half the timestep. '''
    if delta == 1:
        return df
    else:
        ds = df.to_xarray()
        new_index = np.arange(0, len(ds.index), delta)
        ds_interp = ds.interp(index = new_index)
        df_interp = ds_interp.dropna(dim='index').to_dataframe()
        return df_interp

def track_distance_to_box( df, lonmin, lonmax, latmin, latmax ):
    ''' Uses Shapely to check proximity of a storm track to a box.
        The box is defined by specifying lonmin, lonmax, latmin and latmax'''
    
    pol = Polygon( ([lonmin, latmin], [lonmin, latmax],
                    [lonmax, latmax],[lonmax, latmin]) )
    t_points = list(zip( df.longitude, df.latitude) )
    p_list = []
    for t in t_points:
        p_list.append( Point(t) )
    
    ls = LineString(p_list)
    
    return pol.distance(ls)
    
def track_distance_to_grid( df_track, lon1, lat1, radius=100 ):
    ''' Get distances between all points in track dataframe and a grid '''

    # Convert all to radians
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    track_lon = np.radians(df_track.longitude.values)
    track_lat = np.radians(df_track.latitude.values)
    category = df_track.category.values

    n_pts = len(track_lon)
    n_grid = len(lon1)
    category_grid = np.zeros((6, n_grid))

    for ii in range(n_pts):
        dist = dist = _utils.haversine_rad( lon1, lat1, 
                                            track_lon[ii], 
                                            track_lat[ii] )
        distb = dist < radius
        cat_ii = int(category[ii])
        category_grid[ cat_ii ] = category_grid[ cat_ii ] + distb

    return np.clip(category_grid,0,1)
    