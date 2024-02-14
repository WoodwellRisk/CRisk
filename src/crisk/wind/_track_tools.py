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
from shapely.geometry import Point, Polygon, LineString
import xarray as xr
from datetime import datetime, timedelta
from . import _utils
from climada.hazard import TCTracks

def filter_bad_rmw( track_list ):

    output_list = []

    for ii, tr in enumerate(track_list):
        zero_sum = np.sum( tr.radius_max_wind.values != 0 )
        nan_sum = np.sum( ~np.isnan( tr.radius_max_wind.values ) )
        if zero_sum > 0 and nan_sum > 0:
            output_list.append( tr )
    return output_list
    

def shift_tracks_lon( track_list ):
    for ii, track in enumerate(track_list):
        lon = track['lon'].values
        lon[lon>180] = lon[lon>180] - 360
        track['lon'] = (['time'], lon)
    return track_list

def filter_tracks_by_intensity( track_list, min_intensity = 64 ):
    ''' Filter tracks that never reach a minimum intensity '''
    keep_idx = []
    for ii, tr in enumerate(track_list):
        winds_over = tr.max_sustained_wind.values > min_intensity
        if np.sum( winds_over ) > 0:
            keep_idx.append(ii)
    return [track_list[ii] for ii in keep_idx]

def read_one_from_ibtracs( year=None, name=None, basin=None, sid=None):
    if sid is None:
        track = TCTracks.from_ibtracs_netcdf(year_range=[year,year], basin=basin)
        track_info = get_track_info( track.data )
        name_compare = [ _utils.compare_str(name, nii) for nii in track_info.name ]
        storm_idx = np.where( name_compare )[0][0]
        track.data = [track.data[storm_idx]]
    else:
        track = TCTracks.from_ibtracs_netcdf(storm_id=sid)
    return track

def get_track_info( track_list ):
    ''' Get dataframe of basic tc characteristics from list of climada tracks '''

    names = []
    sid = []
    year = []
    category = []
    
    for tt, track in enumerate(track_list):
        names.append(track.name)
        sid.append(track.sid)
        year.append( pd.to_datetime(track.time[0].values).year)
        category.append( track.category )

    df = pd.DataFrame()
    df['name'] = names
    df['sid'] = sid
    df['year'] = year
    df['category'] = category
    return df

def distance_track_to_poly( track_list, pol ):
    ''' Uses Shapely to check minimum proximity of a storm track to a box.
        The box is defined by specifying lonmin, lonmax, latmin and latmax'''
    
    linestrings = tracks_to_linestring( track_list )
    return pol.distance(linestrings)

def pad_track_start( track, start_time ):
    return

def clip_track_to_poly( track, poly, max_dist = 1, round_days=True ):

    n_tracks = len(track.data)
    track_clipped = []
    for ii in range(n_tracks):
        print(ii, n_tracks, end='\r')
        trackii = track.data[ii]
        t_points = list(zip( trackii.lon, trackii.lat) )
        points = [Point(tc) for tc in t_points]
        dist = np.array( [poly.distance(pt) for pt in points] )
        keep_idx = np.where(dist <= max_dist)[0]
        trackii_clipped = trackii.isel(time=slice( np.min(keep_idx), np.max(keep_idx ) ) )

        if round_days: 
            date0 = datetime(*pd.to_datetime(trackii_clipped.time.values[0]).timetuple()[:3])
            date1 = datetime(*pd.to_datetime(trackii_clipped.time.values[-1]).timetuple()[:3])
            date1 = date1 + timedelta(days=1)
            track_clipped.append( trackii.sel(time=slice(date0, date1) ) )
        else:
            track_clipped.append(trackii_clipped)

    track.data = track_clipped
    return track

def tracks_to_linestring( track_list ):
    ls_list =[ LineString(list(zip( tr.lon.values, tr.lat.values) ) )  for tr in track_list ]

    return ls_list

def subset_tracks_in_poly( track_list, pol, buffer = 0):
    ''' Subsets tracks into a geographical box. Tracks should be CLIMADA datasets
        in a list '''

    distances = distance_track_to_poly( track_list, pol )
    keep_idx = np.where( distances <= buffer )[0]
    new_tracks = [track_list[ii] for ii in keep_idx]

    return new_tracks

def subset_tracks_in_year( tracks, year ):
    ''' Subsets tracks into integer year '''
    track_list = tracks.data
    year_list = [ pd.to_datetime(tr.time[0].values).year for tr in track_list ]
    year_list = np.array(year_list)
    keep_idx = np.where(year_list == year)[0]
    new_tracks = TCTracks()
    new_tracks.data = [track_list[ii] for ii in keep_idx]
    return new_tracks

def get_storm_years( tracks ):
    return


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

    