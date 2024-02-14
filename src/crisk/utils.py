import numpy as np
from datetime import datetime, timedelta
from shapely.geometry import Point, Polygon, LineString
import pandas as pd

def subset_tracks_in_poly( track_list, pol, buffer = 0):
    ''' Identifies tracks (dataframe) in a list that pass through a specified
    polygon, with some added buffer. This is not an exact procedure and geographical distances
    are not used for the buffer. Instead, distance in degrees is used.
    
    Args:
        track_list (list): List of track pandas dataframes with lon, lat columns
        pol (shapely.geometry.Polygon): Polygon to compare
        buffer (float): Buffer around polygon to allow (in degrees)
    returns
        Indices of filtered tracks.
    '''
    if len(track_list) == 0:
        return []
    distances = distance_track_to_poly( track_list, pol )
    keep_idx = np.where( distances <= buffer )[0]
    return keep_idx

def track_to_linestring( track ):
    ''' Converts a track (or list of) dataframes into a list of shapely (lon,lat) LineStrings '''
    if type(track) is not list:
        track = [track]
        
    ls_list =[ LineString(list(zip( tr.lon.values, tr.lat.values) ) )  for tr in track ]

    if len(track) > 1:
        return ls_list
    else:
        return ls_list[0]

def get_grid_poly( lon, lat ):
    '''
    For a rectangular grid defined by 2D lon, lat arrays, return a 
    Shapely Polygon object representing its bounds.
    '''

    c1 = (lon[0,0], lat[0,0])      # Top left
    c2 = (lon[0,-1],lat[0,-1])     # Top right
    p1 = (lon[-1,0], lat[-1,0])    # Bottom left
    p2 = (lon[-1,-1], lat[-1,-1])  # Bottom right
    pol = Polygon( (c1, c2, p2, p1) )
    return pol

def distance_track_to_poly( track_list, pol ):
    ''' Uses Shapely to check minimum proximity of a storm track to a polygon.
        This does not calculate geographical distances, result will be in 
        degrees. This function is used for determining if a track passes through
        a polygon.'''
    linestrings = track_to_linestring( track_list )
    return pol.distance(linestrings)

def clip_track_to_poly( track, poly, max_dist = 1, round_days=True ):
    '''
    Takes a track dataframe and clips it to a shapely polygon.

    Args:
        track (pd.dataframe): Pandas dataframe track with lon, lat columns
        poly (shapely Polygon): Poly to clip to in same crs
        max_dist (float): Maximum distance around poly to clip (degrees)
        round_days (bool): If true, round resulting poly to day start
    '''
    t_points = list(zip( track.lon, track.lat) )
    points = [Point(tc) for tc in t_points]
    dist = np.array( [poly.distance(pt) for pt in points] )
    keep_idx = np.where(dist <= max_dist)[0]
    track_clipped = track.iloc[np.min(keep_idx):np.max(keep_idx)]

    if round_days: 
        date0 = datetime(*pd.to_datetime(track_clipped.time.values[0]).timetuple()[:3])
        date1 = datetime(*pd.to_datetime(track_clipped.time.values[-1]).timetuple()[:3])
        date1 = date1 + timedelta(days=1)
        didx = np.logical_and( track.time >= date0, track.time <= date1 )
        track_clipped = track.iloc[ np.where(didx)[0]]
    return track_clipped

def climada_to_dataframe( track, convert_units = True ):
    ''' Converts a climada track xarray dataset into an appropriate dataframe for ParaTC'''

    if 'lon' in track.dims:
        track['lon'] = (['time'], track.lon.values )
    
    df_track = track.rename({'radius_max_wind':'rmw', 
                             'central_pressure':'pcen',
                             'environmental_pressure':'penv',
                             'max_sustained_wind':'vmax'})
    df_track = df_track.to_dataframe().reset_index()
    if convert_units:
        df_track['rmw'] = df_track['rmw']*1.852
        df_track['vmax'] = (df_track['vmax']*0.51444) / 0.9
    return df_track

def compare_str( str1, str2 ):

    # Make both lowercase
    str1 = str1.lower()
    str2 = str2.lower()

    # Strip any whitespace
    str1 = str1.strip()
    str2 = str2.strip()

    return str1 == str2


def make_datetimes( year, month, day, hour ):
    ''' Convert a list of years, months, days and hours into a list of
    datetime objects '''
    dt = [datetime( year[ii], month[ii], 15 ) for ii in range(len(month))]
    for ii in range(len(hour)):
        dt[ii] = dt[ii] + timedelta(hours = hour[ii])
    return dt

def haversine(lon1, lat1, lon2, lat2, radians = True):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    if not radians:
        lon1 = np.radians(lon1)
        lat1 = np.radians(lat1)
        lon2 = np.radians(lon2)
        lat2 = np.radians(lat2)
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km