import numpy as np
from datetime import datetime, timedelta

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

def haversine_rad(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km