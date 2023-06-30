'''
This module contains routines that are specific to the STORM dataset.
URL: https://www.nature.com/articles/s41597-020-0381-2
'''

import pandas as pd

def read_STORM(fp, start_year = 0):
    '''Read data from a STORM formatted csv file into a Pandas dataframe.
       If start_year is specified, then all year values will be offsetted by
       this value.
    '''

    # Define header for STORM file and read file data
    header = header = ['year','month','num','time_step','basin_id','latitude',
                       'longitude','p_min','w_max','r_max','category',
                       'landfall','distance_to_land' ]
    df = pd.read_csv(fp, names = header)

    # Define an indicator variable -- start of storm
    df['ind'] = (df['time_step'] ==0).astype(int)

    # Define hour of storm (same as time_elapsed or similar)
    df['hour'] = df['time_step']*3

    # Convert longitude to range 0 - 360 (matches TCRM)
    dflon = df['longitude'].values
    dflon[dflon > 180] = dflon[dflon>180] - 360
    df['longitude'] = dflon

    # Add start_year to year
    df['year'] = df['year'] + start_year
    return df