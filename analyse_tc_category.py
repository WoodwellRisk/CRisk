import pandas as pd
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import dask
import dask.delayed as delayed
from dask.distributed import Client
from glob import glob
import track_analysis as tana
from track_analysis import radial_analysis
import multiprocessing
import rioxarray

def main():

    # Loop over models in sepcified list
    ds_list = []
    
    # List of models to analyse. Model data read as {dir_input}/{model_name}/*{basin}*
    model_list = ['IBTRACS','CMCC','HADGEM','ECEARTH','CNRM']
    
    for model_name in model_list:
        print(model_name)
    
        # Get list of inputs and read into dataframe
        fp_inputs = glob(f'{dir_input}/{model_name}/*{basin}*')
        tracks_list = [ tana.read_STORM( fp ) for fp in fp_inputs ]
        tracks = pd.concat(tracks_list).reset_index()
        print('    Tracks read.')
    
        # Separate events into years and randomly select n_years
        tracks_years = tana.separate_years( tracks )
        np.random.shuffle(tracks_years)
        tracks = pd.concat(tracks_years[:n_years]).reset_index()
    
        # Separate tracks into individual events
        tracks_events = tana.separate_events( tracks )
        n_events = len(tracks_events)
    
        print('    Tracks separated into events.')
    
        client = Client()
        
        pdel = delayed(radial_analysis.multiple_analysis)
        compute_list = []
        batch_size = 2000
        for idx0 in np.arange(0,n_events,batch_size):
            idx1 = idx0 + batch_size
            tracks_ii_list = tracks_events[idx0:idx1]
            compute_list.append( pdel( tracks_ii_list, margin, 
                                       lonmin, lonmax, latmin, latmax, 
                                       resolution, radius, delta ) )
        
        out = dask.compute(compute_list, scheduler = 'processes')[0]

        # Strip out processes where no storms passed
        out = [i for i in out if i is not None]
    
        print('    Parallel work done.')
        
        ds_concat = xr.concat(out, dim='storm')
        ds_total_count = ds_concat.sum(dim='storm').compute()
        ds_count_per_year = ds_total_count / n_years
    
        attrs = {'title': f'Average number of tropical cyclones of varying categories passing within {radius}km per year.',
                 'radius': f'{radius}km',
                 'years_analysed': f'{n_years}',
                 'track_dataset':'STORM',
                 'model':f'{model_name}'}
        ds_count_per_year.attrs = attrs
        ds_list.append(ds_count_per_year)
    
        # Clean up
        client.close()
        del tracks
        del tracks_list
        del tracks_events

    ds_ibtracs = ds_list[0]
    ds_cmip = xr.concat(ds_list[1:], dim='model').mean(dim='model')
    attrs = {'title': f'Average number of tropical cyclones of varying categories passing within {radius}km per year.',
                 'radius': f'{radius}km',
                 'years_analysed': f'{n_years}',
                 'track_dataset':'STORM',
                 'model':f'Model Mean'}
    ds_cmip.attrs = attrs
    ds_ibtracs = ds_ibtracs.rio.write_crs(4326)
    ds_cmip = ds_cmip.rio.write_crs(4326)
    
    categories = [0, 1, 2, 3, 4, 5]
    for cat in categories:
        ds_cmip[f'category_{cat}'].rio.to_raster( os.path.join( dir_output, 
                                   f'TC_annualprob_STORM_{basin}_{radius}km_cat{cat}_2015_2050.tif' ) )
        ds_ibtracs[f'category_{cat}'].rio.to_raster( os.path.join( dir_output, 
                                   f'TC_annualprob_STORM_{basin}_{radius}km_cat{cat}_1971_2000.tif' ) )
    
    # GREATER THAN DATASETS
    for ii in range(5):
        cc_list = categories[ii:]
        cc = categories[ii]
        ds_proj_ii = [ ds_cmip[f'category_{cc}'] for cc in cc_list]
        ds_past_ii = [ ds_ibtracs[f'category_{cc}'] for cc in cc_list]
    
        ds_proj_ii = 1/xr.concat(ds_proj_ii, dim='category').sum(dim='category')
        ds_proj_ii = ds_proj_ii.where( ds_proj_ii <= 200 ).rename(f'gt_category_{categories[ii]}')
        ds_past_ii = 1/xr.concat(ds_past_ii, dim='category').sum(dim='category')
        ds_past_ii = ds_past_ii.where( ds_past_ii <= 200 ).rename(f'gt_category_{categories[ii]}')
    
        ds_proj_ii.rio.to_raster( os.path.join( dir_output, f'TC_returnperiod_STORM_{basin}_{radius}km_gtcat{cc}_2015_2050.tif' ) )
        ds_past_ii.rio.to_raster( os.path.join( dir_output, f'TC_returnperiod_STORM_{basin}_{radius}km_gtcat{cc}_1971_2000.tif' ) )


if __name__ == "__main__":

    # DIRECTORIES
    dir_input = '../STORM'
    dir_output = './output'
    
    # STORM basin from which to read input files
    basin = 'NA' # NA, SP, WP, EI or WI
    output_name = basin
    
    # ANALYSIS GRID SETTINGS
    resolution = 1/5 # Grid resolution in degrees
    lonmin = -107 # Grid minimum longitude
    lonmax = 2  # Grid maximum longitude
    latmin = 3  # Grid minimum latitude
    latmax = 63  # Grid maxmimum latitude
    margin = 3  # Margin around grid within which to generate storms (degrees)
    
    # Example (lonmin, lonmax, latmin, latmax) for basins
    # NA: (-107, 2, 63, 3)
    # SP: ()
    # WP: (100, 180, 0, 80)
    # SI: (0, 140, -60, 0)
    # WI: ()
    
    # ANALYSIS SETTINGS
    radius = 100
    n_years = 2500
    delta = 1/3
    
    main()