import xarray as xr
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

def compare_windspeed_grid( data_past, data_future, diff_proportion = False ):
    ''' Creates a 3 panel plots showing past windspeeds, future windspeeds
    and the difference between them. Input datasets should be data arrays
    on the same grid and return period. '''

    # Get bounds for the figures
    minlon = data_past.lon.min().values
    maxlon = data_past.lon.max().values 
    minlat = data_past.lat.min().values
    maxlat = data_past.lat.max().values
    ctrlat = (maxlat + minlat) / 2
    ctrlon = (maxlon + minlon) / 2

    # Get color bounds
    cmin = np.nanmin( [np.nanmin(data_past.values), np.nanmin(data_future.values)] )
    cmax = np.nanmax( [np.nanmax(data_past.values), np.nanmax(data_future.values)] )

    # Get orientation
    width = maxlon - minlon
    height = maxlat - minlat

    # 1 = Portrait, 2 = Landscape
    if width > height:
        f,a = plt.subplots(3, 1, subplot_kw={'projection': ccrs.PlateCarree()},
                           sharey=True, sharex=True, figsize=((9,9)))
    else:
        f,a = plt.subplots(1, 3, subplot_kw={'projection': ccrs.PlateCarree()},
                           sharey=True, sharex=True, figsize = ((10,10)))

    for ii in range(3):
        gl = a[ii].gridlines(draw_labels=True, linewidth=1, 
                             color='gray', alpha=0.5, linestyle='--')
        # manipulate `gridliner` object
        gl.xlabels_top = False
        gl.ylabels_right = False
        if ii < 2:
            gl.xlabels_bottom = False
        #gl.xlines = False
        #gl.ylines = False
        a[ii].set_extent([minlon, maxlon, minlat, maxlat], 
                         crs=ccrs.PlateCarree())
        a[ii].coastlines(resolution='10m', color='black')

    # First panel is past
    a[0].set_title('1989 - 2015', fontsize=10)
    im = a[0].pcolormesh(data_past.lon, data_past.lat, data_past, 
                    cmap=plt.get_cmap('Greens',10), 
                    vmin=cmin, vmax = cmax)
    plt.colorbar(im, ax=a[0])

    # Second panel is future
    a[1].set_title('2015 - 2050', fontsize=10)
    im = a[1].pcolormesh(data_future.lon, data_future.lat, 
                    data_future, cmap=plt.get_cmap('Greens',10), 
                    vmin=cmin, vmax = cmax)
    plt.colorbar(im, ax=a[1])

    # Third panel is different
    a[2].set_title('Difference (future - past)', fontsize=10)
    if diff_proportion:
        diff = ((data_future / data_past) - 1)*100
        cminmax = np.nanmax( [np.abs( np.nanmin(diff) ), np.abs( np.nanmax(diff) ) ] )
        im = a[2].pcolormesh(diff.lon, diff.lat, diff, 
                             cmap=plt.get_cmap('RdBu_r',11), 
                             vmin=-10, vmax=10)
    else:
        diff = data_future - data_past
        cminmax = np.nanmax( [np.abs( np.nanmin(diff) ), np.abs( np.nanmax(diff) ) ] )
        im = a[2].pcolormesh(diff.lon, diff.lat, diff, 
                             cmap=plt.get_cmap('RdBu_r',11), 
                             vmin=-cminmax, vmax = cminmax)
    plt.colorbar(im, ax=a[2])

    f.tight_layout()

    rp = data_past.return_period.values
    f.text(0.5, 1.15, f'Windspeeds of the {rp}-year event (ms^-1)', fontweight='bold',
          transform = a[0].transAxes, ha='center', fontsize=12)

    return f, a

    