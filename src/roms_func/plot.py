import xarray as xr
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

def plot_surge( data, track = None, elev_max = 10, surge_max = 2,
               lonbnds=None, latbnds=None):

    elev = data.h
    wdmask = data.wetdry_mask_rho

    # Create figure
    f,a = plt.subplots(1,1)

    # Adjust h so that over the ocean it is 0
    h_adj = xr.where( data.h<=0, 0, data.h)

    # Adjust z so that we mask out all land areas
    z = data.zeta
    z = xr.where(data.h<0, np.nan, z)

    # Plot the elevation where land
    a.pcolormesh( data.lon_rho[1:-1, 1:-1], 
                  data.lat_rho[1:-1, 1:-1], 
                  -h_adj[1:-1,1:-1], cmap=plt.get_cmap('YlOrBr'))

    # Plot the ocean height
    pc = a.pcolormesh( z.lon_rho[1:-1, 1:-1], z.lat_rho[1:-1, 1:-1], 
                       z[1:-1, 1:-1], vmin=-2, vmax=2, 
                       cmap='RdBu')
    
    if track is not None:
        a.scatter(track.track_lon, track.track_lat, c='k', marker='x')
        
    if lonbnds is None:
        a.set_xlim( np.min(data.lon_rho), np.max(data.lon_rho))
    else:
        a.set_xlim( lonbnds[0], lonbnds[1])

    if latbnds is None:
        a.set_ylim( np.min(data.lat_rho), np.max(data.lat_rho))
    else:
        a.set_ylim( latbnds[0], latbnds[1])
    #a.set_title( data.ocean_time.values )
    
    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    f.colorbar(pc, cax=cbar_ax, extend='both')
    
    return f

def plot_flood( data, track = None, elev_max = 10, 
                lonbnds = None, latbnds = None):

    elev = data.h
    wdmask = data.wetdry_mask_rho

    # Create figure
    f,a = plt.subplots(1,1)

    # Adjust h so that over the ocean it is 0
    h_adj = xr.where( data.h<0, 0, data.h)

    # Adjust z so that we mask out solid mask areas with NaN
    z = data.zeta
    z = xr.where(data.mask_rho == 0, np.nan, z)

    # Adjust z by adding h_adj = 0 over ocean
    z_adj = z[1:-1, 1:-1] + h_adj[1:-1, 1:-1]

    # Plot the ocean height
    pc = a.pcolormesh( z.lon_rho[1:-1, 1:-1], z.lat_rho[1:-1, 1:-1], 
                       z_adj, vmin=0, vmax=5, 
                       cmap=plt.get_cmap('Blues',10))

    # Plot the elevation where ocean is dry
    elevii = xr.where(wdmask==1, np.nan, elev)
    a.pcolormesh( data.lon_rho[1:-1, 1:-1], 
                  data.lat_rho[1:-1, 1:-1], 
                  -elevii[1:-1,1:-1], cmap=plt.get_cmap('YlOrBr'), 
                  vmin=-10, vmax=50)

    # a.contour(    data.lon_rho[1:-1, 1:-1], 
    #               data.lat_rho[1:-1, 1:-1], 
    #               data.h[1:-1, 1:-1], zorder=100, 
    #               levels=0, colors=[[1, 102/255, 1]], linewidths=1)
    
    if track is not None:
        a.scatter(track.track_lon, track.track_lat, c='k', marker='x')

    if lonbnds is None:
        a.set_xlim( np.min(data.lon_rho), np.max(data.lon_rho))
    else:
        a.set_xlim( lonbnds[0], lonbnds[1])

    if latbnds is None:
        a.set_ylim( np.min(data.lat_rho), np.max(data.lat_rho))
    else:
        a.set_ylim( latbnds[0], latbnds[1])
        
    a.set_title( data.ocean_time.values )

    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    f.colorbar(pc, cax=cbar_ax, extend='max')
    
    return f