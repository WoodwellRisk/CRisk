import numpy as np
from climada.hazard.tc_tracks import TCTracks, estimate_rmw
from climada.hazard import Centroids, TropCyclone
from climada.hazard import trop_cyclone
import climada.util.coordinates as u_coord
from shapely.geometry import Polygon
import xarray as xr
import pandas as pd
from . import postprocessing
import xesmf as xe
from crisk_surge import forcing
import os
import xesmf as xe

def interp_to_forcing_time( fp_frc, fp_bry ):

    bry = xr.open_dataset( fp_bry )
    frc = xr.open_dataset( fp_frc )
    bry = bry.interp( zeta_time = frc.sms_time.values )
    bry = bry.interp( v2d_time = frc.sms_time.values )
    bry.zeta_time.encoding['units'] = 'days since 1900-01-01'
    bry.v2d_time.encoding['units'] = 'days since 1900-01-01'

    if os.path.exists( 'roms_bry.nc' ):
        os.remove( 'roms_bry.nc' )
    bry.to_netcdf('roms_bry.nc')


def make_boundaries_from_parent( fp_grd_c, fp_grd_f, fp_data, fp_out=None ):

    grd_c = xr.open_dataset(fp_grd_c)
    grd_f = xr.open_dataset(fp_grd_f)
    ds_c = xr.open_dataset(fp_data)

    # Get h on u and v
    hu_c, hv_c = interp_h_to_uv( grd_c )
    hu_f, hv_f = interp_h_to_uv( grd_f )
    ds_c['umask'] = hu_c > 0
    ds_c['vmask'] = hv_c > 0
    ds_c['mask'] = ds_c.h > 0
    ds_c = ds_c.drop(['lon','lat'])

    grd_f['umask'] = hu_f > 0
    grd_f['vmask'] = hv_f > 0
    grd_f['mask'] = grd_f.h > 0
    grd_f = grd_f.drop(['lon','lat'])

    # Interpolate coarse UV onto coarse rho
    ds_r = interp_uv_to_rho( ds_c )
    ds_r = ds_r.set_coords(['lon_rho','lat_rho'])
    
    # Rotate coarse UV to face N/E
    ds_r['ubar'], ds_r['vbar'] = forcing.rotate_winds( ds_r.ubar, ds_r.vbar, 
                                                       grd_c.angle.values )
    ds_r['zeta'] = ds_c.zeta#.where( grd_c.h > 0, 0 )
    #ds_r['ubar'] = ds_r.ubar.where( grd_c.h > 0, 0 )
    #ds_r['vbar'] = ds_r.vbar.where( grd_c.h > 0, 0 )
    ds_r['mask'] = grd_c.h > 0
    
    # Interpolate onto finer grid and rotate to grid
    ds_f = interp_data_to_grid( ds_r, grd_f )
    ds_f['ubar'], ds_f['vbar'] = forcing.rotate_winds( ds_f.ubar, ds_f.vbar, -grd_f.angle.values )
    #ds_f['ubar'] = ds_f.ubar.where( grd_f.h >= 0, 0 )
    #ds_f['vbar'] = ds_f.vbar.where( grd_f.h >= 0, 0 )
    ds_f['h'] = grd_f.h
    return ds_f
    ds_u, ds_v = interp_rho_to_uv( ds_f, grd_f )
    ds_f['ubar'] = ds_u['ubar']
    ds_f['vbar'] = ds_v['vbar']
    return ds_f

    # Rename ocean_time to bry_time
    ds_f = ds_f.rename({'ocean_time':'bry_time'})
    ds_f.bry_time.encoding['units'] = 'days since 1900-01-01'

    # Apply wetdry mask properly
    ds_f['zeta'] = ds_f.zeta.where( grd_f.h >= 0, grd_f.h -0.2 )
    
    # Extract boundary data
    ds_bry = extract_boundaries( ds_f )
    if os.path.exists( 'roms_bry.nc' ):
        os.remove( 'roms_bry.nc' )
    ds_bry.to_netcdf('roms_bry.nc')
    
    return ds_bry, ds_f

def interp_and_rotate( ds, rmask, umask, vmask, angle ):
    
    grd_r = ds[['lon_rho','lat_rho']].rename({'lon_rho':'lon', 'lat_rho':'lat'})
    grd_r['mask'] = rmask

    grd_u = ds[['lon_u','lat_u']].rename({'lon_u':'lon', 'lat_u':'lat'})
    grd_u['mask'] = umask

    grd_v = ds[['lon_v','lat_v']].rename({'lon_v':'lon', 'lat_v':'lat'})
    grd_v['mask'] = vmask

    # Regrid U and V to rho
    rg_u = xe.Regridder( grd_u, grd_r, "bilinear", extrap_method="nearest_s2d" )
    ds_u = rg_u( ds['ubar'] )
    rg_v = xe.Regridder( grd_v, grd_r, "bilinear", extrap_method="nearest_s2d" )
    ds_v = rg_v( ds['vbar'] )
    ds_u, ds_v = forcing.rotate_winds( ds_u, ds_v, angle)

    # Regrid back to RHO
    

def extract_boundaries( ds ):

    ds_bry = xr.Dataset()

    # Extract zeta
    ds_bry['zeta_north'] = ds.zeta.isel(eta_rho=-1).rename({'bry_time':'zeta_time'})
    ds_bry['zeta_south'] = ds.zeta.isel(eta_rho=0).rename({'bry_time':'zeta_time'})
    ds_bry['zeta_west'] = ds.zeta.isel(xi_rho=0).rename({'bry_time':'zeta_time'})
    ds_bry['zeta_east'] = ds.zeta.isel(xi_rho=-1).rename({'bry_time':'zeta_time'})
    
    # Extract U
    ds_bry['ubar_north'] = ds.ubar.isel(eta_u=-1).rename({'bry_time':'v2d_time'})
    ds_bry['ubar_south'] = ds.ubar.isel(eta_u=0).rename({'bry_time':'v2d_time'})
    ds_bry['ubar_west'] = ds.ubar.isel(xi_u=0).rename({'bry_time':'v2d_time'})
    ds_bry['ubar_east'] = ds.ubar.isel(xi_u=-1).rename({'bry_time':'v2d_time'})

    # Extract V
    ds_bry['vbar_north'] = ds.vbar.isel(eta_v=-1).rename({'bry_time':'v2d_time'})
    ds_bry['vbar_south'] = ds.vbar.isel(eta_v=0).rename({'bry_time':'v2d_time'})
    ds_bry['vbar_west'] = ds.vbar.isel(xi_v=0).rename({'bry_time':'v2d_time'})
    ds_bry['vbar_east'] = ds.vbar.isel(xi_v=-1).rename({'bry_time':'v2d_time'})

    return ds_bry.fillna(0)

def extract_h_boundaries( ds ):

    hu, hv = interp_h_to_uv( ds )
    ds_bry = xr.Dataset()
    
    # Extract zeta
    ds_bry['h_north'] = ds.h.isel(eta_rho=-1)
    ds_bry['h_south'] = ds.h.isel(eta_rho=0)
    ds_bry['h_west'] = ds.h.isel(xi_rho=0)
    ds_bry['h_east'] = ds.h.isel(xi_rho=-1)
    
    # Extract U
    ds_bry['hu_north'] = hu.isel(eta_u=-1)
    ds_bry['hu_south'] = hu.isel(eta_u=0)
    ds_bry['hu_west'] = hu.isel(xi_u=0)
    ds_bry['hu_east'] = hu.isel(xi_u=-1)

    # Extract V
    ds_bry['hv_north'] = hv.isel(eta_v=-1)
    ds_bry['hv_south'] = hv.isel(eta_v=0)
    ds_bry['hv_west'] = hv.isel(xi_v=0)
    ds_bry['hv_east'] = hv.isel(xi_v=-1)

    return ds_bry

def interp_h_to_uv( ds ):

    ds = ds.rename({'lon_rho':'lon', 'lat_rho':'lat'})

    # Regrid stresses onto U and V grids
    ds_u = ds[['lon_u','lat_u']].set_coords(['lon_u','lat_u'])
    ds_u = ds_u.rename({'lon_u':'lon', 'lat_u':'lat'})
    regridder = xe.Regridder(ds, ds_u, "bilinear")
    ds_u = regridder( ds.h )
    
    ds_v = ds[['lon_v','lat_v']].set_coords(['lon_v','lat_v'])
    ds_v = ds_v.rename({'lon_v':'lon', 'lat_v':'lat'})
    regridder = xe.Regridder(ds, ds_v, "bilinear")
    ds_v = regridder( ds.h )
    return ds_u, ds_v

def interp_rho_to_uv(ds):

    # U vectors
    ds = ds.rename({'lon_rho':'lon', 'lat_rho':'lat'})

    # Regrid currents onto U and V grids
    ds_u = grd[['lon_u','lat_u']].set_coords(['lon_u','lat_u'])
    ds_u = ds_u.rename({'lon_u':'lon', 'lat_u':'lat'})
    regridder = xe.Regridder(ds, ds_u, "bilinear", extrap_method="nearest_s2d")
    ds_u = regridder( ds[['ubar','h','mask']] )
    ds_u = ds_u.rename({'lon':'lon_u', 'lat':'lat_u'})
    
    ds_v = grd[['lon_v','lat_v']].set_coords(['lon_v','lat_v'])
    ds_v = ds_v.rename({'lon_v':'lon', 'lat_v':'lat'})
    regridder = xe.Regridder(ds, ds_v, "bilinear", extrap_method="nearest_s2d")
    ds_v = regridder( ds[['vbar','h','mask']] )
    ds_v = ds_v.rename({'lon':'lon_v', 'lat':'lat_v'})
    return ds_u, ds_v

def interp_uv_to_rho( ds ):

    ds_r = ds[['lon_rho','lat_rho','mask','zeta']].rename({'lon_rho':'lon','lat_rho':'lat'})

    # Regrid stresses onto U and V grids
    ds_u = ds[['lon_u','lat_u','umask','ubar']].rename({'lon_u':'lon','lat_u':'lat', 
                                                         'umask':'mask'})
    regridder = xe.Regridder( ds_u, ds_r, "bilinear",extrap_method="nearest_s2d" )
    ds_r['ubar'] = regridder( ds_u ).ubar
    
    ds_v = ds[['lon_v','lat_v','vmask','vbar']].rename({'lon_v':'lon','lat_v':'lat', 
                                                         'vmask':'mask'})
    regridder = xe.Regridder(ds_v, ds_r, "bilinear",extrap_method="nearest_s2d")
    ds_r['vbar'] = regridder( ds_v ).vbar

    return ds_r.rename({'lon':'lon_rho','lat':'lat_rho'})

def interp_data_to_grid( data, new_grid ):

    ''' Interpolate variables on rho grid to UV grid '''

    # Regrid zeta
    data = data.rename({'lon_rho':'lon','lat_rho':'lat'})
    new_grid = new_grid.rename({'lon_rho':'lon','lat_rho':'lat'})[['lon','lat','mask']]
    regridder = xe.Regridder(data, new_grid, method='bilinear', extrap_method="nearest_s2d")
    ds_new = regridder( data )
    ds_new['lon_rho'] = new_grid.lon
    ds_new['lat_rho'] = new_grid.lat
    ds_new = ds_new.set_coords(['lon_rho','lat_rho'])
    
    return ds_new

def interp_onto_finer_grid( ds_c, ds_f ):

    # RHO points
    grd_rho_c = ds_c[['lon_rho','lat_rho','mask']]
    grd_rho_c = grd_rho_c.rename({'lon_rho':'lon', 'lat_rho':'lat'})
    grd_rho_f = ds_f[['lon_rho','lat_rho','mask']].rename({'lon_rho':'lon', 'lat_rho':'lat'})
    regridder = xe.Regridder(grd_rho_c, grd_rho_f, 
                             method='bilinear', extrap_method="nearest_s2d")
    rho_f = regridder( ds_c[['zeta']] )

    # U points
    grd_u_c = ds_c[['lon_u','lat_u','umask']]
    grd_u_c = grd_u_c.rename({'lon_u':'lon', 'lat_u':'lat', 'umask':'mask'})
    grd_u_f = ds_f[['lon_u','lat_u','umask']].rename({'lon_u':'lon', 'lat_u':'lat', 'umask':'mask'})
    regridder = xe.Regridder(grd_u_c, grd_u_f, 
                             method='bilinear', extrap_method="nearest_s2d")
    u_f = regridder( ds_c[['ubar']] )

    # V points
    grd_v_c = ds_c[['lon_v','lat_v','vmask']]
    grd_v_c = grd_v_c.rename({'lon_v':'lon', 'lat_v':'lat', 'vmask':'mask'})
    grd_v_f = ds_f[['lon_v','lat_v','vmask']].rename({'lon_v':'lon', 'lat_v':'lat', 'vmask':'mask'})
    regridder = xe.Regridder(grd_v_c, grd_v_f, 
                             method='bilinear', extrap_method="nearest_s2d")
    v_f = regridder( ds_c[['vbar']] )

    coords_f = ds_f[['lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v']]
    ds_f = xr.merge( [rho_f, u_f, v_f, coords_f] )
    return ds_f
    

    
    
