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

    # Interpolate coarse UV onto coarse rho
    ds_r = interp_uv_to_rho( ds_c.ubar, ds_c.vbar, grd_c )
    ds_r = ds_r.set_coords(['lon_rho','lat_rho'])
    
    # Rotate coarse UV to face N/E
    ds_r['ubar'], ds_r['vbar'] = forcing.rotate_winds( ds_r.ubar, ds_r.vbar, 
                                                       grd_c.angle.values )
    ds_r['zeta'] = ds_c.zeta.where( grd_c.h > 0, 0 )
    ds_r['ubar'] = ds_r.ubar.where( grd_c.h > 0, 0 )
    ds_r['vbar'] = ds_r.vbar.where( grd_c.h > 0, 0 )
    
    # Interpolate onto finer grid and rotate to grid
    ds_f = interp_data_to_grid( ds_r, grd_f )
    ds_f['ubar'], ds_f['vbar'] = forcing.rotate_winds( ds_f.ubar, ds_f.vbar, -grd_f.angle.values )
    ds_f['ubar'] = ds_f.ubar.where( grd_f.h >= 0, 0 )
    ds_f['vbar'] = ds_f.vbar.where( grd_f.h >= 0, 0 )
    ds_f['ubar'], ds_f['vbar'] = forcing.interpolate_rho_to_uv( grd_f, ds_f, 'ubar', 'vbar')

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
    
    return ds_bry

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

def interp_uv_to_rho( ds_u, ds_v, ds_grd ):

    ds_r = ds_grd[['lon_rho','lat_rho']].rename({'lon_rho':'lon','lat_rho':'lat'})

    # Regrid stresses onto U and V grids
    ds_u = ds_u.rename({'lon_u':'lon','lat_u':'lat'})
    regridder = xe.Regridder(ds_u, ds_r, "bilinear")
    ds_r['ubar'] = regridder( ds_u )
    
    ds_v = ds_v.rename({'lon_v':'lon','lat_v':'lat'})
    regridder = xe.Regridder(ds_v, ds_r, "bilinear")
    ds_r['vbar'] = regridder( ds_v )

    return ds_r.rename({'lon':'lon_rho','lat':'lat_rho'})

def interp_data_to_grid( data, new_grid ):

    ''' Interpolate variables on rho grid to UV grid '''

    # Regrid zeta
    data = data.rename({'lon_rho':'lon','lat_rho':'lat'})
    ds_new = xr.Dataset()
    ds_new['lon'] = new_grid.lon_rho #(['eta_u','xi_u'], new_grid.lon_u.values)
    ds_new['lat'] = new_grid.lat_rho #(['eta_u','xi_u'], new_grid.lat_u.values)
    regridder = xe.Regridder(data, ds_new, method='bilinear')
    ds_new = regridder( data )
    ds_new['lon_rho'] = new_grid.lon_rho
    ds_new['lat_rho'] = new_grid.lat_rho
    ds_new = ds_new.set_coords(['lon_rho','lat_rho'])
    
    return ds_new
