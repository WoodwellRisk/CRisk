import numpy as np
from climada.hazard.tc_tracks import TCTracks, estimate_rmw
from climada.hazard import Centroids, TropCyclone
import climada.util.coordinates as u_coord
from shapely.geometry import Polygon
import xarray as xr
import pandas as pd
from . import postprocessing

def make_era5_forcing( ds_grd, ds_era5 ):
    return

def get_forcing_info( fp_frc = './roms_frc.nc' ):

    ds_frc = xr.open_dataset(fp_frc)
    min_time = pd.to_datetime( ds_frc.sms_time.min().values )
    max_time = pd.to_datetime( ds_frc.sms_time.max().values )
    n_seconds = (max_time - min_time).total_seconds()
    return min_time, max_time, n_seconds

def get_grid_poly( lon, lat ):

    c1 = (lon[0,0], lat[0,0])      # Top left
    c2 = (lon[0,-1],lat[0,-1])     # Top right
    p1 = (lon[-1,0], lat[-1,0])    # Bottom left
    p2 = (lon[-1,-1], lat[-1,-1])  # Bottom right
    pol = Polygon( (c1, c2, p2, p1) )
    return pol

def rotate_winds( x,y, angle ):
    ''' Rotate wind vectors counter clockwise '''

    cosang = np.cos(angle)
    sinang = np.sin(angle)
    new_x = x * cosang - y * sinang
    new_y = x * sinang + y * cosang

    return new_x, new_y

def make_uniform_windfield( lon, lat, time, u, v, angle=None ):

    wind_u = np.ones( (len(time), lon.shape[0], lon.shape[1]) )*u
    wind_v = np.ones( (len(time), lon.shape[0], lon.shape[1]) )*v

    if angle is not None:
        wind_u, wind_v = rotate_winds(wind_u, wind_v, -angle)
    
    return wind_u, wind_v

def make_pressure_field(lon, lat, ds_track, 
                         model = 'H1980'):

    # Make distance array
    n_nodes = len( ds_track.time )
    pressure = np.zeros((n_nodes, lon.shape[0], lon.shape[1]))
    rmw = estimate_rmw(ds_track.radius_max_wind.values,
                       ds_track.central_pressure.values )
                       
    for ii in range(n_nodes):
        distkm = postprocessing.haversine( lon, lat, 
                                           ds_track.lon.values[ii], 
                                           ds_track.lat.values[ii], 
                                           radians=False )
        B_param = _B_holland_1980( ds_track.max_sustained_wind.values[ii],
                                   ds_track.environmental_pressure.values[ii],
                                   ds_track.central_pressure.values[ii])
        
        pii = holland_pressure( distkm, rmw[ii],
                                     B_param, 
                                     ds_track.environmental_pressure.values[ii],
                                     ds_track.central_pressure.values[ii], 
                                     ds_track.lat.values[ii] )
            
        pressure[ii] = pii
    return pressure

def _B_holland_1980(
    wmax: np.ndarray,
    penv: np.ndarray,
    pcen: np.ndarray,
):
    """Holland's 1980 B-value computation for gradient-level winds.
    """
    GRADIENT_LEVEL_TO_SURFACE_WINDS = 0.9
    RHO_AIR = 1.15

    wmax = wmax / GRADIENT_LEVEL_TO_SURFACE_WINDS
    
    # the factor 100 is from conversion between mbar and pascal
    pdelta = 100 * (penv - pcen)
    hol_b = wmax**2 * np.exp(1) * RHO_AIR / np.fmax(np.spacing(1), pdelta)
    return np.clip(hol_b, 1, 2.5)

def holland_pressure( d_centr: np.ndarray,
                      r_max: np.ndarray,
                      hol_b: np.ndarray,
                      penv: np.ndarray,
                      pcen: np.ndarray,
                      lat: np.ndarray,
    ) -> np.ndarray:
    """Symmetric and static wind fields (in m/s) according to Holland 1980.

     surface winds by adjusting the parameter `hol_b` (see function `_bs_holland_2008`).
    """
    pres = np.zeros_like(d_centr)
    too_close = d_centr < 0.1
    
    r_max_norm = (r_max / d_centr)**hol_b
    pres = pcen + (penv - pcen)*np.exp( -r_max_norm )
    pres[too_close] = pcen
    return pres

def _coriolis_parameter(lat: np.ndarray) -> np.ndarray:
    """Compute the Coriolis parameter from latitude."""
    V_ANG_EARTH = 7.29e-5
    return 2 * V_ANG_EARTH * np.sin(np.radians(np.abs(lat)))

def make_windfield( lon, lat, track, max_dist_eye_km = 400, angle=None,
                   model = 'H1980'):
    ''' Make windfield xarray dataset from a climada track. Wind vectors are relative to N/E '''

    # Make centroids from grid file
    cent = Centroids( lat=lat.flatten(), lon=lon.flatten() )
    cent.check()
    
    # Make u and v windfields
    wind = TropCyclone.from_tracks(track, cent, store_windfields=True,
                                   max_dist_eye_km = max_dist_eye_km, 
                                   model=model)
    
    # Reshape wind output into time x lat x lon
    wind2 = wind.windfields[0].toarray()
    n_time = wind2.shape[0]
    
    wind_v = np.array([wind2[ii,::2].reshape(lon.shape) for ii in range(n_time)])
    wind_u = np.array([wind2[ii,1::2].reshape(lon.shape) for ii in range(n_time)])

    if angle is not None:
        wind_u, wind_v = rotate_winds(wind_u, wind_v, -angle)

    return wind_u, wind_v
    
def tau_large_pond( U, V, rho = 1.290 ):
    ''' Wind stress vectors according to Large & Pond '''
    s = np.sqrt(U**2 + V**2)
    Cd = 0.001*(0.49+0.065*s)
    tau = s * rho * Cd
    tau_u = U * s * rho * Cd
    tau_v = V * s * rho * Cd
    return tau, tau_u, tau_v

def tau_andreas( U, V, rho = 1.290 ):
    ''' Wind stress vectors according to Andreas '''
    s = np.sqrt(U**2 + V**2)
    ws1 = s - 9.271
    ws2 = np.sqrt( 0.12*ws1**2 + 0.181 )
    ws = .239 + .0433*( ws1 + ws2 )
    Cd = (ws / s)**2
    Cd[Cd > 1] = 0
    tau = s * rho * Cd
    tau_u = U * s * rho * Cd
    tau_v = V * s * rho * Cd
    return tau, tau_u, tau_v

def make_forcing_dataset( ds_grd, time, 
                          sustr = None, svstr=None, 
                          press = None, 
                          track_lon = None, track_lat = None ):

    ds_tmp = ds_grd[['lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v']]
    ds_tmp['sms_time'] = (['sms_time'], time)
    ds_tmp['pair_time'] = (['pair_time'], time)

    if sustr is not None:
        ds_tmp['sustr'] = (['sms_time','eta_u', 'xi_u'], sustr)
        ds_tmp['sustr'].attrs = {'long_name':"surface u-momentum stress",
                                 'units':"m/s",
                                 'field':"surface u-momentum stress",
                                 'time':"sms_time",
                                 'coordinates':'lon_u lat_u' }
    if svstr is not None:
        ds_tmp['svstr'] = (['sms_time','eta_v', 'xi_v'], svstr)
        ds_tmp['svstr'].attrs = {'long_name':"surface v-momentum stress",
                                 'units':"Newton meter-2",
                                 'field':"surface v-momentum stress",
                                 'time':"sms_time",
                                 'coordinates':'lon_v lat_v' }

    if press is not None:
        ds_tmp['Pair'] = (['pair_time','eta_rho', 'xi_rho'], press)
        ds_tmp['Pair'].attrs = {'long_name':"surface air pressure",
                                 'units':"millibar",
                                 'field':"surface air pressure",
                                 'time':"pair_time",
                                 'coordinates':'lon_rho lat_rho' }
        
    ds_tmp.sms_time.encoding['units'] = 'days since 1900-01-01'
    ds_tmp.pair_time.encoding['units'] = 'days since 1900-01-01'

    if track_lon is not None:
        ds_tmp['track_lon'] = (['sms_time'], track_lon)
        ds_tmp['track_lat'] = (['sms_time'], track_lat)

    return ds_tmp