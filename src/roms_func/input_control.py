import xarray as xr
import subprocess
import pandas as pd

def make_infile_from_files( fp_template = './roms.in.template',
                            fp_out = './roms.in',
                            fp_grd = './roms_grd.nc',
                            fp_frc = './roms_frc.nc',
                            dt = 10):

    # Open datasets
    ds_grd = xr.open_dataset(fp_grd)
    ds_frc = xr.open_dataset(fp_frc)

    # Figure out grid shape
    Lm = len(ds_grd.xi_rho) - 2
    Mm = len(ds_grd.eta_rho) - 2

    # Figure out forcing times
    time_ref = pd.to_datetime(ds_frc.sms_time[0].values)
    time_end = pd.to_datetime(ds_frc.sms_time[-1].values)
    time_diff = (time_end - time_ref).total_seconds()
    ntimes = int( time_diff / dt )

    # Make the infile
    make_infile( fp_template, fp_out, Lm, Mm, time_ref, ntimes )
    

def make_infile( fp_template = './roms.in.template',
                 fp_out = './roms.in',
                 Lm = None, Mm = None,
                 time_ref = None, ntimes = None):

    # Copy template to output file. Changes will be inplace
    subprocess.run( f'cp {fp_template} {fp_out}', shell=True )
    
    if Lm is not None:
        sed_str = f"\'s/DBREPLACE_Lm/{Lm}/g\'"
        subprocess.run(f'sed -i -e {sed_str} {fp_out} {fp_out}', shell=True)
    if Mm is not None:
        sed_str = f"\'s/DBREPLACE_Mm/{Mm}/g\'"
        subprocess.run(f'sed -i -e {sed_str} {fp_out} {fp_out}', shell=True)
    if time_ref is not None:
        time_ref = time_ref.strftime("%Y%m%d") + '.00'
        sed_str = f"\'s/DBREPLACE_TIMEREF/{time_ref}/g\'"
        subprocess.run(f'sed -i -e {sed_str} {fp_out} {fp_out}', shell=True)
    if ntimes is not None:
        sed_str = f"\'s/DBREPLACE_NTIMES/{ntimes}/g\'"
        subprocess.run(f'sed -i -e {sed_str} {fp_out} {fp_out}', shell=True)

    return