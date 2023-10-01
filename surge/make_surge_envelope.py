import argparse
import numpy as np
import xarray as xr
import warnings
import os
from pathlib import Path
from roms_func import postprocessing

# Filter warnings
warnings.filterwarnings("ignore")

description = (' Create surge envelope from roms history file ')
                 
# Parse input arguments
parser = argparse.ArgumentParser(
                prog='make_surge_envelope',
                description= description,
                epilog='Author: David Byrne, Woodwell Climate Research Center',
                formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument('-proj', type=str, help='Project name/directory', default='.')
parser.add_argument('-f', type=str, help='Name of ROMS history output file', default='roms_his.nc')
parser.add_argument('-o', type=str, help='Name of output .nc file', default='roms_zenv.nc')
args = parser.parse_args()

#- - - - - - - - - - - - - - - -#
# MAKE PLOTS
#- - - - - - - - - - - - - - - -#

postprocessing.surge_envelope_from_file( os.path.join( args.proj, args.f) , 
                                         os.path.join( args.proj, args.o ) )