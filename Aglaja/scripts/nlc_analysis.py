import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import argparse
from matplotlib import pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
import pdb
import netCDF4 as nc



###############
## TEST NLC ##
###############
# python get_data.py --channel UV2 --start_time 2023 2 11 0 0 --stop_time 2023 2 12 0 0 --ncdf_out UV2_test.nc --version 0.6
# read .nc file
data = nc.Dataset('UV2_test.nc', 'r')

print(data.variables.keys())

