#%%
from pyarrow import fs, schema, string
from pyarrow.dataset import FilenamePartitioning
import pyarrow.dataset as ds
import boto3
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
import matplotlib.pylab as plt
import matplotlib.colors as colors
from scipy import optimize
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation as R
import io
from PIL import Image, ImageFilter
from mats_l1_processing.read_parquet_functions import *
from mats_l1_processing.L1_calibration_functions import *
from mats_l1_processing import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_utils.geolocation.coordinates import *
from mats_utils.geolocation import satellite as satellite
from mats_utils.rawdata.calibration import calibrate_dataframe
from skyfield import api as sfapi
import skyfield.sgp4lib as sgp4lib
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from pandas import concat, DataFrame, read_pickle

#%%
starttime = DT.datetime(2023,3, 28, 0, 0, 0)   
endtime = DT.datetime(2023, 3, 28, 10, 0, 0)
#starttime = datetime(2023,2, 13, 0, 30, 0)
#endtime = datetime(2023, 2, 13, 12, 45, 0)
ccd_data = read_MATS_data(starttime, endtime,filter=None,level='1a',version='0.5')

# #%%
#calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/Donal_stuff/calibration_data.toml"    
#instrument = Instrument(calibration_file)

# #%%
#l1b_data = calibrate_dataframe(ccd_data,instrument)
# %%
ccd_data_cal = read_MATS_data(starttime, endtime,filter=None,level='1b',version='0.4')

#%%

print('tmp')