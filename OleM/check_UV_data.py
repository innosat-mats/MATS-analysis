#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items

# #%% 
# calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-L1-processing/scripts/calibration_data.toml"    
# instrument = Instrument(calibration_file)

# #%% Select on explicit time
# start_time = DT.datetime(2023, 5, 8, 0, 0)
# stop_time = DT.datetime(2023, 5, 8, 23, 0)

# df = read_MATS_data(start_time,stop_time,version='0.6',level='1a')

#%%
start_time=DT.datetime(2023,5,7,1,0,0)
stop_time=DT.datetime(2023,5,7,2,0,0)
filter_UV2={'CCDSEL': [6,6]}
df=read_MATS_data(start_time, stop_time, filter_UV2, level='1b')

#%%
start_time=DT.datetime(2023,5,7,1,0,0)
stop_time=DT.datetime(2023,5,7,2,0,0)
filter_UV2={'CCDSEL': [6,6]}
df=read_MATS_data(start_time, stop_time, filter_UV2, level='1a')
# %%
