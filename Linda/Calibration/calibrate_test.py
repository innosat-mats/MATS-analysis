#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items


#%% 

def get_partial_calibration(CCDitems,instrument,varname='image_calibrated'):    
    for CCDitem in CCDitems:
        image_lsb, image_se_corrected, image_hot_pixel_corrected, image_bias_sub, image_linear,image_desmeared, image_dark_sub, image_calib_nonflipped, image_calib_flipped, image_calibrated, errors = L1_calibrate(CCDitem,instrument, return_steps=True)
    return locals()[varname]

#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument = Instrument(calibration_file)

#%% Select on explicit time
start_time = DT.datetime(2023, 5, 5, 1, 40)
stop_time = DT.datetime(2023, 5, 5, 1, 45)

#%%
df = read_MATS_data(start_time,stop_time,version='0.6',level='1a',dev=False)

#df_ph = read_MATS_PM_data(start_time,stop_time,version='0.1',level='1b')

#%%
CCDitems = dataframe_to_ccd_items(df)

#%%

caloutput=L1_calibrate(CCDitems[4],instrument)

# %%

outputlist= get_partial_calibration(CCDitems[0:3],instrument,varname='image_calibrated')
# %%
