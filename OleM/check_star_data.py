#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from matplotlib import pyplot as plt 
#%% 
calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-L1-processing/tests/calibration_data_test.toml"    
instrument = Instrument(calibration_file)

#%% Select on explicit time
start_time = DT.datetime(2023, 10, 18, 0, 0)
stop_time = DT.datetime(2023, 10, 18, 23, 0)

df = read_MATS_data(start_time,stop_time,version='0.3',level='0/CCD')

df.to_pickle('star3.pkl')
#%%
df = pd.read_pickle('star3.pkl')
df = df[df.NCBINCCDColumns == 1]
# %%
for i in range(len(df)):
    plt.imshow(df.iloc[i].IMAGE)
    plt.clim([0,1000])
    plt.title(str(df.iloc[i].CCDSEL))
    plt.colorbar()
    plt.show()
# %%
