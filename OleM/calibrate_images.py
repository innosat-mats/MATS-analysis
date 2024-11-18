#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data, add_ccd_item_attributes
import pandas as pd
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items, read_ccd_data, remove_faulty_rows,convert_image_data
from pyarrow import fs
import pyarrow.dataset as ds
import boto3
#%% 
#calibration_file ="/home/olemar/Projects/Universitetet/MATS/MATS-L1-processing/scripts/calibration_data.toml"    
#instrument = Instrument(calibration_file)

#%% Select on explicit time
start_time = DT.datetime(2022, 11, 27, 0, 0)
stop_time = DT.datetime(2022, 11,30, 0, 0)

#%%
#df = read_MATS_data(start_time,stop_time,version='0.2',level='1a',dev=False)


session = boto3.session.Session(profile_name="mats")
credentials = session.get_credentials()

s3 = fs.S3FileSystem(
    secret_key=credentials.secret_key,
    access_key=credentials.access_key,
    region=session.region_name,
    connect_timeout=10,
    session_token=credentials.token)

#%%
main_level = '1a'
version = '0.2'
subdir = '2022/11/28'
path = f"ops-payload-level{main_level}-v{version}" + "/" + subdir

dataframe = read_ccd_data(path,s3)

if dataframe.index.name == 'EXPDate':
    dataframe.reset_index(inplace=True)
    dataframe.set_index('TMHeaderTime', inplace=True)
    dataframe.sort_index(inplace=True)
    dataframe.reset_index(inplace=True)
else:
    dataframe.reset_index(drop=True, inplace=True)
    dataframe.set_index('TMHeaderTime', inplace=True)
    dataframe.sort_index(inplace=True)
    dataframe.reset_index(inplace=True)

add_ccd_item_attributes(dataframe)
convert_image_data(dataframe)
remove_faulty_rows(dataframe)    



#%%
#df_ph = read_MATS_PM_data(start_time,stop_time,version='0.1',level='1b')

#%%
CCDitems = dataframe_to_ccd_items(df)

#%%

L1_calibrate(CCDitems[0],instrument)
L1_calibrate(CCDitems[1],instrument)
L1_calibrate(CCDitems[2],instrument)
L1_calibrate(CCDitems[3],instrument)
L1_calibrate(CCDitems[4],instrument)
L1_calibrate(CCDitems[5],instrument)
L1_calibrate(CCDitems[6],instrument)



#%%
# for i in range(len(CCDitems)):
#     image_lsb, image_bias_sub, image_desmeared, image_dark_sub, image_calib_nonflipped, image_calibrated, errors = L1_calibrate(CCDitems[i],instrument)
# %%
#for i in range(len(CCDitems)):
I = 1
image_lsb, image_bias_sub, image_desmeared, image_dark_sub, image_calib_nonflipped, image_calibrated, errors = L1_calibrate(CCDitems[I],instrument)

# %%
from matplotlib import pyplot as plt

plt.imshow(image_lsb,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"]+ " lsb")
plt.colorbar()
plt.show()

plt.imshow(image_bias_sub,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"] + " bias subtracted")
plt.colorbar()
plt.show()

plt.imshow(image_desmeared,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"]+ " desmeared")
plt.colorbar()
plt.show()

plt.imshow(image_dark_sub,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"]+ " dark subtracted")
plt.colorbar()
plt.show()

plt.imshow(image_calibrated,origin='lower')
plt.title('Channel ' + CCDitems[I]["channel"] + " flatfied and abs calibrated")
plt.colorbar()
plt.show()

# %%
