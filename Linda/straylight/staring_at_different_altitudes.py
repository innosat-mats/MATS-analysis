#%%
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
from database_generation.experimental_utils import plot_CCDimage
import matplotlib.pyplot as plt
from mats_utils.rawdata.calibration import calibrate_dataframe

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
print(path)

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
# Calibrate the data
calibrate=False
if calibrate:
    if not 'instrument' in locals():
        instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

    dataframe_cal=calibrate_dataframe(dataframe, instrument, debug_outputs=True)

#%%


# start_time = DT.datetime(2022, 11, 28, 12, 0, 0)
# stop_time = DT.datetime(2022, 11, 29, 0, 0, 0)
# dataframe = read_MATS_data(start_time,stop_time,version='0.7',level='1a')

# %% [markdown]
# 14:00-17:00	2022.11.28	Rayleigh step	Get statistics on rayleigh atmosphere for calibration	Do images binned calibration looking at 60 km	Take WDW 7 images with binned calibration bu 12 sec exposure intervall		3080		Look at 60 km			20
# 17:00-20:00	2022.11.28	Rayleigh step		Do images binned calibration looking at 90 km	Take WDW 7 images with binned calibration bu 12 sec exposure intervall		3080		Look at 90 km			20
# 20:00-22:00	2022.11.28	Rayleigh step		Do images binned calibration looking at 120 km	Take WDW 7 images with binned calibration bu 12 sec exposure intervall		3080		Look at 120 km			20
# %%



# Define your datetime objects and localize them to UTC
start_60km = pd.Timestamp(DT.datetime(2022, 11, 28, 14, 0, 0), tz='UTC')
end_60km = pd.Timestamp(DT.datetime(2022, 11, 28, 17, 0, 0), tz='UTC')
start_90km = pd.Timestamp(DT.datetime(2022, 11, 28, 17, 0, 0), tz='UTC')
end_90km = pd.Timestamp(DT.datetime(2022, 11, 28, 20, 0, 0), tz='UTC')
start_120km = pd.Timestamp(DT.datetime(2022, 11, 28, 20, 0, 0), tz='UTC')
end_120km = pd.Timestamp(DT.datetime(2022, 11, 28, 22, 0, 0), tz='UTC')
df=dataframe

read_with_mats_read=True
if read_with_mats_read:
    df = read_MATS_data(start_60km, end_120km, version='0.7', level='1a')

#%%
# Convert your dataframe's datetime column to timezone-aware
df['TMHeaderTime_'] = pd.to_datetime(df['TMHeaderTime'], utc=True)

# Filter the dataframe
df_60km = df[(df['TMHeaderTime_'] > start_60km) & (df['TMHeaderTime_'] < end_60km)]
df_90km = df[(df['TMHeaderTime_'] > start_90km) & (df['TMHeaderTime_'] < end_90km)]
df_120km = df[(df['TMHeaderTime_'] > start_120km) & (df['TMHeaderTime_'] < end_120km)]
df_altitudes= [df_60km, df_90km, df_120km]
titles_km = ['60 km', '90 km', '120 km']

# Select based on SZA and latitude
#df_60km = df_60km[(df_60km['SZA'] > 90) & (df_60km['SZA'] < 100)]
#df_90km = df_90km[(df_90km['SZA'] > 90) & (df_90km['SZA'] < 100)]
#df_120km = df_120km[(df_120km['SZA'] > 90) & (df_120km['SZA'] < 100)]
# Select based on SZA and latitude
df_altitudes= df_altitudes

# %%

channels = ['IR1', 'IR2', 'IR3', 'IR4','UV1', 'UV2']

for channel in channels:

    fig, ax=plt.subplots(3,1)
    fig.suptitle(channel)
    for index, df in enumerate(df_altitudes):
        CCD=df[df.channel==channel].iloc[0]
        unitfix=CCD.TEXPMS/1000/10
        plot_CCDimage(CCD.IMAGE/unitfix, axis=ax[index],fig=fig,  title=titles_km[index])
        


    
#%%
for index, CCD in df_120km[:6].iterrows():
    fig, ax=plt.subplots(2,1)
    unitfix=CCD.TEXPMS/1000/10
    plot_CCDimage(CCD.IMAGE/unitfix, axis=ax[0],fig=fig, title=df.iloc[index].channel)

# %%
