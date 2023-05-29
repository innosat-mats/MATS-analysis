#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from mats_utils.statistiscs.images_functions import create_imagecube
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, timezone
import warnings
import sys
import boto3
import re
import pyarrow.parquet as pq  # type: ignore
import pyarrow.dataset as ds
from pyarrow import fs
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from monitoring_functions import timeline_stat,timeline_plot,multi_timeline,temperatureCRBD_plot,temperatureHTR_plot,read_MATS_payload_data,PWRT_plot,PWRC_plot,PWRV_plot

pd.set_option('display.max_rows', 500)

#%%

# defining some parameters

# folders to store monitoring data
monitoring_folder = "/home/louis/MATS/MATS-Data/Monitoring"

sampling = 'custom'

custom_period = timedelta(minutes=2)


#%%

args = sys.argv
start_time = datetime.strptime(args[1],'%Y:%m:%d:%H:%M:%S')
stop_time = datetime.strptime(args[2],'%Y:%m:%d:%H:%M:%S')

if type(start_time) == type(None) or type(stop_time) == type(None):
    dt = datetime.now()
    dt = datetime(2023,4,25)
    stop_time = datetime(dt.year, dt.month, dt.day, 0, 0, 0)
    start_time = stop_time - timedelta(days=1)


print('===========================================')
print(f"Monitoring from {start_time} to {stop_time}")

# folders to store figures
data_folder = f"{monitoring_folder}/daily_monitoring_{start_time.strftime('%Y_%m_%d')}"
if not os.path.exists(data_folder):
        os.mkdir(data_folder)


#%%

if sampling == 'custom':
    start=start_time
    end=stop_time
    time_sampling = pd.date_range(start=start,
                    end=end,
                    periods=(end-start).total_seconds()/custom_period.total_seconds() + 1,tz=timezone.utc)
        
elif sampling == 'orbit':
    orbit_filter = {'CCDSEL': [1,7],'satlat':[-1.0,+1.0]}
    df = read_MATS_data(start_time, stop_time,pltfilter=orbit_filter,level='1a',version='0.5')
    df = df[~np.isnan(df['satlat'])].sort_values('EXPDate')
    time_sampling = [start_time.replace(tzinfo=timezone.utc)]
    for i in range(len(df)-1):
        if df.iloc[i]['satlat']<0.0 and df.iloc[i+1]['satlat']>0.0:
            time_sampling.append(df.iloc[i]['EXPDate'])  
    time_sampling.append(stop_time.replace(tzinfo=timezone.utc)) 
    
dataframes = []
dataframe_labels = []

try :
    print("Importing level 1b data")
    df1b = read_MATS_data(start_time, stop_time,level='1b',version='0.4')
    dataframes.append(df1b)
    dataframe_labels.append('l1b v0.4')
except :
    print('No level 1b data')

df1b = df1b.drop('ImageCalibrated', axis=1)


try :
    print("Importing level 1a data")
    df1a = read_MATS_data(start_time, stop_time,level='1a',version='0.5')
    dataframes.append(df1a)
    dataframe_labels.append('l1a v0.5')
except :
    print('No level 1a data')

df1a = df1a.drop(columns=['IMAGE','ImageData','id'], axis=1)


try :
    print("Importing level 0 data")
    df0 = read_MATS_data(start_time, stop_time,level='0',version='0.3')
    dataframes.append(df0)
    dataframe_labels.append('l0 v0.3')
except :
    print('No level 0 data')

df0 = df0.drop('ImageData', axis=1)
    
    
if len(dataframes)>0:
    multi_timeline(dataframes,dataframe_labels,time_sampling,data_folder=data_folder)

try:
    print(f"Plotting CRB-D temperatures")
    start = min(df0['EXPDate'])
    end = max(df0['EXPDate'])
    file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_CRBD_temp.png"
    temperatureCRBD_plot(df0,title='',file=file_path)
except:
    print(f"Unable to plot CRB-D temperatures from l0 v0.3")



try:
    print(f"Importing HTR temperature data")
    HTR_df = read_MATS_payload_data(start_time,stop_time,data_type='HTR')
    file_path = f"{data_folder}/HTR_temp.png"
    print(f"Plotting HTR temperatures")
    temperatureHTR_plot(HTR_df,file=file_path)
except:
    print(f"Unable to plot HTR temperatures")


try:
    print(f"Importing PWR data")
    PWR_df = read_MATS_payload_data(start_time,stop_time,data_type='PWR')
    print(f"Plotting PWR voltages")
    PWRV_plot(PWR_df,file=f"{data_folder}/PWR_voltage.png")
    print(f"Plotting PWR temperature")
    PWRT_plot(PWR_df,file=f"{data_folder}/PWR_temp.png")
    print(f"Plotting PWR currents")
    PWRC_plot(PWR_df,file=f"{data_folder}/PWR_current.png")
except:
    print(f"Unable to plot PWR temperatures")


# %%
