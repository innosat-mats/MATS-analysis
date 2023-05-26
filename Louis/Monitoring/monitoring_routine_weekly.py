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

dt = datetime.now() - timedelta(days=7*6)
day = datetime(dt.year, dt.month, dt.day, 0, 0, 0)
week_start = day - timedelta(days=day.weekday()+7)
week_end = week_start + timedelta(days=7)


data_folder = f"{monitoring_folder}/weekly_monitoring{week_start.strftime('%Y_%m_%d')}_{week_end.strftime('%Y_%m_%d')}"
if not os.path.exists(data_folder):
        os.mkdir(data_folder)



#%%

try:
    print(f"Importing HTR temperature data")
    HTR_df = read_MATS_payload_data(week_start,week_end,data_type='HTR')
    file_path = f"{data_folder}/HTR_temp.png"
    print(f"Plotting HTR temperatures")
    temperatureHTR_plot(HTR_df,file=file_path)
except:
    print(f"Unable to plot HTR temperatures")


try:
    print(f"Importing PWR data")
    PWR_df = read_MATS_payload_data(week_start,week_end,data_type='PWR')
    print(f"Plotting PWR voltages")
    PWRV_plot(PWR_df,file=f"{data_folder}/PWR_voltage.png")
    print(f"Plotting PWR temperature")
    PWRT_plot(PWR_df,file=f"{data_folder}/PWR_temp.png")
    print(f"Plotting PWR currents")
    PWRC_plot(PWR_df,file=f"{data_folder}/PWR_current.png")
except:
    print(f"Unable to plot HTR temperatures")
# %%
