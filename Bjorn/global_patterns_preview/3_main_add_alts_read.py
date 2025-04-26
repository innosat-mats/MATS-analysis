# Reads and saves output from main_adda_alts
#
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from datetime import date, timedelta
from mats_utils.plotting.plotCCD import all_channels_plot
from mats_utils.daily_preview.temp_nadirs import NADIR_geolocation, average_stacking
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import ceil
import cartopy.crs as ccrs
from multiprocessing import Manager
import xarray as xr
from mats_utils.geolocation import coordinates
from mats_utils.geolocation.altitude_correction import rows_to_altitudes
import os
import multiprocessing


# load in 
monthstrings=["022023_1","022023_2","022023_3","022023_4"]
ccdnos=[3]
#pathname='/media/waves/AVAGO/data/MATS/pandas_csv/7070/'
pathname='/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/'

for ccdno in ccdnos:
    for month in monthstrings:

        # to check path
        path = f'{pathname}CCDSEL{ccdno}/alts_80to90/finished/'

        # Check whether the specified path exists or not
        isExist = os.path.exists(path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path)
            print("! DIRECTORY CREATED !")
        
        location=f'{pathname}CCDSEL{ccdno}/alts_80to90/{month}/'
        CCDs = pd.concat([pd.read_pickle(f'{location}/'+x) for x in os.listdir(f'{location}')])

        # save
        CCDs.to_pickle(f'{pathname}CCDSEL{ccdno}/alts_80to90/finished/{month}.pk1')

        print(f'{month} finished')