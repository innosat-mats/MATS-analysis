#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd 
import sys
import numpy as np
import matplotlib.pyplot as plt 
from mats_utils.geolocation.coordinates import TPpos
#%%
from Keogram import plotKeogram
# %%
start_time = DT.datetime(2023,2,16,18,30,0)
stop_time = DT.datetime(2023,2,18,18,30,0)
timedelta = stop_time-start_time #number of days
print(timedelta.days)

# %%

for day in range(timedelta.days):
    fig, axs = plt.subplots(nrows=len(timedelta.days), ncols=1)
    df = read_MATS_data(start_time,stop_time)
