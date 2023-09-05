#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd
import matplotlib.pyplot as plt
import time
from Keogram import CenterStrip
start_time = DT.datetime(2023,2,15,18,45,0)
stop_time = DT.datetime(2023,2,15,19,6,0)
channel = 'IR1'
numdays = stop_time-start_time #number of days
#%%
df = read_MATS_data(start_time,stop_time,filter={"TPlat":[50,90],'NROW': [0,400]},version='0.5',level='1b')
df.to_pickle('15feborb11')
"change latitude filter depending on if you want to look at north or south pole."

# %%
items = pd.read_pickle('15feborb11')
items = items[items['channel'] == channel]
pic = items.iloc[194]


ccdimage = pic['ImageCalibrated']
ccd_strip = ccdimage[:,22]

new_strip = CenterStrip(pic) #creates strip object
new_strip.makeVerticalStrip()


# %%
