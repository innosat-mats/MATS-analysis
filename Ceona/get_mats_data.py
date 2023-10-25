# %%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd 
import sys
import numpy as np
import matplotlib.pyplot as plt 
from mats_utils.geolocation.coordinates import TPpos, satpos
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot, all_channels_plot
from datetime import timedelta

start_time = DT.datetime(2023,3,2,18,58,0)
stop_time = DT.datetime(2023,3,2,19,22,0)
numdays = stop_time - start_time
# %%
df = read_MATS_data(start_time,stop_time,filter={"TPlat":[45,90], 'NROW': [0,400]},version=0.5,level='1b')
# If you want to select only IR1 :  'CCDSEL': [1, 1]
# filter={'afsTPLongLatGeod'[1],:[-90,90]}

#%%
df = df[df['channel']=='IR1']
df.to_pickle('3marstest')
#%%
df = pd.read_pickle('data')
IR1 = df[df['channel'] == 'IR1']
IR2 = df[df['channel']== 'IR2']
#print(IR1.keys)
#mask_night = (df['EXPDate']>= pd.to_datetime(timedelta(hours=19))) & (df['EXPDate']<= pd.to_datetime(timedelta(hours=23)))
#CCD_nighttime = IR2.loc[mask_night]

#CCD_nighttime = IR2.between_time('19:00','20:00')

#print(type(IR2['IMAGE'].iloc[0]))

# %%
pic = IR1.iloc[39]
#plt.pcolormesh(IR2['ImageCalibrated'].iloc[0],vmin=-20, vmax=260)
plt.pcolormesh(pic['ImageCalibrated'])

# variables

# Use TLE to get some positions from time of measurement
#satlat,satlon,satLT,nadir_sza,nadir_mza,TPlat,TPlon,TPLT,TPsza,TPssa = satellite.get_position(CCDitems['EXPDate'][4])


# %%
