# script for studying the number of saturated pixels in NADIR images as a function of solare zenith angle 

#%% Import modules
#%matplotlib qt5
from Nadir_artifact_correction import *
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import warnings
  








# %%
# import data
start_time_mask = DT.datetime(2023, 4, 13, 3, 30, 0)
stop_time_mask = DT.datetime(2023, 4, 13, 4, 30, 0)

start_time_mask_1day = DT.datetime(2023, 4, 13, 0, 0, 0)
stop_time_mask_1day = DT.datetime(2023, 4, 14, 0, 0, 0)

start_time = DT.datetime(2023, 4, 14, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 14, 2, 0, 0)

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}

df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
df1a_tot.sort_values('EXPDate')
df1a_tot = df1a_tot[~np.isnan(df1a_tot['satlat'])]
print(len(df1a_tot))
df1a = df1a_tot

df1a_tot_mask= read_MATS_data(start_time_mask, stop_time_mask,filter,level='1a',version='0.5')
df1a_tot_mask.sort_values('EXPDate')
df1a_tot_mask = df1a_tot_mask[~np.isnan(df1a_tot_mask['satlat'])]
print(len(df1a_tot_mask))
df1a_mask = df1a_tot_mask

df1a_tot_mask_1day= read_MATS_data(start_time_mask_1day, stop_time_mask_1day,filter,level='1a',version='0.5')
df1a_tot_mask_1day.sort_values('EXPDate')
df1a_tot_mask_1day = df1a_tot_mask_1day[~np.isnan(df1a_tot_mask_1day['satlat'])]
print(len(df1a_tot_mask_1day))
df1a_mask_1day = df1a_tot_mask_1day

# displaying keys
pd.set_option('display.max_rows', 500)

#%%
# calculate masks

azimuths = np.concatenate((np.linspace(-105,-85.5,41),np.linspace(-85,-74,111),np.linspace(-73.5,-70,8)))

azimuth_masks = azimuth_bias_mask2(df1a_mask,-6000,azimuths)

azimuth_masks_5deg = azimuth_bias_mask2(df1a_mask,-6000,azimuths)

azimuth_masks_1day = azimuth_bias_mask2(df1a_mask_1day,-6000,azimuths)

#%%
# correcting images
df_corr = azimuth_corr_mask(df1a,azimuth_masks)

df_corr_5deg = azimuth_corr_mask(df1a,azimuth_masks_5deg)

df_corr_1day = azimuth_corr_mask(df1a,azimuth_masks_1day)


#%%
# plotting masks

angles = [-100,-90,-85,-80,-79,-78,-77-76,-75,-70]
angles = azimuths[::10]

azimuth_masks_plot(azimuth_masks,angles)
azimuth_masks_plot(azimuth_masks_5deg,angles)
azimuth_masks_plot(azimuth_masks_1day,angles)

# %%
vmax = 32000.0
vmin = 0.0

i = 568


plt.figure()
plt.imshow(df1a.iloc[i]['IMAGE'],vmin=vmin,vmax=vmax,origin='lower')
plt.show()

plt.figure()
plt.imshow(df_corr.iloc[i]['IMAGE'],vmin=vmin,vmax=vmax,origin='lower')
plt.show()

plt.figure()
plt.imshow(df_corr_5deg.iloc[i]['IMAGE'],vmin=vmin,vmax=vmax,origin='lower')
plt.show()

plt.figure()
plt.imshow(df_corr_1day.iloc[i]['IMAGE'],vmin=vmin,vmax=vmax,origin='lower')
plt.show()


# %%
azimuth_masks_1day = azimuth_bias_mask2(df1a_mask,-6000,np.linspace(-80,-74,61))
# %%
