#%% Import modules
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.geolocation import satellite as satellite
import matplotlib.pyplot as plt

#%%
# Select Time

start_time=DT.datetime(2023,1,11,18,0,0)
stop_time=DT.datetime(2023,1,12,6,0,0)

df=read_MATS_data(start_time, stop_time)
print(df.columns.tolist())


fig, ax=plt.subplots(1)
sensors=['HTR8B']#, 'HTR1B','HTR2A','HTR2A','HTR8A','HTR8B']
for sensor in sensors:
    line, =ax.plot(df.EXPDate, df[sensor])
    line.set_label(sensor)
plt.xlabel('Time')
plt.ylabel('Temperature [C]')
ax.legend()
# %%
