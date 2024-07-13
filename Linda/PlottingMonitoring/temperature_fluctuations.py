
#%% Import modules
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.geolocation import satellite as satellite
import matplotlib.pyplot as plt

#%%
# Select Time

starttime=DT.datetime(2023,2,11,14,0,0)
endtime=DT.datetime(2023,2,15,14,30,0)


df = read_MATS_data(starttime, endtime,filter=None,level='1a',version='0.6')

print(df.columns.tolist())

#%%

fig, ax=plt.subplots(1)
sensors=['HTR1A','HTR1B','HTR2A','HTR2B','HTR8A','HTR8B']
for sensor in sensors:
    line, =ax.plot(df.EXPDate, df[sensor])
    line.set_label(sensor)
plt.xlabel('Time')
plt.ylabel('Temperature [C]')
plt.ylim([-25, -17])
ax.legend()
# %%
df['HTR8A'].plot()

# %%
