#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_l1_processing.read_in_functions import read_CCDitems 
import numpy as np
import pickle
import matplotlib.pyplot as plt
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.rawdata.calibration import calibrate_dataframe

#%%

# # #%% Select on explicit time NLC
#start_time = DT.datetime(2023, 2, 2, 19, 38, 0) 
#stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)
start_time = DT.datetime(2023, 2, 9, 18, 59, 39) 
stop_time = DT.datetime(2023, 2, 9, 19, 1, 44)
#start_time = DT.datetime(2023, 3, 14, 19, 38, 0) 
#stop_time = DT.datetime(2023, 3, 14, 19, 50, 0)
level='1b'
dfold = read_MATS_data(start_time,stop_time,version='0.9',level=level,dev=False)
print('Old version read')
dfnew = read_MATS_data(start_time,stop_time,version='1.0',level=level,dev=False)
print('New version read')
if level=='1b':
    ImageName='ImageCalibrated'
else:
    ImageName='IMAGE'



ir1old=dfold[dfold.channel=='IR1']
ir2old=dfold[dfold.channel=='IR2']
ir3old=dfold[dfold.channel=='IR3']
ir4old=dfold[dfold.channel=='IR4']
uv1old=dfold[dfold.channel=='UV1']
uv2old=dfold[dfold.channel=='UV2']

ir1new=dfnew[dfnew.channel=='IR1']
ir2new=dfnew[dfnew.channel=='IR2']
ir3new=dfnew[dfnew.channel=='IR3']
ir4new=dfnew[dfnew.channel=='IR4']
uv1new=dfnew[dfnew.channel=='UV1']
uv2new=dfnew[dfnew.channel=='UV2']
#create tuples with the dataframes
ir1=(ir1old,ir1new)
ir2=(ir2old,ir2new)
ir3=(ir3old,ir3new)
ir4=(ir4old,ir4new)
uv1=(uv1old,uv1new)
uv2=(uv2old,uv2new)



#%%
#Plot and compare the two versions
channels = [ir1, ir2, ir3, ir4, uv1, uv2]
fig, axs = plt.subplots(6, 4, figsize=(14, 10))
for ind, channel in enumerate(channels):
    dfentryold = channel[0].iloc[0]
    dfentrynew = channel[1].iloc[0]

    if dfentryold['TMHeaderTime'] != dfentrynew['TMHeaderTime']:
        print('The two dataframes are not aligned')
        print('Time old:', dfentryold['TMHeaderTime'])
        print('Time new:', dfentrynew['TMHeaderTime'])
    else:
        plot_CCDimage(dfentryold[ImageName], title=dfentryold.channel+ ' old', fig=fig, axis=axs[ind, 0])
        plot_CCDimage(dfentrynew[ImageName], title=dfentrynew.channel + ' new', fig=fig, axis=axs[ind, 1])
        plot_CCDimage(dfentryold[ImageName]-dfentrynew[ImageName], title='Difference old-new', fig=fig, axis=axs[ind, 2])
        meanratio=(dfentrynew[ImageName]/dfentryold[ImageName]).mean()
        plot_CCDimage(dfentrynew[ImageName]/dfentryold[ImageName], title = f' new/old, Mean={meanratio:.3f}', fig=fig, axis=axs[ind, 3])
        

plt.tight_layout()


# %%
