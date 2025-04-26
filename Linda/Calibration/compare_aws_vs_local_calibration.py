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
#start_time = DT.datetime(2023, 3, 14, 19, 38, 0) 
#stop_time = DT.datetime(2023, 3, 14, 19, 40, 0)
start_time = DT.datetime(2023, 2, 9, 18, 59, 39) 
stop_time = DT.datetime(2023, 2, 9, 19, 3, 44)


df_loc_l1a = read_MATS_data(start_time,stop_time,version='1.0',level='1a',dev=False)
print('Level 1a read')

#%%
if not 'instrument' in locals():
    instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')
df_loc=calibrate_dataframe(df_loc_l1a, instrument, debug_outputs=False)
print('Level 1a calibrated')
#%%
df_aws = read_MATS_data(start_time,stop_time,version='1.0',level='1b',dev=False)
print('Level 1b read')

ImageName='ImageCalibrated'

#%%

ir1_loc=df_loc[df_loc.channel=='IR1']
ir2_loc=df_loc[df_loc.channel=='IR2']
ir3_loc=df_loc[df_loc.channel=='IR3']
ir4_loc=df_loc[df_loc.channel=='IR4']
uv1_loc=df_loc[df_loc.channel=='UV1']
uv2_loc=df_loc[df_loc.channel=='UV2']

ir1_aws=df_aws[df_aws.channel=='IR1']
ir2_aws=df_aws[df_aws.channel=='IR2']
ir3_aws=df_aws[df_aws.channel=='IR3']
ir4_aws=df_aws[df_aws.channel=='IR4']
uv1_aws=df_aws[df_aws.channel=='UV1']
uv2_aws=df_aws[df_aws.channel=='UV2']
#create tuples with the dataframes
ir1=(ir1_loc,ir1_aws)
ir2=(ir2_loc,ir2_aws)
ir3=(ir3_loc,ir3_aws)
ir4=(ir4_loc,ir4_aws)
uv1=(uv1_loc,uv1_aws)
uv2=(uv2_loc,uv2_aws)



#%%
#Plot and compare the two versions
channels = [ir1, ir2, ir3, ir4, uv1, uv2]
fig, axs = plt.subplots(6, 4, figsize=(18, 10))
for ind, channel in enumerate(channels):
    dfentry_loc = channel[0].iloc[0] 
    dfentry_aws = channel[1].iloc[0]
    

    if dfentry_loc['TMHeaderTime'] != dfentry_aws['TMHeaderTime']:
        print('The two dataframes are not aligned')
        print('Time _loc:', dfentry_loc['TMHeaderTime'])
        print('Time _aws:', dfentry_aws['TMHeaderTime'])
    else:
        
        plot_CCDimage(dfentry_loc[ImageName], title=dfentry_loc.channel+ ' _loc '+str(dfentry_loc.TMHeaderTime), fig=fig, axis=axs[ind, 0])
        plot_CCDimage(dfentry_aws[ImageName], title=dfentry_aws.channel + ' _aws '+str(dfentry_aws.TMHeaderTime), fig=fig, axis=axs[ind, 1])
        plot_CCDimage(dfentry_loc[ImageName]-dfentry_aws[ImageName], title='Difference _loc-_aws', fig=fig, axis=axs[ind, 2])
        meanratio=(dfentry_aws[ImageName]/dfentry_loc[ImageName]).mean()
        plot_CCDimage(dfentry_aws[ImageName]/dfentry_loc[ImageName], title = f' _aws/_loc, Mean={meanratio:.3f}', fig=fig, axis=axs[ind, 3])
        

plt.tight_layout()


# %%
