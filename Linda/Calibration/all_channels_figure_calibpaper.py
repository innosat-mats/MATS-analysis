#%%
# Plot the cahnnels in for figure in calibration paper
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
from mats_utils.plotting.plotCCD import plot_image

def converttostarcal(df_oldcalib, field= 'ImageCalibrated'):
    # convert to full star calibration
    # Star calibration values
    import pandas as pd
    star_calibration = [10.1, 2.99, 21.2, 26.9, 60.9, 23.6]

    # Absolute star and relative lab calibration values
    absolute_star_calibration = [10.2, 2.99, 17.8, 25.7, 57.8, 23.6]
    # Channels
    channels = ['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']

    # Calculate the ratios
    ratios = [star_calibration[i] / absolute_star_calibration[i] for i in range(len(star_calibration))]

    # Create a DataFrame
    df_calibconvert = pd.DataFrame([ratios], columns=channels)

    #Convert the calibration
    df=df_oldcalib.copy()
    print(' mean value:'+ df_oldcalib.iloc[3].channel, df_oldcalib.loc[3][field].mean())
    #
    for index, row in df_oldcalib.iterrows():
        channel = row['channel']
        df.at[index, field] = row[field] * df_calibconvert[channel].values[0]
    print(' mean value:'+ df.iloc[3].channel, df.loc[3].ImageCalibrated.mean())

    print('WARNING: Calibration is altered, this should be removed when data is reprocessed')
    return df


#%%

#
# # #%% Select on explicit time NLC
#start_time = DT.datetime(2023, 2, 2, 19, 38, 0) 
#stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)
print('**************** The figure in the paper must be changed. It was from (2023, 2, 9, 19, 38, 0) to (2023, 2, 9, 19, 50, 0), when CROPC was used. ****************')

start_time = DT.datetime(2023, 2, 9, 20, 30, 0) 
stop_time = DT.datetime(2023, 2, 9, 20, 35, 0)

#df_oldcalib= read_MATS_data(start_time,stop_time,version='0.9',level='1b',dev=False)
#df=converttostarcal(df_oldcalib)
df= read_MATS_data(start_time,stop_time,version='1.0',level='1b',dev=False)
print('schedule',df.schedule_name.unique())

#%%
n=0
fig, axs = plt.subplots(3, 2, figsize=(16, 10))
sp=plot_CCDimage(df[df.channel=='IR1'].iloc[n].ImageCalibrated,fig=fig, axis=axs[0,0], title='IR1')
clim = sp.get_clim()
sp=plot_CCDimage(df[df.channel=='IR2'].iloc[n].ImageCalibrated,fig=fig, axis=axs[0,1], title='IR2', clim=clim)

sp=plot_CCDimage(df[df.channel=='IR3'].iloc[n].ImageCalibrated,fig=fig, axis=axs[1,0], title='IR3')
clim = sp.get_clim()
sp=plot_CCDimage(df[df.channel=='IR4'].iloc[n].ImageCalibrated,fig=fig, axis=axs[1,1], title='IR4', clim=clim)

sp=plot_CCDimage(df[df.channel=='UV1'].iloc[n].ImageCalibrated,fig=fig, axis=axs[2,0], title='UV1')
clim = sp.get_clim()
sp=plot_CCDimage(df[df.channel=='UV2'].iloc[n].ImageCalibrated,fig=fig, axis=axs[2,1], title='UV2', clim=clim)


fig.suptitle('Observation of dayglow and NLC ')
plt.tight_layout()
#fig.savefig('../output/observation_dayglow_NLC_20230202.png')


# %%
# def plot_image(CCD, ax=None, fig=None, outpath=None,
#                nstd=2, cmap='inferno', ranges=None,
#                optimal_range=False, format='png', save=True,
#                fontsize=10, TPheights=True, image_field='None')
n=42
fig, axs = plt.subplots(3, 2, figsize=(8, 7))
ax, subfig, im=plot_image(df[df.channel=='IR1'].iloc[n],ax=axs[0,0],fig=fig,  save=False, image_field='ImageCalibrated', title='IR1')
clim = im.get_clim()
fig.colorbar(im, ax=axs[0,0])
axs[0,0].set_xlabel('Column number')
axs[0,0].set_ylabel('Row number')

ax, subfig, im=plot_image(df[df.channel=='IR2'].iloc[n],ax=axs[0,1],fig=fig,  save=False, image_field='ImageCalibrated', ranges=clim, title='IR2')  
fig.colorbar(im, ax=axs[0,1])
axs[0,1].set_xlabel('Column number')
axs[0,1].set_ylabel('Row number')

ax, subfig, im=plot_image(df[df.channel=='IR3'].iloc[n],ax=axs[1,0],fig=fig,  save=False, image_field='ImageCalibrated', title='IR3')
clim = im.get_clim()
fig.colorbar(im, ax=axs[1,0])
axs[1,0].set_xlabel('Column number')
axs[1,0].set_ylabel('Row number')

ax, subfig, im=plot_image(df[df.channel=='IR4'].iloc[n],ax=axs[1,1],fig=fig,  save=False, image_field='ImageCalibrated', ranges=clim, title='IR4')  
fig.colorbar(im, ax=axs[1,1])
axs[1,1].set_xlabel('Column number')
axs[1,1].set_ylabel('Row number')

ax, subfig, im=plot_image(df[df.channel=='UV2'].iloc[n],ax=axs[2,1],fig=fig,  save=False, image_field='ImageCalibrated', ranges=[0, 1000], title='UV2') 
fig.colorbar(im, ax=axs[2,1])
axs[2,1].set_xlabel('Column number')
axs[2,1].set_ylabel('Row number')

clim = im.get_clim()
ax, subfig, im=plot_image(df[df.channel=='UV1'].iloc[n],ax=axs[2,0],fig=fig,  save=False, image_field='ImageCalibrated', ranges=clim, title='UV1')
fig.colorbar(im, ax=axs[2,0])
axs[2,0].set_xlabel('Column number')
axs[2,0].set_ylabel('Row number')


fig.suptitle('Example of dayglow and NLC observation [10¹² photons m⁻² sr⁻¹ s⁻¹ nm⁻¹]')
plt.tight_layout()
fig.savefig('../output/observation_dayglow_NLC_20230209.png')

## %%

# %%
