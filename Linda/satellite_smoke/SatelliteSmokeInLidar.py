#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
import numpy as np
import matplotlib.pyplot as plt

from database_generation.experimental_utils import plot_CCDimage
from mats_utils.retrieval.error_flags import remove_flagged_images
from mats_utils.retrieval.hot_pix import create_hot_pix_map_one_channel, create_all_hot_pix_maps, hot_pix_removal_several_channels, hot_pix_removal_one_channel

import pickle



#%%
# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

readinrawdata=True
readindfwithhotpix=False
readincleandataframe=False
createhotpixmap=True

start_time = DT.datetime(2025, 2, 19, 0, 0, 0)
stop_time = DT.datetime(2025, 2, 20, 23, 59, 59)
level='1b'
version='1.0'
name='df_all'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')+'_'+level+'v'+version
# filter
#filter={'CCDSEL': [1,4], 'TPlon': [1,4]}

#%%

    # read in measurements
if readinrawdata:
    df = read_MATS_data(start_time, stop_time,level=level,version=version)
    #df = read_MATS_data(start_time, stop_time,filter,level='1b',version='1.0')
    
    df.to_pickle(data_folder+name+'.pkl')
else:
    if readindfwithhotpix:
        df = pd.read_pickle(data_folder+'df_hot_pix_'+name+'.pkl')
    else:
        df = pd.read_pickle(data_folder+name+'.pkl')


#%%
if level=='1b':
    df_clean, mask=remove_flagged_images(df, [8, 9], return_mask=True) # Remove images where desmearing malfunctioned
    df_flagged = df[mask].reset_index(drop=True)
else:
    df_clean = df.copy()
    

#%%

if createhotpixmap:
    hot_pix_maps_dict=create_all_hot_pix_maps(df_clean, level=level)
    with open(data_folder+'hot_pixel_maps.pkl', 'wb') as f:
        pickle.dump(hot_pix_maps_dict, f)
else:
    #read in picked data
    with open(data_folder+'hot_pixel_maps.pkl', 'rb') as f:
        hot_pix_maps_dict = pickle.load(f)

#%%
#plot hot pixel maps for all channels
for channel, hot_pix_map in hot_pix_maps_dict.items():
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_CCDimage(hot_pix_map, fig=fig, axis=ax, title=f"Hot Pixel Map - {channel}")
    plt.show()

#%%
# fix ir2 hot pixel map since it is not working properly
def apply_running_mean_filter(img, filter_rows, filter_cols):
    from scipy.ndimage import uniform_filter
    running_mean = uniform_filter(img, size=(filter_rows, filter_cols))
    image = img.copy().astype(float)
    img_fixed = image - running_mean
    return img_fixed

hot_pix_map_new=apply_running_mean_filter(hot_pix_maps_dict['IR2'], filter_rows=2, filter_cols=5)
plot_CCDimage(hot_pix_map_new, title='Fixed Hot Pixel Map - IR2')
#hot_pix_maps_dict['IR2'] = hot_pix_map_new

#%%
channel='IR2'
dfchannel=df[df.channel==channel]
#hot_pix_map=hot_pix_maps_dict[channel]
# hot_pix_map=create_hot_pix_map_one_channel(dfchannel)
dfchannel=hot_pix_removal_one_channel(dfchannel, hot_pix_map_new)

#%%
df_clean=hot_pix_removal_several_channels(df_clean,hot_pix_maps_dict)
#%%
for idx, row in df_clean[df_clean.channel=='IR2'][1020:1030].iterrows():
    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    plot_CCDimage(row['ImageCalibrated'],fig=fig, axis=ax[0], title=row.channel)
    plot_CCDimage(row['ImageCalibrated_HPremoved'],fig=fig, axis=ax[1], title=row.channel)



# # plot the first 10 images
# for i in range(10):
#     fig, ax = plt.subplots(2, 1, figsize=(8, 8))
#     plot_CCDimage(dfchannel['ImageCalibrated'].iloc[i],fig=fig, axis=ax[0], title=f"{channel} - Image {i+1}")
#     plot_CCDimage(dfchannel['ImageCalibrated_HPremoved'].iloc[i],fig=fig, axis=ax[1], title=f"{channel} - Image {i+1} (HP Removed)")


#%%
# pickle df_clean
#df_clean.to_pickle(data_folder+'df_clean_'+name+'.pkl')
#if readincleandataframe:
#df_clean=pd.read_pickle(data_folder+'df_clean_'+name+'.pkl')


#%% 

#channel='IR1'
#df_channel = df[df['channel'] == channel].copy()
#select area TPlat between 20 and 70, TPlon between -30 and 30
if 'df' in locals():
    old_df = df.copy()
df = df_clean.copy()
print('Number of images before selection:', len(df))
df_sel = df[(df['TPlat'] >= 20) & (df['TPlat'] <= 70) & (df['TPlon'] >= -30) & (df['TPlon'] <= 30)].copy()

#Cut on time
# UTC-aware datetimes
t_start = DT.datetime(2025, 2, 19, 0, 0, 0, tzinfo=DT.timezone.utc)
t_stop  = DT.datetime(2025, 2, 20, 4, 59, 59, tzinfo=DT.timezone.utc)

df_sel = df_sel[(df_sel['EXPDate'] >= t_start) & (df_sel['EXPDate'] <= t_stop)]
print('Number of images after selection:', len(df_sel))


#%%
#Plot all channels
from mats_utils.plotting.plotCCD import all_channels_plot
import os

subfolder='SatSmokeLidarAllIncFlagsHPremoved'
imagedir=data_folder+subfolder+'/'
print('Images are to be found in '+imagedir)
# create folder if it does not exist
if not os.path.exists(imagedir):
    os.makedirs(imagedir)
all_channels_plot(df_sel,imagedir, field_of_choise='ImageCalibrated_HPremoved')
#all_channels_plot(df_sel,imagedir)




#cmap='YlGnBu_r'
# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2
channel='ALL'
if channel == 'IR1':
    generate_gif(imagedir+'CCDSEL1/', data_folder+'orbit_IR1_'+subfolder+'.gif')
elif channel == 'IR2':
    generate_gif(imagedir+'CCDSEL4/', data_folder+'orbit_IR2_'+subfolder+'.gif')
elif channel == 'IR3':
    generate_gif(imagedir+'CCDSEL2/', data_folder+'orbit_IR3_'+subfolder+'.gif')
elif channel == 'IR4':
    generate_gif(imagedir+'CCDSEL3/', data_folder+'orbit_IR4_'+subfolder+'.gif')
elif channel == 'NADIR':
    generate_gif(imagedir+'CCDSEL7/', data_folder+'orbit_NADIR_'+subfolder+'.gif')
elif channel == 'UV1':
    generate_gif(imagedir+'CCDSEL5/', data_folder+'orbit_UV1_'+subfolder+'.gif')
elif channel == 'UV2':
    generate_gif(imagedir+'CCDSEL6/', data_folder+'orbit_UV2_'+subfolder+'.gif')
elif channel=='ALL':
    generate_gif(imagedir+'ALL/', data_folder+'orbit_ALL_'+subfolder+'.gif')

print('Gif is to be found in '+data_folder+'orbit_'+channel+'_'+subfolder+'.gif')

# %%
plt.plot( df.TMHeaderTime, df.satheight)
plt.xlabel('Time')
plt.ylabel('Satellite Height')
plt.title('Satellite Height over Time')
plt.grid()
plt.show()

plt.plot(df.satlat, df.satheight)
plt.xlabel('Satellite Latitude')
plt.ylabel('Satellite Height')
plt.title('Satellite Height over Latitude')
plt.grid()
plt.show()

# %%
