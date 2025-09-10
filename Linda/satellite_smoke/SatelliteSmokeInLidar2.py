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

readinrawdata=False

event_start_time = DT.datetime(2025, 2, 19, 0, 0, 0)
event_stop_time = DT.datetime(2025, 2, 20, 23, 59, 59)
deltatime = event_stop_time - event_start_time
level='1b'
version='1.0'


# set start times every 10 days prior to the event so that each time slot covers the same time as the event duration
starttimes = [event_start_time-DT.timedelta(days=30),
            event_start_time-DT.timedelta(days=20),
            event_start_time-DT.timedelta(days=10),
            event_start_time]
endtimes = [starttime + deltatime for starttime in starttimes]

name='df_SatSmokeInLidar2'+event_start_time.strftime('%Y%m%d_%H%M%S')+'_'+event_stop_time.strftime('%Y%m%d_%H%M%S')+'_'+level+'v'+version
# filter
filter={'satlat': [20,70], 'satlon': [-20,30]} 

if readinrawdata:
    df_list = []
    for start_time, stop_time in zip(starttimes, endtimes):
        df = read_MATS_data(start_time, stop_time,filter=filter,level=level,version=version)
        df_list.append(df)

    pickle.dump(df_list, open(data_folder+'df_list_'+name+'.pkl', 'wb'))
else:
    df_list = pickle.load(open(data_folder+'df_list_'+name+'.pkl', 'rb'))   


#%%
#concatenate all dataframes
dflong=pd.concat(df_list, axis=0)
#%%
createhotpixmap=True

if createhotpixmap:
    hot_pix_maps_dict=create_all_hot_pix_maps(dflong, level=level)
    with open(data_folder+name+'hot_pixel_maps.pkl', 'wb') as f:
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
#Fix the IR2 hot pix map, 
hot_pix_maps_dict_old = pickle.load(open(data_folder+'hot_pixel_maps.pkl', 'rb'))
def apply_running_mean_filter(img, filter_rows, filter_cols):
    from scipy.ndimage import uniform_filter
    av_dark_signal= np.nanmean(img[int(img.shape[0]*0.2):int(img.shape[0]*0.8)])
    running_mean = uniform_filter(img, size=(filter_rows, filter_cols))
    image = img.copy().astype(float)
    img_fixed = image - running_mean 
    return img_fixed

hot_pix_map_new=apply_running_mean_filter(hot_pix_maps_dict_old['IR2'], filter_rows=2, filter_cols=5)
plot_CCDimage(hot_pix_map_new, title='Fixed Hot Pixel Map - IR2')
hot_pix_maps_dict['IR2'] = hot_pix_map_new



#%%

for df in df_list:
    df = hot_pix_removal_several_channels(df, hot_pix_maps_dict)
#%%
df = df_list[-1]
for idx, row in df[df.channel=='IR1'][120:130].iterrows():
    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    plot_CCDimage(row['ImageCalibrated'],fig=fig, axis=ax[0], title=row.channel)
    plot_CCDimage(row['ImageCalibrated_HPremoved'],fig=fig, axis=ax[1], title=row.channel)

#%%
# pickle df_list
pickle.dump(df_list, open(data_folder+'df_list_hpremoved_'+name+'.pkl', 'wb'))
#%%

df_list = pickle.load(open(data_folder+'df_list_hpremoved_'+name+'.pkl', 'rb'))


#%% 
#Now check how many pixels are much higher than the running mean for the different dataframes df

def apply_running_mean_filter(img, filter_rows, filter_cols):
    from scipy.ndimage import uniform_filter
    running_mean = uniform_filter(img, size=(filter_rows, filter_cols))
    image = img.copy().astype(float)
    img_fixed = image - running_mean
    return img_fixed

def count_pixels_above_threshold(df, threshold):
    """Counts the number of hot pixels in the DataFrame."""
    nr_hot_pixels = df['ImageCalibrated_RMfiltered'].apply(lambda img: np.sum(img > threshold))
    return nr_hot_pixels


threshold = 50
for df in df_list:
    df['ImageCalibrated_RMfiltered'] = df['ImageCalibrated'].apply(lambda img: apply_running_mean_filter(img, filter_rows=2, filter_cols=5))
    df['nr_strong_pixels'] = count_pixels_above_threshold(df, threshold=threshold)
    df['mean_strong_pixel_value'] = df['ImageCalibrated_RMfiltered'].apply(lambda img: np.nanmean(img[img > threshold]))


#%%
# plot histograms
channel='IR1'
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
for df in df_list:

    time = df.iloc[0].TMHeaderTime
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    df['nr_strong_pixels'].hist(bins=50, ax=ax[0])
    ax[0].set_title('Histogram of Number of Strong Pixels (>50)')
    ax[0].set_xlabel('Number of Strong Pixels')
    ax[0].set_ylabel('Frequency')
    ax[0].text(0.5, 0.9, f"Mean nr of strong pixels: {df['nr_strong_pixels'].mean():.2f}", color='r', ha='center', transform=ax[0].transAxes)   



    df['mean_strong_pixel_value'].hist(bins=50, ax=ax[1])
    ax[1].set_title('Histogram of Mean Strong Pixel Value')
    ax[1].set_xlabel('Mean Strong Pixel Value')
    ax[1].set_ylabel('Frequency')
    ax[1].text(0.5, 0.9, f"Mean: {df['mean_strong_pixel_value'].mean():.2f}", color='r', ha='center', transform=ax[1].transAxes)

    fig.suptitle(f"Time: {time}, Channel: {channel}")

    plt.tight_layout()
    plt.show()

# %%
