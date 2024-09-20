#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from database_generation.experimental_utils import plot_CCDimage
import matplotlib.pyplot as plt
import pickle

def plot_correlations(df, time_of_day='', color='black'):




    #Define the center of the image
    midrow=round(df['IMAGE'][0].shape[0]/2) #94
    midcol=round(df['IMAGE'][0].shape[1]/2) #22
    maxrow = df['IMAGE'][0].shape[0]
    maxcol = df['IMAGE'][0].shape[1]

    df['image_center_mean'] = df['IMAGE'].apply(lambda x: x[midrow-20:midrow+20, midcol-5:midcol+5].mean())
    df['image_corner_mean'] = df['IMAGE'].apply(lambda x: x[maxrow-10:maxrow, maxcol-3:maxcol].mean())
    df['image_upper_mean'] = df['IMAGE'].apply(lambda x: x[maxrow-10:maxrow, midcol-5:midcol+5].mean())

    # Plot histograms
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    ax[0].hist(df['image_center_mean'], bins=50, color=color)
    ax[0].set_xlabel('Signal in center of image')
    ax[0].set_ylabel('Frequency')
    #ax[0].set_title('Histogram of signal in center of image')

    ax[1].hist(df['image_corner_mean'], bins=50, color=color)
    ax[1].set_xlabel('Signal in corner of image')
    ax[1].set_ylabel('Frequency')
    #ax[1].set_title('Histogram of signal in corner of image')

    ax[2].hist(df['image_upper_mean'], bins=50, color=color)
    ax[2].set_xlabel('Signal in upper part of image')
    ax[2].set_ylabel('Frequency')
    #ax[2].set_title('Histogram of signal in upper part of image')
    fig.suptitle(time_of_day)


    #make figure with subplots
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    ax[0].scatter(df['image_corner_mean'], df['image_center_mean'],  label=time_of_day, color=color)
    ax[0].set_xlabel('Signal in corner of image')
    ax[0].set_ylabel('Signal in center of image')
    #ax[0].set_title('Correlation between corner and center of image')

    ax[1].scatter(df['image_corner_mean'], df['image_upper_mean'], label=time_of_day, color=color)
    ax[1].set_xlabel('Signal in corner of image')
    ax[1].set_ylabel('Signal in upper part of image')
    #ax[1].set_title('Correlation between corner and upper part of image')

    ax[2].scatter(df['image_upper_mean'], df['image_center_mean'], label=time_of_day, color=color)
    ax[2].set_xlabel('Signal in upper part of image')
    ax[2].set_ylabel('Signal in center of image')
    #ax[2].set_title('Correlation between upper part and center of image')
    fig.suptitle(time_of_day)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # plot the correlation with solar zenith angle
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    ax[0].scatter(df['TPsza'], df['image_center_mean'],  label=time_of_day, color=color)
    ax[0].set_xlabel('Solar zenith angle')
    ax[0].set_ylabel('Signal in center of image')
    #ax[0].set_title('Correlation between corner and center of image')

    ax[1].scatter(df['TPsza'], df['image_corner_mean'], label=time_of_day, color=color)
    ax[1].set_xlabel('Solar zenith angle')
    ax[1].set_ylabel('Signal in corner of image')
    #ax[1].set_title('Correlation between corner and upper part of image')

    ax[2].scatter(df['TPsza'], df['image_upper_mean'], label=time_of_day, color=color)
    ax[2].set_xlabel('Solar zenith angle')
    ax[2].set_ylabel('Signal in upper part of image')
    #ax[2].set_title('Correlation between upper part and center of image')
    fig.suptitle(time_of_day)
    


#%%
output_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'
datafolder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-data/'
name=straylight_correlation
# times for start and stop
start_time = DT.datetime(2023, 1, 4, 0, 0, 0)
stop_time = DT.datetime(2023, 1, 8, 0, 0, 0)

# filter
filter={'CCDSEL': [1,1], 'schedule_name': 'MODE1y'}

#%%
# read in measurements
dfo = read_MATS_data(start_time, stop_time,filter=filter, level='1a',version='0.7')
print('length of df', len(dfo))

#%%
#%%
pickle.dump(dfo, open(datafolder+'df_'+start_time.strftime('%Y%m%d')+'_'+stop_time.strftime('%Y%m%d')+'.pkl', 'wb'))

#%%
# Separate data into day and night
df_day = dfo[(dfo['TPsza'] <= 90)]
#re-index the dataframe
df_day.reset_index(drop=True, inplace=True)
df_night = dfo[(dfo['TPsza'] >= 100)]
#re-index the dataframe
df_night.reset_index(drop=True, inplace=True)


# %%
if len(df_day) > 20:
    plot_correlations(df_day, time_of_day='Day', color='orange')
#%%    
if len(df_night) > 20:
    plot_correlations(df_night, time_of_day='Night', color='blue')

# %%
