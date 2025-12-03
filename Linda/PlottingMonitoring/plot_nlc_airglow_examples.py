#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
import numpy as np

from database_generation.experimental_utils import plot_CCDimage


# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'


stringn='peter'
level='1b'  # '1a' or '1b'
if stringn=='nlc':
    start_time = DT.datetime(2023, 2, 2, 19, 41, 0)
    stop_time = DT.datetime(2023, 2, 2, 19, 43, 0)
    filter={'CCDSEL': [5,6]}
elif stringn=='dayglow':
    start_time = DT.datetime(2023, 2, 11, 12, 38, 0)
    stop_time = DT.datetime(2023, 2, 11, 13, 28, 0)
    filter={'CCDSEL': [1,4]}
elif stringn=='nightglow':
    start_time = DT.datetime(2023, 2, 11, 13, 48, 0)
    stop_time = DT.datetime(2023, 2, 11, 14, 1, 0)
    filter={'CCDSEL': [1,4]}
elif stringn=='peter':
    #start_time = DT.datetime(2024, 6, 25, 10, 47, 30)
    #stop_time = DT.datetime(2024, 6, 25, 10, 49, 30)
    start_time = DT.datetime(2024, 6, 25, 10, 47, 50)
    stop_time = DT.datetime(2024, 6, 25, 10, 48, 20)
    df = read_MATS_data(start_time, stop_time,level=level,version='0.6')
    data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/PeterDalin5/'
if stringn != 'peter':
    df = read_MATS_data(start_time, stop_time,filter,level=level,version='1.0')
name='df_'+stringn+'_'+level+'_'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
#%%

from matplotlib.colors import LinearSegmentedColormap
# Define a custom colormap from dark blue → light blue → white
colors = [
    # (0.0, '#001f66'),  # Deep blue
    # (0.3, '#004c99'),  # Medium blue
    # (0.7, '#3399ff'),  # Light blue
    # (1.0, '#ffffff')   # White
    (0.0, '#001f66'),  # Deep blue
    (0.4, '#004c99'),  # Medium blue
    (0.7, '#00bfff'),  # Bright sky blue
    (0.8, '#00ffff'),   # Cyan 
    (1.0, '#ffffff')   # White

]

#df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
if stringn in ['airglow', 'dayglow','nightglow']:
    channel='IR2'
    colorscheme='viridis'
elif stringn == 'nlc': 
    channel='UV2'
    colorscheme= LinearSegmentedColormap.from_list("blue_to_white", colors)
elif stringn == 'peter':
    channel='UV1'
    colorscheme= LinearSegmentedColormap.from_list("blue_to_white", colors)

else:
    Exception('Unknown stringn: {}'.format(stringn))

df_channel_hot = df[df['channel'] == channel].copy()
#df_channel = df.copy()


#%%
hotpixelcorr=True
if hotpixelcorr:
    #Hot pixel removal
    from mats_utils.retrieval.hot_pix import create_hot_pix_map_one_channel, create_all_hot_pix_maps, hot_pix_removal_several_channels, hot_pix_removal_one_channel
    start_time_hotmap = start_time - DT.timedelta(hours=2)
    stop_time_hotmap = start_time - DT.timedelta(hours=0.5)
    filter = {'CCDSEL': [5, 5]} # 5 is UV1
    dflong_channel = read_MATS_data(start_time_hotmap, stop_time_hotmap,level=level,version='0.6', filter=filter)

    hot_pix_map=create_hot_pix_map_one_channel(dflong_channel, thresholdalt=100000, remove_background=False)
    plot_CCDimage(hot_pix_map)
    df_channel=hot_pix_removal_one_channel(df_channel_hot, hot_pix_map)
    image_field='ImageCalibrated_HPremoved'
else:
    image_field='ImageCalibrated'
    df_channel=df_channel_hot.copy()

#%%


#plot every image in df_channel
for i in range(0, len(df_channel)):

    fig, ax, img = plot_image(df_channel.iloc[i], cmap=colorscheme,  save=False, image_field=image_field, ranges=[0 ,400])
    #set xlabel and y label
    ax.set_xlabel('CCD Column')
    ax.set_ylabel('CCD Row')
    #set figure size
    fig.set_size_inches(10, 3.5)
    #add colorbar
    cbar = fig.colorbar(img, ax=ax, orientation='vertical')
    # set lable for colorbar to 'Brightness in 10^12 photons/m^2/s/nm/sr'
    cbar.set_label('Brightness in 10^12 photons/m^2/s/nm/sr', rotation=270, labelpad=20)
    #add text to current title
    ax.set_title('{} image at {}, TPlat: {:.2f}, TP lon: {:.2f}'.format(
        df_channel.iloc[i]['channel'],
        df_channel.iloc[i]['EXPDate'].strftime('%Y-%m-%d %H:%M:%S'),
        df_channel.iloc[i]['TPlat'],
        df_channel.iloc[i]['TPlon']
    ))
    #save figure
    #fig.savefig(data_folder + name + '_image_' + str(i) + '.png', dpi=300, bbox_inches='tight')

#  %%
