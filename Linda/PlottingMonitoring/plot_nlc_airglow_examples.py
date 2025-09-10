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
    start_time = DT.datetime(2024, 6, 25, 10, 47, 30)
    stop_time = DT.datetime(2024, 6, 25, 10, 49, 30)
    df = read_MATS_data(start_time, stop_time,level=level,version='0.6')
    data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/PeterDalin4/'
if stringn != 'peter':
    df = read_MATS_data(start_time, stop_time,filter,level=level,version='1.0')
name='df_'+stringn+'_'+level+'_'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
#%%

#df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
if stringn in ['airglow', 'dayglow','nightglow']:
    channel='IR2'
    colorscheme='viridis'
elif stringn == 'nlc': 
    channel='UV2'
    colorscheme='Blues_r'
elif stringn == 'peter':
    channel='IR1'
    colorscheme='inferno'
else:
    Exception('Unknown stringn: {}'.format(stringn))

df_channel = df[df['channel'] == channel].copy()
df_channel = df.copy()

#%%

#plot every image in df_channel
for i in range(0, len(df_channel)):


    fig, ax, img = plot_image(df_channel.iloc[i], cmap=colorscheme,  save=False)
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
    ax.set_title('{} image at {}, TPlat: {}, TP lon{}'.format(df_channel.iloc[i]['channel'], df_channel.iloc[i]['EXPDate'].strftime('%Y-%m-%d %H:%M:%S'), df_channel.iloc[i]['TPlat'], df_channel.iloc[i]['TPlon']))
    #save figure
    fig.savefig(data_folder + name + '_image_' + str(i) + '.png', dpi=300, bbox_inches='tight')

#  %%
