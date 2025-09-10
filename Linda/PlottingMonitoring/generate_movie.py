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

readindata=True
if readindata:

    start_time = DT.datetime(2025, 2, 19, 2, 0, 0)
    stop_time = DT.datetime(2025, 2, 20, 23, 59, 59)

    # filter
    #filter={'CCDSEL': [1,4]}
    #filter={'TPlat': [1,4]}

    # read in measurements

    df = read_MATS_data(start_time, stop_time,level='1b',version='1.0')
    #df = read_MATS_data(start_time, stop_time,filter,level='1b',version='1.0')
    name='df_all'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(data_folder+name+'.pkl')
else:
    satsmokename='df_satellite_smoke_1b_20250219_020000_20250220_235959.pkl'
    df = pd.read_pickle(data_folder+satsmokename)


#%% 
subfolder='movie28'
channel='IR2'
df_channel = df[df['channel'] == channel].copy()
imagedir=data_folder+subfolder+'/'+channel+'/'
orbit_plot(df_channel,imagedir,nbins=7, cmap='YlGnBu_r', plothistogram=False)   
print('Images are to be found in '+imagedir)
#cmap='YlGnBu_r'
# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2

if channel == 'IR1':
    generate_gif(imagedir+'CCDSEL1/', data_folder+'orbit_IR1_'+subfolder+'.gif')
elif channel == 'IR2':
    generate_gif(imagedir+'CCDSEL4/', data_folder+'orbit_IR2_'+subfolder+'.gif')
elif channel == 'UV2':
    generate_gif(imagedir+'CCDSEL6/', data_folder+'orbit_UV2_'+subfolder+'.gif')
print('Gif is to be found in '+data_folder+'orbit_'+channel+'_'+subfolder+'.gif')
# %%
