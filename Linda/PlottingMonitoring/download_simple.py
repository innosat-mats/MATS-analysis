#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.rawdata.cropping import make_crop_filter
import os


# head data folder
head_folder = '/Users/lindamegner/MATS/MATS-retrieval/data/rawdata/'
subdir='level1a_v1.0/'
directory = head_folder+subdir
#directory = '/Users/lindamegner/MATS/MATS-retrieval/data/for_movie/'
os.makedirs(directory, exist_ok=True)
#directory = os.path.join(data_folder, 'satellite_smoke'+anomalymethod+'_'+channel+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'_'+crop)  # Create a directory to save the images in




#%%
# read in measurements
#df = read_MATS_data(start_time, stop_time,filter,level='1a',version='0.6')

start_time = DT.datetime(2026, 7, 6, 22, 0, 0)
stop_time = DT.datetime(2026, 7, 9, 3, 0, 0)

tstarttime = start_time
while tstarttime < stop_time:
    tstoptime = tstarttime + DT.timedelta(hours=2)
    df = read_MATS_data(tstarttime, tstoptime,level='1a',version='1.0')
    name='df_'+tstarttime.strftime('%Y%m%d_%H%M%S')+'_'+tstoptime.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(directory+name+'.pkl')
    tstarttime = tstarttime + DT.timedelta(days=1)
    print(name)
    
#%%
# df = pd.concat([df1, df2], ignore_index=True) 
df = read_MATS_data(start_time, stop_time,level='1a',version='1.0')
#save the data locally
#%%





#%%
#df.iloc


