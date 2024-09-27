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
subdir='nonfiltered/'
directory = head_folder+subdir
directory = '/Users/lindamegner/MATS/MATS-retrieval/data/for_movie/'
os.makedirs(directory, exist_ok=True)
#directory = os.path.join(data_folder, 'satellite_smoke'+anomalymethod+'_'+channel+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'_'+crop)  # Create a directory to save the images in



# times for start and stop
starttime = DT.datetime(2023, 2, 10, 0, 0, 0)
endtime = DT.datetime(2023, 5, 10, 0, 0, 0)

# filter
#filter={'CCDSEL': [5,6]}

#%%
# read in measurements
#df = read_MATS_data(start_time, stop_time,filter,level='1a',version='0.6')

start_time = DT.datetime(2023, 2, 11, 12, 0, 0)
stop_time = DT.datetime(2023, 2, 11, 15, 0, 0)

tstarttime = start_time
while tstarttime < stop_time:
    tstoptime = tstarttime + DT.timedelta(hours=1)
    df = read_MATS_data(tstarttime, tstoptime,level='1b',version='0.6')
    name='df_'+tstarttime.strftime('%Y%m%d_%H%M%S')+'_'+tstoptime.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(directory+name+'.pkl')
    tstarttime = tstoptime
    print(name)
    
#%%
# df = pd.concat([df1, df2], ignore_index=True) 
df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
#save the data locally
#%%





#%%
#df.iloc


