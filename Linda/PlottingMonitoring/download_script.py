#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.rawdata.cropping import make_crop_filter
import os
from matplotlib import pyplot as plt
from mats_utils.plotting.plotCCD import plot_CCDimage


# head data folder
head_folder = '/Users/lindamegner/MATS/MATS-retrieval/data/rawdata/'
subdir='nonfiltered/'
directory = head_folder+subdir
directory = '/Users/lindamegner/MATS/MATS-retrieval/data/for_movie/'
directory='/Users/lindamegner/MATS/MATS-retrieval/data/daily_20230209-20230214/'
os.makedirs(directory, exist_ok=True)
#directory = os.path.join(data_folder, 'satellite_smoke'+anomalymethod+'_'+channel+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'_'+crop)  # Create a directory to save the images in





# filter
#filter={'CCDSEL': [5,6]}

#%%
# read in measurements
#df = read_MATS_data(start_time, stop_time,filter,level='1a',version='0.6')
# times for start and stop
start_time = DT.datetime(2023, 1, 11, 0, 0, 0)
stop_time = DT.datetime(2023, 1, 11, 6, 0, 0)
headdirectory = '/Users/lindamegner/MATS/MATS-retrieval/data/'
version='1.0'
level='1a'
directory = os.path.join(headdirectory, level+'_'+version+'_'+start_time.strftime('%Y%m%d')+'-'+stop_time.strftime('%Y%m%d')+'/')
os.makedirs(directory, exist_ok=True)


tstarttime = start_time
while tstarttime < stop_time:
    tstoptime = tstarttime + DT.timedelta(hours=1)
    # try to read in data for this hour
    try:
        df = read_MATS_data(tstarttime, tstoptime,level=level,version=version)
        name='df_'+tstarttime.strftime('%Y%m%d_%H%M%S')+'_'+tstoptime.strftime('%Y%m%d_%H%M%S')
        df.to_pickle(directory+name+'.pkl')
        tstarttime = tstoptime
        print(name)
    except Exception as e:
        print(f"Error reading data for {tstarttime} to {tstoptime}: {e}")
        tstarttime = tstoptime

#%%
# df = pd.concat([df1, df2], ignore_index=True) 
df = read_MATS_data(start_time, stop_time,level=level,version=version)
#save the data locally
#%%
# read in measurements
#df= pd.read_pickle(directory+'df_'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')+'.pkl') 


if level=='1a':
    plotfield='IMAGE'
elif level=='1b':
    plotfield='ImageCalibrated'

for index, CCD in df[:10].iterrows():
    plot_CCDimage(CCD[plotfield])


#%%
#plot mean of all images in df
image_mean = df[plotfield].mean()
plot_CCDimage(image_mean)



#%%
#df.iloc


