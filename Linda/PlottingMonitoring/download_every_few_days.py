#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.rawdata.cropping import make_crop_filter
import os
import pickle


# head data folder

directory = '/Users/lindamegner/MATS/MATS-retrieval/data/every_few_days/'
os.makedirs(directory, exist_ok=True)



#%%


start_time = DT.datetime(2023, 10, 1, 1, 0, 0)
stop_time = DT.datetime(2024, 10, 1, 1, 0, 0)
deltatime = DT.timedelta(days=3.1)

tstarttime = start_time
while tstarttime < stop_time:
    tstoptime = tstarttime + DT.timedelta(minutes=96)
    try:
        df= read_MATS_data(tstarttime, tstoptime, level='1a', version='0.7', pritfilesys=False)
    except:
        print('found no data for intervalstarttime', tstarttime)
    else: #if the dataframe is not empty, ie MATS data was found
        name='df_'+tstarttime.strftime('%Y%m%d_%H%M%S')+'_'+tstoptime.strftime('%Y%m%d_%H%M%S')
        df.to_pickle(directory+name+'.pkl')
        print('pickling:'+name)
    tstarttime = tstoptime+deltatime






#%%




# %%
