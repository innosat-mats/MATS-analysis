#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot

# data folder
data_folder = '/home/louis/MATS/MATS-Data/'

# times for start and stop
start_time = DT.datetime(2023, 2, 12, 6, 0, 0)
stop_time = DT.datetime(2023, 2, 12, 6, 0, 30)

# filter
# filter={'CCDSEL': [7,7]}

#%%
# read in measurements
df = read_MATS_data(start_time, stop_time)

#%%
#df.iloc

for index, CCD in df.iterrows():
    plot_image(CCD, outpath=data_folder)
    

simple_plot(df,data_folder)

#%% testing orbit_plot changes

orbit_plot(df,data_folder+'test1/',nbins=7.9)
orbit_plot(df,data_folder+'test2/',nbins=3)
orbit_plot(df,data_folder+'test3/',nbins=100)
orbit_plot(df,data_folder+'test4/',nbins=None)
orbit_plot(df,data_folder+'test5/')
# %%
