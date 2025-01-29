#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif


# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

# times for start and stop
start_time = DT.datetime(2023, 2, 4, 9, 0, 0)
stop_time = DT.datetime(2023, 2, 4, 10, 0, 30)

# filter
filter={'channel': 'IR1', 'TPsza': [90,200]}

#%%
# read in measurements
#df = read_MATS_data(start_time, stop_time,filter,level='1a',version='0.6')
df = read_MATS_data(start_time, stop_time,filter, level='1b',version='0.9')

#%%
#plot TPsza
from matplotlib import pyplot as plt
from mats_utils.plotting.plotCCD import plot_CCDimage
df['TPsza'].plot()

plot_CCDimage(df['ImageCalibrated'][0]) 

#plot mean of all images in df
image_mean = df['ImageCalibrated'].mean(1)
plot_CCDimage(image_mean)

#%%
#df.iloc

for index, CCD in df[:4].iterrows():
    plot_image(CCD, outpath=data_folder)
    

simple_plot(df,data_folder)

#%% testing orbit_plot changes

orbit_plot(df,data_folder+'test2/',nbins=7)

# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2
generate_gif('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2/CCDSEL5/', data_folder+'orbit.gif')
# %%
