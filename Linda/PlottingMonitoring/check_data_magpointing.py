#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from matplotlib import pyplot as plt
from database_generation.experimental_utils import plot_CCDimage



#Check available data with 
#aws s3 ls ops-payload-level1b-v0.6/2023/5/ --profile mats 

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

# times for start and stop

#start_time = DT.datetime(2023, 2, 12, 4, 50, 0)
#stop_time = DT.datetime(2023, 2, 12, 4, 54, 0)
start_time = DT.datetime(2024, 6, 25, 0, 0, 0)
stop_time = DT.datetime(2024, 6, 26, 0, 0, 0)

# filter
#filter={'channel': ['IR1', 'UV1', 'UV2']} 63 79
filter={'CCDSEL': [1,6], 'TPlat': [63, 79], 'TPlon': [290, 340]}
#%%
# read in measurements
df = read_MATS_data(start_time, stop_time,level='1b',version='0.6',filter=filter)



#%%
#df.iloc


for index, CCD in df[:6].iterrows():
    plot_CCDimage(df.iloc[index].ImageCalibrated,  title=df.iloc[index].channel)

#%%
df=df[:6]
simple_plot(df,data_folder)

#%% testing orbit_plot changes
picdir='PeterDalin'
orbit_plot(df,data_folder+picdir+'/',nbins=7)

# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2
for i in range(1,7):

    generate_gif('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'+picdir+'/CCDSEL'+str(i)+'/', data_folder+picdir+str(i)+'.gif')
# %%
