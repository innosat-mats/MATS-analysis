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

start_time = DT.datetime(2023, 2, 22, 4, 50, 0)
stop_time = DT.datetime(2023, 2, 22, 4, 51, 0)
#start_time = DT.datetime(2025, 1, 24, 0, 0, 0)
#stop_time = DT.datetime(2025, 1, 25, 0, 0, 0)

# filter
CCDfilter={'CCDSEL': [5, 6]}

#%%
# read in measurements

df102 = read_MATS_data(start_time, stop_time,level='1b',filter=CCDfilter,version='1.0.2')
# %%
