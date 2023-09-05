# Code to test the hot pixel desmearing algorithm, as we expect it to overcorrect the hot pixels.


#%%
from mats_utils.rawdata.read_data import read_MATS_data
from database_generation.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
from mats_utils.rawdata.calibration import  calibrate_dataframe
import datetime as DT
import copy
import matplotlib.pyplot as plt
import numpy as np
from mats_utils.statistiscs.images_functions import create_imagecube
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot

#%%
# image folder
image_folder = 'images/'

 
starttime=DT.datetime(2023,2,23,0,0,0)
endtime=DT.datetime(2023,2,24,0,0,0)

#Fliter on full frame to get the star images.
filter_fullframe={'NRBIN': [1, 1], 'NCBINCCDColumns': [1, 1], 'NCOL':[2047,2048], 'NROW':[511,512]} 
filter_nightglow={'NRBIN': [2, 2], 'TPsza': [98.5, 150], 'CCDSEL':[1,1]}  

df=read_MATS_data(starttime, endtime,filter_fullframe, level='1a',version='0.5')


print(len(df))



#%%
dfl1a=df
for index, CCDitem in dfl1a.iterrows():
    plot_image(CCDitem, save=False)


#%%

#Calibrate
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
#calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/calibration_data.toml'
instrument=Instrument(calibration_file)

#dfl1b=calibrate_dataframe(dfl1a, instrument)
dfl1b=calibrate_dataframe(dfl1a, instrument, debug_outputs=True)
#%%
xstart=0
xstop=2027 #200
ystart=0 #100
ystop=512 #300

for i, CCDitem in dfl1b.iterrows():
    fig, ax=plt.subplots(3,1)
    plot_CCDimage( CCDitem['image_desmeared'][xstart:xstop, ystart:ystop], fig, ax[2], title='after desmearing')
    plot_CCDimage( CCDitem['image_bias_sub'][xstart:xstop, ystart:ystop], fig, ax[0], title='prior to desmearing')
    plot_CCDimage( (CCDitem['image_bias_sub']-CCDitem['image_desmeared'])[xstart:xstop, ystart:ystop], fig, ax[1], title='desmearing correction', clim=[0, 0.5])



# %%
