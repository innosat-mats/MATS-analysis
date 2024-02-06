# Author: Linda Megner
# Created: October 26, 2023
# Script to investigate the flickering in UV2

#%%
from mats_utils.rawdata.read_data import read_MATS_data
from database_generation.experimental_utils import plot_CCDimage
from database_generation.flatfield import read_flatfield
from mats_l1_processing.instrument import Instrument
import datetime as DT
from mats_utils.rawdata.calibration import  calibrate_dataframe
from mats_l1_processing.L1_calibration_functions import calculate_flatfield, bin_image_with_BC, meanbin_image_with_BC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from lindas_own_functions import rename_CCDitem_entries
from database_generation import flatfield as flatfield
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items, CCDitems_to_dataframe
import pickle
from mats_utils.plotting.animate import generate_gif
import os


from mats_l1_processing.L1_calibration_functions import (
    get_true_image,
    desmear_true_image,
    subtract_dark,
    flatfield_calibration,
    get_linearized_image,
    combine_flags,
    make_binary,
    flip_image,
    handle_bad_columns,
    artifact_correction,
)

#%%

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'
starttime=DT.datetime(2023,5,7,1,0,0)
endtime=DT.datetime(2023,5,7,2,0,0)


filter_UV2={'CCDSEL': [6,6]} 
df=read_MATS_data(starttime, endtime,filter_UV2, level='1a',version='0.6')
print(len(df))
print(df.columns)


#pickle.dump(CCDitems, open('testdata/.pkl', 'wb'))

#%%
# Setup calibration
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)

def wrap_get_CCD(CCDitem, instrument):
    CCDunit=instrument.get_CCD(CCDitem["channel"])
    return CCDunit  

def wrap_get_true_image(CCDitem):
    image_bias_sub, error=get_true_image(CCDitem)
    return image_bias_sub
def wrap_get_linearized_image(CCDitem):    
    image_linear, error = get_linearized_image(CCDitem, CCDitem["image_bias_sub"], force_table=True)
    return image_linear
def wrap_desmear_true_image(CCDitem):
    image_desmeared, error=desmear_true_image(CCDitem, CCDitem["image_linear"])
    return image_desmeared
def wrap_subtract_dark(CCDitem):
    image_dark_sub, error=subtract_dark(CCDitem, CCDitem["image_desmeared"])
    return image_dark_sub
def wrap_flatfield_calibration(CCDitem):
    image_calib_nonflipped, error=flatfield_calibration(CCDitem, CCDitem["image_dark_sub"])
    return image_calib_nonflipped
def wrap_flip_image(CCDitem):
    image_calib_flipped=flip_image(CCDitem, CCDitem["image_calib_nonflipped"])
    return image_calib_flipped

CCDitems = dataframe_to_ccd_items(df[:50])

#%%
# Calibrate the images
def calibration_steps(CCDitem, instrument):
    CCDitem["CCDunit"] = wrap_get_CCD(CCDitem, instrument)
    CCDitem["image_bias_sub"] = wrap_get_true_image(CCDitem)
    CCDitem["image_linear"] = wrap_get_linearized_image(CCDitem)
    CCDitem["image_linear"] = np.flipud(CCDitem["image_linear"])
    CCDitem["image_desmeared"] = wrap_desmear_true_image(CCDitem)
    CCDitem["image_dark_sub"] = wrap_subtract_dark(CCDitem)
    CCDitem["image_calib_nonflipped"] = wrap_flatfield_calibration(CCDitem)
    CCDitem["ImageCalib"] = wrap_flip_image(CCDitem)
    return CCDitem



for CCDitem in CCDitems:
    CCDitem = calibration_steps(CCDitem, instrument) 

#pickle.dump(CCDitems, open('output/CCDitemsCalib.pkl', 'wb'))


# # Turn CCDitems back into a dataframe
#df = pd.DataFrame(CCDitems)
#reverse_rename_ccd_item_attributes(df)
#%%
df=CCDitems_to_dataframe(CCDitems)




# %%
