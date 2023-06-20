#%%
#from types import NoneType
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate,calibrate_all_items
import datetime as DT
from mats_l1_processing.instrument import Instrument
import numpy as np
import matplotlib.pyplot as plt
import time
import os


def calibrate(CCDitem, instrument):
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated_flipped,
        image_calibrated,
        errors
    ) = L1_calibrate(CCDitem, instrument)
    return image_calibrated
#%%
os.chdir('/home/louis/MATS')
calibration_file='/home/louis/MATS/tests/calibration_data_test.toml'
instrument=Instrument(calibration_file)



# start_time=DT.datetime(2023,2,13,9,5,0)
# stop_time=DT.datetime(2023,2,13,9,15,0)

start_time=DT.datetime(2023,4,7,6,0,0)
stop_time=DT.datetime(2023,4,7,7,0,0)

df=read_MATS_data(start_time, stop_time,level='1a',version='0.6')


CCDitems=df

CCDitems['NCBIN CCDColumns'] = CCDitems['NCBINCCDColumns'] 
CCDitems['GAIN Truncation'] = CCDitems['GAINTruncation'] 
CCDitems['NCBIN FPGAColumns'] = CCDitems['NCBINFPGAColumns']

mylist=[]
CCDitems['BC'] = [mylist for i in CCDitems.index]
CCDitems['GAIN Mode'] = CCDitems['GAINMode'] 


#%%

for index,CCDitem in CCDitems.iterrows():
    print(index)
    im_cal = calibrate(CCDitem,instrument)


# calibrate_all_items(CCDitems,instrument)





# image_calibrated = calibrate(CCDitem_nadir,instrument)
# %%
