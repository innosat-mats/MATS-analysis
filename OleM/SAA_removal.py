#%%
import numpy as np
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from matplotlib import pyplot as plt

#%%
starttime = DT.datetime(2023, 3, 1, 14, 0)   
endtime = DT.datetime(2023, 3, 1, 20, 0)
ccd_data_cal = read_MATS_data(starttime, endtime,filter=None,level='1b',version='0.4')


#%%
ccd_data_cal

#%%
# data = ccd_data_cal.ImageCalibrated[ccd_data_cal.CCDSEL==1]