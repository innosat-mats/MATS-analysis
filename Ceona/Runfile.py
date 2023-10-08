#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from skyfield.api import load
from skyfield.units import Distance
from skyfield.toposlib import wgs84
from skyfield.positionlib import Geocentric
from aacgmv2 import get_aacgm_coord

start_time = DT.datetime(2023,2,15,00,0,0)
stop_time = DT.datetime(2023,2,17,00,0,0)
channel = 'IR1'
numdays = stop_time-start_time #number of days
#%%
dfSH = read_MATS_data(start_time,stop_time,filter={"TPlat":[-90,-45],'NROW': [0,400]},version='0.5',level='1b')
dfSH = dfSH[dfSH['channel'] == channel]
#dfSH.to_pickle('22to28febIR1SH')
#%%
dfNH = read_MATS_data(start_time,stop_time,filter={"TPlat":[45,90],'NROW': [0,400]},version='0.5',level='1b')
dfNH = dfNH[dfNH['channel'] == channel]
#dfNH.to_pickle('22to28febIR1NH')

#%%
items = pd.concat([dfNH,dfSH], ignore_index=True)
sortitems = items.sort_values(['EXPDate'])
sortitems.to_pickle('15to16febIR1')

# %%
ccd = pd.read_pickle('15to16febIR1')

#%% Plot a chosen image
def saveSpecIm(ccditem):
    ccdimage = ccditem['ImageCalibrated']
    #saves the matrix to matlab-file
    scipy.io.savemat('aurorapic2',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
    return

# %%
image = ccd.iloc[400].ImageCalibrated
plt.pcolormesh(image,rasterized = True, vmin=0, vmax=480)

# %%
