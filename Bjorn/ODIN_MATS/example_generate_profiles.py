#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
import matplotlib.pyplot as plt
#import matplotlib.pylab as plt

#%%
def prepare_profile(ch):
    # This function averages some columns and
    # calculates tangent heights
    
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    # multiply with factor to get right values (depends on version?)
    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e13)

    # set heights
    common_heights = np.arange(60000,110250,250)
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile

#%%
starttime=datetime(2023,3,1,0,0)
stoptime=datetime(2023,3,1,10,0)

channel='IR1'
l1b_version="0.6"

# some filter settings - dayglow (ALL SZA)
dmin,dmax = 0, 95
tplat0, tplat1 = -90, 90
dftop=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1]})
df = dftop[dftop['channel'] == channel].dropna().reset_index() # not sure if necessary
#%%
# Call and plot
#for i in range(len(df)):
for i in range(0,3):
    heights, profile = prepare_profile(df.iloc[i])
    plt.plot(profile,heights)
# %%
