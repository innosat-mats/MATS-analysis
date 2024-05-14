import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import argparse
from matplotlib import pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
import pdb

# MATS NLC Profil
'''
The MATS image should be reduced from 2D to 1D by integrating over an appropriate number of columns 
(e.g. corresponding the the 40 km across-track-FOV of OSIRIS).  
'''
'''
starttime=datetime(2023,2,11,0,0)
stoptime=datetime(2023,2,12,0,0)
df=read_MATS_data(starttime,stoptime, level="1b", version="0.6", filter={"CCDSEL":[5,5]}) # df.shape[0] --> number of rows; df.shape[1] --> number of columns
# CCD sensor number (1-7). Corresponds to channel name (see 'channel' variable) as: (1: 'IR1'), (2: 'IR4'), (3: 'IR3'), (4: 'IR2'), (5: 'UV1'), (6: 'UV2') and  (7: 'NADIR')

# save dataframe (pickle)
df.to_pickle('NLC_test_file.pkl')
'''
'''
# load dataframe
df = pd.read_pickle('NLC_test_file.pkl')
UV2_df = df[df['channel'] == 'UV2']
UV2 = df[df['channel'] == 'UV2'].dropna().reset_index(drop=True)

# Technical data
columns = np.arange(0, UV2["NCOL"][0], 1)
ref_col = columns[int(len(columns) / 2)] # reference column
rows = np.arange(0, UV2["NROW"][0] - 10, 1)
num_images = len(UV2)

# plot iamge
plot_image(df.iloc[3], save=False)
plt.show(block=True)

plt.pcolormesh(df['ImageCalibrated'].iloc[3])
plt.show()

plt.pcolormesh(UV2['ImageCalibrated'].iloc[3])
plt.show()
'''
# Code Björn
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

starttime=datetime(2023,2,11,0,0)
stoptime=datetime(2023,2,12,0,0)

channel='UV1'
l1b_version="0.6"

# some filter settings - dayglow (ALL SZA)
dmin,dmax = 0, 95
tplat0, tplat1 = -90, 90
dftop=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1]})
df = dftop[dftop['channel'] == channel].dropna().reset_index() # not sure if necessary

# Call and plot
#for i in range(len(df)):
for i in range(0,3):
    heights, profile = prepare_profile(df.iloc[i])
    plt.plot(profile,heights)
    plt.show()


# pdb.set_trace()

