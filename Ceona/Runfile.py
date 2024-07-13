#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from aurora_analysis import  get_aurora_max
from Auroradetection import get_strips
from Keogram import get_stripRow
from orbitoverview import overview_points
from gradientOverview import overview_grad

start_time = DT.datetime(2023,3,15,0,0,0)
stop_time = DT.datetime(2023,3,23,0,0,0)
channel = 'IR1'
numdays = stop_time-start_time
#%% Read in data for Southern hemisphere, using certain filter
dfSH = read_MATS_data(start_time,stop_time,filter={"TPlat":[-90,-45],'NROW': [0,400]},version='0.5',level='1b')
dfSH = dfSH[dfSH['channel'] == channel]

#%% Read in data for Northern hemisphere, using certain filter
dfNH = read_MATS_data(start_time,stop_time,filter={"TPlat":[45,90],'NROW': [0,400]},version='0.5',level='1b')
dfNH = dfNH[dfNH['channel'] == channel]

#%% Concatenate the two hemispheres, and sort the CCD on time, saves in pickle
items = pd.concat([dfNH,dfSH], ignore_index=True)
sortitems = items.sort_values(['EXPDate'])
sortitems.to_pickle('MatsData/mar3W')

#%% Saves a spcific image into matlab matrix
def saveSpecIm(ccditem):
    ccdimage = ccditem['ImageCalibrated']
    #saves the matrix to matlab-file
    plt.pcolormesh(ccdimage,rasterized = True) # vmin=0, vmax=180
    #scipy.io.savemat('feb15_1901IR4.mat',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
    return

# %% Main file to run scripts
def Main():
    start_time = DT.datetime(2023,3,15,0,0,0)
    stop_time = DT.datetime(2023,3,23,0,0,0)
    numdays = stop_time-start_time #number of days
    filedate = 'mar3W'
    items = pd.read_pickle(r'MatsData\mar3W')

    #returns aurorastrips and saves the other strips in pickles.
    aurorastrips = get_strips(items,numdays,filedate)  
    get_aurora_max(aurorastrips,filedate)

    #load in all the strips from file to plot red dots in overviews
    allstrips = pd.read_pickle(r'MatsData\mar3Wallstrips')
    allrows = get_stripRow(allstrips)

    #Choose this option to create and save overviews using gradient keograms
    overview_grad(items,channel,allrows,filedate + '_IR1grad.pdf',numdays)
    #Choose this option for original keogram overviews but with dots
    #overview_points(items,channel,allrows,filedate + '_IR1.pdf',numdays)

    return
# %%
