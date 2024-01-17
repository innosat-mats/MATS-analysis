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

start_time = DT.datetime(2023,3,22,0,0,0)
stop_time = DT.datetime(2023,4,1,0,0,0)
channel = 'IR1'
numdays = stop_time-start_time #number of days
#%%
dfSH = read_MATS_data(start_time,stop_time,filter={"TPlat":[-90,-45],'NROW': [0,400]},version='0.5',level='1b')
dfSH = dfSH[dfSH['channel'] == channel]
#%%
#dfSH.to_pickle('26to27apr')
#%%
dfNH = read_MATS_data(start_time,stop_time,filter={"TPlat":[45,90],'NROW': [0,400]},version='0.5',level='1b')
dfNH = dfNH[dfNH['channel'] == channel]
#dfNH.to_pickle('20feborb2')

#%%
items = pd.concat([dfNH,dfSH], ignore_index=True)
sortitems = items.sort_values(['EXPDate'])
sortitems.to_pickle('22to31marIR1')

# %%
ccd = pd.read_pickle(r'MatsData\15to16febIR1')

#%% Plot a chosen image
def saveSpecIm(ccditem):
    ccdimage = ccditem['ImageCalibrated']
    #saves the matrix to matlab-file
    scipy.io.savemat('aurorapic2',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
    return

# %%
image = ccd.iloc[400].ImageCalibrated
plt.pcolormesh(image,rasterized = True, vmin=0, vmax=480)

# %% Main file to run scripts
def Main():
    start_time = DT.datetime(2023,4,26,0,0,0)
    stop_time = DT.datetime(2023,4,28,0,0,0)
    numdays = stop_time-start_time #number of days
    filedate = 'apr4W'

    items = pd.read_pickle(r'MatsData\23to30aprIR1')

    #returns aurorastrips and saves the other strips in pickles.
    aurorastrips = get_strips(items,numdays,filedate)  
    get_aurora_max(aurorastrips,filedate)

    #load in all the strips from file to plot overviews
    #allstrips = pd.read_pickle(r'MatsData\apr4Wallstrips')
    #allrows = get_stripRow(allstrips)

    #overview_grad(items,channel,allrows,filename,numdays)
    #overview_points(items,channel,allrows,filename,numdays)

    return