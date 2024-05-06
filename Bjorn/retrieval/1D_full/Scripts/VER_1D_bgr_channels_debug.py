# DEBUG IMAGES
#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT

import xarray as xr
from numpy.linalg import inv
from fast_histogram import histogramdd
from bisect import bisect_left
import matplotlib.pyplot as plt
import proplot as pplt
from mats_utils.rawdata.read_data import read_MATS_data

from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
#%% ISOLATE SMALL REGION

starttime=datetime(2023,3,5,0,0)
stoptime=datetime(2023,3,5,23,59) 
dmin,dmax = 0, 95

dftop=read_MATS_data(starttime,stoptime,level="1b",version='0.5', filter={'TPsza': [dmin, dmax], 'TPlat': [22, 23]})

# %% 
ir1 = dftop[dftop['channel'] == 'IR1'].dropna().reset_index()#[0:10]
ir2 = dftop[dftop['channel'] == 'IR2'].dropna().reset_index()#[0:10]
ir3 = dftop[dftop['channel'] == 'IR3'].dropna().reset_index()#[0:10]
ir4 = dftop[dftop['channel'] == 'IR4'].dropna().reset_index()#[0:10]
# %% changes with time
fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=2,abc='a.',sharex=0)
fig.format(suptitle='Image means',ylabel='counts')

axs[0].scatter(ir1.EXPDate.values,np.mean(np.stack(ir1.ImageCalibrated)[:,:,:], axis=(1,2)),s=1, c=np.stack(ir1.afsTangentH_wgs84),cmap='viridis')
axs[0].format(title='IR1')

axs[1].scatter(ir2.EXPDate.values,np.mean(np.stack(ir2.ImageCalibrated)[:,:,:], axis=(1,2)),s=1, c=np.stack(ir2.afsTangentH_wgs84),cmap='viridis')
axs[1].format(title='IR2')

axs[2].scatter(ir3.EXPDate.values,np.mean(np.stack(ir3.ImageCalibrated)[:,:,:], axis=(1,2)),s=1, c=np.stack(ir3.afsTangentH_wgs84),cmap='viridis')
axs[2].format(title='IR3')

m=axs[3].scatter(ir4.EXPDate.values,np.mean(np.stack(ir4.ImageCalibrated)[:,:,:], axis=(1,2)),s=1, c=np.stack(ir4.afsTangentH_wgs84),cmap='viridis')
axs[3].format(title='IR4')

fig.colorbar(m, loc='b',label='afsTangentH_wgs84')

# %%
# %% as function of pointing
fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=2,abc='a.',sharex=0)
fig.format(suptitle='Image means',ylabel='counts',xlabel='afsTangentH_wgs84')

axs[0].scatter(ir1.afsTangentH_wgs84.values,np.mean(np.stack(ir1.ImageCalibrated)[:,:,:], axis=(1,2)),s=1,cmap='viridis')
axs[0].format(title='IR1')

axs[1].scatter(ir2.afsTangentH_wgs84.values,np.mean(np.stack(ir2.ImageCalibrated)[:,:,:], axis=(1,2)),s=1,cmap='viridis')
axs[1].format(title='IR2')

axs[2].scatter(ir3.afsTangentH_wgs84.values,np.mean(np.stack(ir3.ImageCalibrated)[:,:,:], axis=(1,2)),s=1,cmap='viridis')
axs[2].format(title='IR3')

m=axs[3].scatter(ir4.afsTangentH_wgs84.values,np.mean(np.stack(ir4.ImageCalibrated)[:,:,:], axis=(1,2)),s=1,cmap='viridis')
axs[3].format(title='IR4')


# %% So, the way we do the bgr profiles, how does that adapt to changes in pointing
fig, axs = pplt.subplots(figwidth='20cm',ncols=2, nrows=1,abc='a.',sharex=0)
fig.format(suptitle='limb radiance (22N to 23N) 5th March')
for i in range(0,len(ir3)):
    ch=ir1.iloc[i]
    image = image = np.stack(ch.ImageCalibrated)/ch.TEXPMS*1000
    
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[:, col-2:col+2].mean(axis=1))*1e13
    axs[0].scatter(profile,heights,c=ch.afsTangentH_wgs84*np.ones(len(profile)), vmin=86, vmax=90, s=1,alpha=1,cmap='viridis')
    axs[0].format(title='IR1')
    #axs[0].scatter(profile,heights,c=ch.TPsza*np.ones(len(profile)), vmin=81.5,vmax=82,s=0.7,alpha=0.5,cmap='viridis')
    #fig.colorbar(loc='b')

    ch=ir2.iloc[i]
    image = image = np.stack(ch.ImageCalibrated)/ch.TEXPMS*1000
    
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[:, col-2:col+2].mean(axis=1))*1e13
    axs[1].scatter(profile,heights,c=ch.afsTangentH_wgs84*np.ones(len(profile)), vmin=86, vmax=90, s=1,alpha=1,cmap='viridis')
    axs[1].format(title='IR2')
    #axs[0].scatter(profile,heights,c=ch.TPsza*np.ones(len(profile)), vmin=81.5,vmax=82,s=0.7,alpha=0.5,cmap='viridis')
    #fig.colorbar(loc='b')

    """
    ch=ir1.iloc[i]
    image = image = np.stack(ch.ImageCalibrated)/ch.TEXPMS*1000
    
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[:, col-2:col+2].mean(axis=1))*1e13
    axs[2].scatter(profile,heights,c=ch.afsTangentH_wgs84*np.ones(len(profile)), vmin=86, vmax=90, s=1,alpha=1,cmap='viridis')
    axs[2].format(title='IR1')
    #axs[0].scatter(profile,heights,c=ch.TPsza*np.ones(len(profile)), vmin=81.5,vmax=82,s=0.7,alpha=0.5,cmap='viridis')
    #fig.colorbar(loc='b')
    """
import matplotlib as mpl
y = np.stack(ir1.afsTangentH_wgs84)[:,0]
cmap = mpl.cm.get_cmap('viridis')
norm = mpl.colors.Normalize(vmin=86, vmax=90)
colors = cmap(norm(y))
labels = np.linspace(86, 90, num=5, endpoint=True)

cbar=fig.colorbar(cmap, loc='r',values=np.linspace(86, 90, endpoint=True, num=100),label='afsTangentH_wgs84')
cbar.ax.set_yticklabels(labels)
fig.savefig('bgr_pointing_IR1IR2.png',format='png')

# %%
fig, axs = pplt.subplots(figwidth='20cm',ncols=2, nrows=2,abc='a.',sharex=0)
fig.format(suptitle='middle profiles',ylabel='counts')

for i in range(0,len(ir1)):
    axs[0].plot(np.stack(ir1.ImageCalibrated)[i,:,4])
    #axs[2].plot(np.mean(np.stack(ir3.ImageCalibrated)[i,:,4:7],axis=2))
    #axs[3].plot(np.mean(np.stack(ir4.ImageCalibrated)[i,:,4:7],axis=2))

# %%
