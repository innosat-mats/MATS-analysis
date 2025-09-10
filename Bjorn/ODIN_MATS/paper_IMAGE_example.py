#%%
from mat4py import loadmat
import proplot as pplt
import pandas as pd
import numpy as np
import xarray as xr
from mats_utils.rawdata.read_data import read_MATS_data
from astropy.time import Time
from datetime import datetime, timezone
import datetime as DT
import calendar
from geopy.distance import geodesic
import cartopy.crs as ccrs
from time import strftime
from mats_utils.geolocation.coordinates import col_heights, satpos
import os
import pysolar
from mats_utils.geolocation import coordinates as coordinates
#%%
def make_ths(CCD):
    # code borrowed from Donal to produce ths
    xpixels = np.linspace(0, CCD['NCOL'], 5)
    ypixels = np.linspace(0, CCD['NROW'], 10)

    print(ypixels)
    ths = np.zeros([xpixels.shape[0], ypixels.shape[0]])
    #print (ths.shape)
    for i,col in enumerate(xpixels): 
        ths[i,:]=coordinates.col_heights(CCD,col,40,spline=True)(ypixels)
    return 0.5+xpixels,ypixels,ths.T


#%%

# daytime examples
starttime=datetime(2023,2,11,12,46)
stoptime=datetime(2023,2,11,12,47)
l1b_version="0.9"
download=False

df_day = read_MATS_data(starttime, stoptime, version=l1b_version, level='1b')
df_day = df_day[df_day['channel'] == 'IR1'].dropna().reset_index()#[0:10]

# nighttime examples
starttime=datetime(2023,1,25,3,25)
stoptime=datetime(2023,1,25,3,26)
l1b_version="0.9"
download=False

df_night = read_MATS_data(starttime, stoptime, version=l1b_version, level='1b')
df_night = df_night[df_night['channel'] == 'IR1'].dropna().reset_index()#[0:10]
#%%
df_d = df_day.iloc[3]
df_n = df_night.iloc[3]

#%%
# add heights


# %%


# plot ImageCalibrated
cmapp = 'Rocket'
fig, axs = pplt.subplots(ncols=2, nrows=1, figheight='6cm', sharey=0, figwidth='20cm',abc='a.')
cbard=axs[0].contourf(df_d['ImageCalibrated']*10**12, cmap=cmapp, levels=100,vmin=0*10**12,vmax=660*10**12)

ths = make_ths(df_d)

CS = axs[0].contour(ths[0],ths[1], ths[2],
                colors='light purple', alpha=1, levels=6)
axs[0].clabel(CS, inline=True)

axs[0].format(title='Dayglow (2023-02-11)')
axs[0].colorbar(cbard, ticks=1*10**14,label=r'Limb radiance [$m^{-2} s^{-1} str^{-1} nm^{-1}$]')


ths = make_ths(df_n)
CS = axs[1].contour(ths[0],ths[1], ths[2],
                colors='light purple', alpha=0.8, levels=8)
axs[1].clabel(CS, inline=True)

cbarn=axs[1].contourf(df_n['ImageCalibrated']*10**12, cmap=cmapp,levels=100, vmin=0*10**12, vmax=160*10**12)
axs[1].colorbar(cbarn,ticks=0.2*10**14,label=r'Limb radiance [$m^{-2} s^{-1} str^{-1} nm^{-1}$]')
axs[1].format(title='Nightglow (2023-01-25)')
fig.suptitle('MATS IR1 images')

fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/imageexamples_v09.png',format='png')
# %%
