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
import copy

#%%
def plot_background(p3, p4, label_suffix):
    axs[bgr_index].plot(p3, dstest.altitude.values[10:endcut], color='blue', label=f'IR3{label_suffix}')
    axs[bgr_index].plot(p4, dstest.altitude.values[10:endcut], color='purple', label=f'IR4{label_suffix}')

def flatten_and_reverse(array, flip):
    """Flatten a nested array and reverse it if flip is True."""
    flat_list = [x for xs in array for x in xs]  # Flatten the array
    return flat_list[::-1] if flip else flat_list

def process_channel(df, channel):
    """Filter, drop duplicates, and set index for the given channel."""
    filtered_df = df[df['channel'] == channel].dropna().reset_index()
    filtered_df = filtered_df.drop_duplicates(subset=['EXPDate']).set_index('EXPDate')
    return filtered_df

def calc_distance(coord1, coord2):
    '''function to calculate distance between the two satellite positions'''
    distance = geodesic(coord1, coord2).kilometers
    return distance # km

def prepare_measurment(ch,ir3,ir4,offset,subtract=True,endcut=-25):

    z1,p1=prepare_profile(ch,offset)
    _,p3=prepare_profile(ir3,offset=0)
    _,p4=prepare_profile(ir4,offset=0)

    p3 = p3[10:endcut]
    p4 = p4[10:endcut]

    p3=p3-p3[-4:].mean()/1.05
    p4=p4-p4[-4:].mean()/1.05

    if subtract:

        p1=p1[10:endcut]-(p3+p4)/2

    return np.array(z1),np.array(p1),np.array(p3),np.array(p4)

def prepare_profile(ch,offset):
    # This function averages some columns and
    # calculates tangent heights
    
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2) - offset
    #cs = col_heights(ch, col, 10, spline=True) # TEST BETTER 
    cs = col_heights(ch, col, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    # multiply with factor to get right values (depends on version?)
    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e12)
    #profile = np.array(image[:, col]*1e12)


    # set heights
    common_heights = np.arange(60,110,0.25)*1000
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile

#%%

# NLC
NLC=False
DAYGLOW=False
NIGHTGLOW=True

# CHANNEL
chan_str='IR1' # DO NOT CHANGE

OSIRIS_stray=False
#bgr_removal=False # FOR BOTH OSIRIS AND MATS
plt_advanced=True # plot
plt_bgr_channels=False

if NLC:
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/OSIRIS_NLCscans.mat')
    savestr='NLC'
# AIRGLOW
else:    
    data = loadmat('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/Odin/sm_ABand_2022_2023_v3.mat')
    if DAYGLOW:
        savestr='DAYGLOW'
        txt_fpath=f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/DAYGLOW/'
    if NIGHTGLOW:
        savestr='NIGHTGLOW'
        txt_fpath=f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/AIRGLOW/{chan_str}/NIGHTGLOW/'

# %%
# sattime lat lon
if NLC:
    lat = np.array(data['lat']).flatten()
    lon = np.array(data['lon']).flatten()
    timeO = np.array(data['t']).flatten()
else:
    lat = np.array(data['sm']['lat90'])
    lon = np.array(data['sm']['lon90'])
    timeO = np.array(data['sm']['time90'])
    time = np.array(timeO)

jd = timeO + 2400000.5
t = Time(jd, format='jd')
timeO = t.to_datetime()

# SZA at 90
sza = np.zeros(len(timeO))
for i in range(len(sza)):
    sza[i] = 90 - pysolar.solar.get_altitude(
        lat[i],
        lon[i],
        pd.Timestamp(timeO[i], tzinfo=timezone.utc).to_pydatetime()
    )

#%%
#  GENERATE MEANS FROM THE ISOLATED DATA
# OSIRIS NEEDS TO BE PUT ON A COMMON GRID
#OSheights = np.arange(59.25,109.25,0.25) # SHIFT
OSheights = np.arange(60,100,0.25)


# messy OSIRIS stuff
# Initialize OSIRIS output arrays
n_mjd = len(data['sm']['mjd'])
OSinter = np.zeros((n_mjd, len(OSheights)))
OSinter2 = np.zeros((n_mjd, len(OSheights)))
OSinter3 = np.zeros((n_mjd, len(OSheights)))
OSinter4 = np.zeros((n_mjd, len(OSheights)))
dz = np.zeros(n_mjd)
zref = np.zeros(n_mjd)
exptime = np.zeros(n_mjd)
zmax = np.zeros(n_mjd)
top_noise_L1 = np.zeros(n_mjd)
top_noise_L2 = np.zeros(n_mjd)
top_noise_L3 = np.zeros(n_mjd)
top_noise_L4 = np.zeros(n_mjd)


# Process each time step
for i in range(n_mjd):
    zt = data['sm']['zt'][i]
    if type(data['sm']['exposure_time'][i]) is list:
        exptime[i] = exptime[i] 
    else: 
        exptime[i] = data['sm']['exposure_time'][i]
    flip = False  # For sorting

    # Flatten zt and check if reversal is needed
    zt = flatten_and_reverse(zt, flip)
    if zt[0] > zt[1]:
        flip = True
        zt = zt[::-1]

    # approx dz at 90 (sometimes)
    dz[i] = zt[-15] - zt[-16]
    zref[i] = zt[-15]

    # SELECT CHANNEL FROM ODIN
    Li1 = flatten_and_reverse(data['sm']['L1'][i], flip)
    Li2 = flatten_and_reverse(data['sm']['L2'][i], flip)

    # remove mean of top altitudes if NIGHTGLOW (straylight):
    #if NIGHTGLOW:
    #    print(np.mean(Li2[-3:])/10**12)
    #    Li1 = Li1 - np.mean(Li1[-3:])
    #    Li2 = Li2 - np.mean(Li2[-3:])

    zmax[i] = zt[-1]        

    #top_noise_L1[i] = Li1[-1]
    #top_noise_L2[i] = Li2[-1]

    top_noise_L1[i] = np.mean(Li1[-3:])
    top_noise_L2[i] = np.mean(Li2[-3:])

    # BACKGROUND CHANNELS
    Li3 = flatten_and_reverse(data['sm']['L3'][i], flip)
    Li4 = flatten_and_reverse(data['sm']['L4'][i], flip)

    top_noise_L3[i] = np.mean(Li3[-3:])
    top_noise_L4[i] = np.mean(Li4[-3:])

    #top_noise_L3[i] = Li3[-1]
    #top_noise_L4[i] = Li4[-1]

    # Interpolate to OSheights
    OSinter[i, :] = np.interp(OSheights, zt, Li1)
    OSinter2[i, :] = np.interp(OSheights, zt, Li2)
    OSinter3[i, :] = np.interp(OSheights, zt, Li3)
    OSinter4[i, :] = np.interp(OSheights, zt, Li4)

    # Apply stray correction if needed
    OSIRIS_stray = False
    if OSIRIS_stray:
        OSinter3[i, :] -= OSinter3[i, -4:].mean()
        OSinter4[i, :] -= OSinter4[i, -4:].mean()


# SUBTRACT BACKGROUND
#OSheights = OSheights + 0.75 # SHIFT
common_heights = OSheights*1000

ds_sza = xr.Dataset(
    data_vars=dict(
        dz = (["time"], dz),
        zref = (["time"], zref),
        zmaxs = (["time"], zmax),
        exptime = (["time"], exptime),
        lat=(["time"], lat),
        lon=(["time"], lon),
        sza=(["time"], sza),
        L1=(["time", "altitude"], OSinter),
        L2=(["time", "altitude"], OSinter2),
        L3=(["time", "altitude"], OSinter3),
        L4=(["time", "altitude"], OSinter4),
        top_noise_L1=(["time"], top_noise_L1),
        top_noise_L2=(["time"], top_noise_L2),
        top_noise_L3=(["time"], top_noise_L3),
        top_noise_L4=(["time"], top_noise_L4),

    ),
    coords=dict(
        time=timeO,
        altitude=OSheights
    ),
    attrs=dict(description="ODIN"),
)

#%%
ALL=False
# Apply conditions for NIGHTGLOW and DAYGLOW
if NIGHTGLOW:
    ds = ds_sza.where(ds_sza.sza > 100)
if DAYGLOW:
    ds = ds_sza.where(ds_sza.sza < 90)
if ALL:
    ds = ds_sza

#%%
ds = ds.sel(time=slice("2022-01-01", "2023-02-01")).dropna(dim='time')

ds.L1[:,:].T.plot.line(y='altitude',add_legend=False, hue='time')
# isolate data where top altitude is between 100 and 103 km
#ds = ds.where(ds.zmaxs > 100).where(ds.zmaxs < 110)

# isolate data for latitudes below -70
#ds = ds.where(ds.lat > 70)

# %%

# plot correlation between OSIRIS signal at 90 km and top of OSIRIS measurement
# for L1 and L2

fig, axs = pplt.subplots(ncols=2, nrows=1, share=0)
axs.format(suptitle='OSIRIS signal at 90 km vs top of INTERPOLATED OSIRIS measurement')
axs.format(xlabel='OSIRIS signal at 90 km [1e12]')
axs.format(ylabel='Top of OSIRIS measurement [km]')

ax = axs[0]
for i in range(0,len(ds.time)):
    axs[0].scatter(np.mean(ds['L1'][i,110:125].values), np.mean(ds['L1'][i,-35:-30].values), s=4, color='red', alpha=0.7)
    axs[1].scatter(np.mean(ds['L2'][i,110:125].values), np.mean(ds['L2'][i,-35:-30].values), s=4, color='red', alpha=0.7)

axs[0].format(ylim=[0, 0.2*1e14])
axs[1].format(ylim=[0, 0.2*1e14])

axs[0].format(title='L1')
axs[1].format(title='L2')
# %%
# %%

# plot correlation between OSIRIS signal at 90 km and top of OSIRIS measurement
# for L1 and L2

fig, axs = pplt.subplots(ncols=2, nrows=1, share=0)
axs.format(suptitle='OSIRIS signal at 90 km vs top of OSIRIS measurement')
axs.format(xlabel='OSIRIS signal at 90 km')
axs.format(ylabel='Top of OSIRIS measurement ')

ax = axs[0]

m=axs[0].scatter(np.mean(ds['L1'][:,125:140].values, axis=1), ds.top_noise_L1[:], s=3, c=ds.zmaxs[:], alpha=0.9,cmap='viridis')
axs[1].scatter(np.mean(ds['L2'][:,125:140].values,axis=1), ds.top_noise_L2[:], s=3, c=ds.zmaxs[:], alpha=0.9,cmap='viridis')

fig.colorbar(m,loc='r', label='Top of OSIRIS measurement [km]')
axs[0].format(ylim=[-1*1e13, 2*1e13],xlim=[0,1e14])
axs[1].format(ylim=[-1*1e13, 2*1e13],xlim=[0,1e14])

axs[0].format(title='L1')
axs[1].format(title='L2')
# %%

# plot correlation between OSIRIS signal at 90 km and top of OSIRIS measurement
# for L1 and L2

fig, axs = pplt.subplots(ncols=2, nrows=1, share=0)
axs.format(suptitle='top signal / signal at 90 km')
axs.format(xlabel='OSIRIS signal at 90 km')
#axs.format(ylabel='')

ax = axs[0]
for i in range(0,len(ds.time)):
    axs[0].scatter(np.mean(ds['L1'][i,110:125].values), ds.top_noise_L1[i]/(np.mean(ds['L1'][i,110:115].values)), s=4, color='blue', alpha=0.7)
    axs[1].scatter(np.mean(ds['L2'][i,110:125].values), ds.top_noise_L2[i]/(np.mean(ds['L2'][i,110:115].values)), s=4, color='blue', alpha=0.7)

axs[0].format(ylim=[0, 0.6])
axs[1].format(ylim=[0, 0.6])

axs[0].format(title='L1')
axs[1].format(title='L2')
# %%

# plot correlation between OSIRIS signal at 90 km and top of OSIRIS measurement
# for L3 and L4

fig, axs = pplt.subplots(ncols=2, nrows=1, share=0)
axs.format(suptitle='OSIRIS signal at 90 km vs top of OSIRIS measurement')
axs.format(xlabel='OSIRIS signal at 90 km ')
axs.format(ylabel='OSIRIS top signal')

ax = axs[0]

m=axs[0].scatter(np.mean(ds['L3'][:,110:125].values, axis=1), ds.top_noise_L3[:], s=3, c=ds.zmaxs[:], alpha=0.9,cmap='viridis')
axs[1].scatter(np.mean(ds['L4'][:,110:125].values,axis=1), ds.top_noise_L4[:], s=3, c=ds.zmaxs[:], alpha=0.9,cmap='viridis')

fig.colorbar(m,loc='r', label='Top height')
axs[0].format(ylim=[0, 1.2*1e13])
axs[1].format(ylim=[0, 1.2*1e13])

axs[0].format(xlim=[0, 1.2*1e13])
axs[1].format(xlim=[0, 1.2*1e13])

# plot 1:1 line
axs[0].plot([0, 1.2*1e13], [0, 1.2*1e13], color='black')
axs[1].plot([0, 1.2*1e13], [0, 1.2*1e13], color='black')

axs[0].format(title='L3')
axs[1].format(title='L4')
# %%
# plot L1 and L2 signal at top as a function of time

fig, axs = pplt.subplots(ncols=2, nrows=1, share=0)

# remove values above 2 * 10^13
ds['top_noise_L1'] = ds['top_noise_L1'].where(ds['top_noise_L1'] < 2*1e13)
ds['top_noise_L2'] = ds['top_noise_L2'].where(ds['top_noise_L2'] < 2*1e13)

# apply rolling means top ds top noise
#ds['top_noise_L1'] = ds['top_noise_L1'].rolling(time=10, center=True).mean()
#ds['top_noise_L2'] = ds['top_noise_L2'].rolling(time=10, center=True).mean()

m=axs[0].scatter(ds.time, ds.top_noise_L1[:], color=ds.zmaxs, label='L1',s=3,cmap='viridis',vmax=106, vmin=100)
axs[1].scatter(ds.time, ds.top_noise_L2[:], color=ds.zmaxs, label='L2',s=3, cmap='viridis',vmax=106, vmin=100)

fig.colorbar(m,loc='r', label='Top height')
axs[0].format(ylim=[0, 0.25*1e14])
axs[1].format(ylim=[0, 0.25*1e14])
# %%

# plot L1, L2, L3 and L4 signal at top as a function of time

fig, axs = pplt.subplots(ncols=1, nrows=1, sharex=1)

# remove values above 2 * 10^13
#ds['top_noise_L1'] = ds['top_noise_L1'].where(ds['top_noise_L1'] < 4*1e13)
#ds['top_noise_L2'] = ds['top_noise_L2'].where(ds['top_noise_L2'] < 4*1e13)
#ds['top_noise_L3'] = ds['top_noise_L3'].where(ds['top_noise_L3'] < 4*1e13)
#ds['top_noise_L4'] = ds['top_noise_L4'].where(ds['top_noise_L4'] < 4*1e13)

# apply rolling means top ds top noise
#ds['top_noise_L1'] = ds['top_noise_L1'].rolling(time=10, center=True).mean()
#ds['top_noise_L2'] = ds['top_noise_L2'].rolling(time=10, center=True).mean()

axs[0].scatter(ds.time, np.mean(ds['L1'][:,-40:-39].values, axis=1), color='r', label='L1',s=3,cmap='viridis',vmax=106, vmin=100)
axs[0].scatter(ds.time, np.mean(ds['L2'][:,-40:-39].values, axis=1), color='b',s=3, label='L2', cmap='viridis',vmax=106, vmin=100)
axs[0].scatter(ds.time, np.mean(ds['L3'][:,-40:-39].values, axis=1), color='g',s=3,label='L3', cmap='viridis',vmax=106, vmin=100)
axs[0].scatter(ds.time, np.mean(ds['L4'][:,-40:-39].values, axis=1), color='y',s=3, label='L4',cmap='viridis',vmax=106, vmin=100)


axs[0].format(ylim=[0, 0.2*1e14])
#axs[1].format(ylim=[0, 0.25*1e14])
axs[0].legend(ncol=4,loc='bottom')
axs[0].format(xlabel=' ')
fig.format(title='Signal at 100 km (100 < Top < 110 km)',fontsize=8)

# %%

if NIGHTGLOW:
    ds = ds_sza.where(ds_sza.sza > 100)
if DAYGLOW:
    ds = ds_sza.where(ds_sza.sza < 90)

ds = ds.sel(time=slice("2021-10-01", "2023-03-01")).dropna(dim='time')

# isolate data where top altitude is between 100 and 103 km
dszmin, dszmax = 101, 105
ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)

# isolate data for latitudes below -70
#ds = ds.where(ds.lat < 70)

# plot L1, L2, L3 and L4 signal at top_noise as a function of time

fig, axs = pplt.subplots(ncols=1, nrows=1, sharex=1)
fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)
axs.scatter(ds.time, ds.top_noise_L1[:], color='r', label='L1',s=3)
axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')
axs.scatter(ds.time, ds.top_noise_L3[:], color='g',s=3,label='L3')
axs.scatter(ds.time, ds.top_noise_L4[:], color='y',s=3, label='L4')
axs.format(ylabel='Top signal')
axs[0].legend(ncol=4,loc='bottom')
axs.format(ylim=[-0.2*1e14, 0.2*1e14])
# %%
if NIGHTGLOW:
    ds = ds_sza.where(ds_sza.sza > 100)
if DAYGLOW:
    ds = ds_sza.where(ds_sza.sza < 90)

ds = ds.sel(time=slice("2021-10-01", "2023-03-01")).dropna(dim='time')

# isolate data where top altitude is between 100 and 103 km
dszmin, dszmax = 70, 110
ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)

# isolate data for latitudes below -70
#ds = ds.where(ds.lat < 70)

fig, axs = pplt.subplots(ncols=2, nrows=1, sharex=0,sharey=0)
fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)
axs[0].scatter(ds.top_noise_L1[:], ds.top_noise_L2[:], color='r', label='L1',s=3)
#axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

axs[0].format(ylabel='Top signal')
#axs[0].legend(ncol=4,loc='bottom')
# plot straight line
axs[0].plot([0, 0.2*1e14], [0, 0.2*1e14], color='black')
axs[0].format(xlim=[0, 0.15*1e14], ylim=[0, 0.15*1e14],
           xlabel='Top signal L1', ylabel='Top signal L2')

axs[1].scatter(ds.top_noise_L3[:], ds.top_noise_L4[:], color='r', label='L3',s=3)
axs[1].format(xlim=[0, 0.15*1e14], ylim=[0, 0.15*1e14],
           xlabel='Top signal L3', ylabel='Top signal L4')
# plot straight line
axs[1].plot([0, 0.2*1e14], [0, 0.2*1e14], color='black')

#axs.format(ylim=[0, 0.2*1e14])

# %%
if NIGHTGLOW:
    ds = ds_sza.where(ds_sza.sza > 100)
if DAYGLOW:
    ds = ds_sza.where(ds_sza.sza < 90)

ds = ds.sel(time=slice("2021-10-01", "2023-03-01")).dropna(dim='time')

# isolate data where top altitude is between 100 and 103 km
dszmin, dszmax = 70, 110
ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)

# isolate data for latitudes below -70
#ds = ds.where(ds.lat < 70)

fig, axs = pplt.subplots(ncols=2, nrows=1, sharex=0,sharey=0)
fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)
axs[0].scatter(ds.top_noise_L1[:], ds.top_noise_L3[:], color='r', label='L1',s=3)
#axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

axs[0].format(ylabel='Top signal')
#axs[0].legend(ncol=4,loc='bottom')
# plot straight line
axs[0].plot([0, 0.2*1e14], [0, 0.2*1e14], color='black')
axs[0].format(xlim=[0, 0.15*1e14], ylim=[0, 0.15*1e14],
           xlabel='Top signal L1', ylabel='Top signal L2')

axs[1].scatter(ds.top_noise_L2[:], ds.top_noise_L4[:], color='r', label='L3',s=3)
axs[1].format(xlim=[0, 0.15*1e14], ylim=[0, 0.15*1e14],
           xlabel='Top signal L2', ylabel='Top signal L4')
# plot straight line
axs[1].plot([0, 0.2*1e14], [0, 0.2*1e14], color='black')

#axs.format(ylim=[0, 0.2*1e14])
# %%
fig, axs = pplt.subplots(ncols=2, nrows=2,sharex=0,sharey=0)
fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)

ds = ds_sza.where(ds_sza.sza > 90)

dszmin, dszmax = 100, 120
ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)

axs[0].scatter(ds.time, ds.top_noise_L3, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')
axs[1].scatter(ds.time, ds.top_noise_L4, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')

axs[2].scatter(ds.time, ds.top_noise_L1, c=ds.sza,label={sza},vmin=87,vmax=110, s=1, extend='both')
m=axs[3].scatter(ds.time, ds.top_noise_L2, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')

#axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

axs[0].format(ylabel='O3')
axs[1].format(ylabel='O4')
axs[2].format(ylabel='O1')
axs[3].format(ylabel='O2')
#axs[0].legend(ncol=4,loc='bottom')
# plot straight line
axs[0].format(ylim=[0,0.1*1e14])
axs[1].format(ylim=[0,0.1*1e14])

axs[2].format(ylim=[0,0.3*1e14])
axs[3].format(ylim=[0,0.3*1e14])

fig.colorbar(m,loc='b', label='SZA')

#axs.format(xlim=[0,0.5*1e14])



# %%

dszmin, dszmax = 100, 120
fig, axs = pplt.subplots(ncols=2, nrows=2,sharex=0,sharey=0)
fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)

ds = ds_sza.where(ds_sza.sza > 100)
ds = ds.sel(time=slice("2023-01-01", "2023-02-01")).dropna(dim='time')

ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)

axs[0].scatter(ds.zmaxs, ds.top_noise_L3, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')
axs[1].scatter(ds.zmaxs, ds.top_noise_L4, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')

axs[2].scatter(ds.zmaxs, ds.top_noise_L1, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')
m=axs[3].scatter(ds.zmaxs, ds.top_noise_L2, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')
#axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

# horizontal lines in axs[2] and axs[3] at 0.3*1e13
axs[2].hlines(0.65*1e13, 100, 120, color='black',alpha=0.3)
axs[3].hlines(0.65*1e13, 100, 120, color='black',alpha=0.3)


axs[2].format(ylabel='O1')
axs[3].format(ylabel='O2')

#axs[0].legend(ncol=4,loc='bottom')
# plot straight line
axs[0].format(ylim=[0,0.2*1e14], xlabel='Top height [km]')
axs[1].format(ylim=[0,0.2*1e14], xlabel='Top height [km]')

axs[2].format(ylim=[0,0.5*1e14], xlabel='Top height [km]')
axs[3].format(ylim=[0,0.5*1e14], xlabel='Top height [km]')

fig.colorbar(m,loc='b', label='SZA')

#axs.format(xlim=[0,0.5*1e14])

fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/correlate_top_OSIRIS_ng_TOPHEIGHT.png')
# %%

dszmin, dszmax = 100, 120
fig, axs = pplt.subplots(ncols=2, nrows=1,sharex=0,sharey=0,figwidth='12cm')

fig.format(suptitle=f'Top of measurement signal',fontsize=8)
axs[0].format(title=f'OSIRIS IR1',fontsize=8)
axs[1].format(title=f'OSIRIS IR2',fontsize=8)
ds = ds_sza.where(ds_sza.sza > 100)
ds = ds.sel(time=slice("2022-12-12", "2023-02-01")).dropna(dim='time')

ds = ds.where(ds.zmaxs > dszmin).where(ds.zmaxs < dszmax)
#axs[0].scatter(ds.zmaxs, ds.top_noise_L3, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')
#axs[1].scatter(ds.zmaxs, ds.top_noise_L4, c=ds.sza,vmin=87,vmax=110, s=1, extend='both')

axs[0].scatter(ds.zmaxs, ds.top_noise_L1,c='blue',vmin=87,vmax=110, s=1, extend='both')
m=axs[1].scatter(ds.zmaxs, ds.top_noise_L2,c='blue',vmin=87,vmax=110, s=1, extend='both')
#axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

# horizontal lines in axs[2] and axs[3] at 0.3*1e13
axs[0].hlines(0.6*1e13, 100, 120, color='red',alpha=0.4)
axs[1].hlines(0.6*1e13, 100, 120, color='red',alpha=0.4)

axs[0].format(ylabel='Top signal '+r'$[m^{-2} s^{-1} str^{-1} nm^{-1}]$',)
#axs[1].format(ylabel='Top signal '+r'$[m^{-2} s^{-1} str^{-1} nm^{-1}]$',)

axs[0].format(ylim=[0,0.5*1e14], xlabel='Tangent point height [km]')
axs[1].format(ylim=[0,0.5*1e14], xlabel='Tangent point height [km]')

#fig.colorbar(m,loc='b', label='SZA')

#axs.format(xlim=[0,0.5*1e14])

fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/correlate_top_OSIRIS_ng_TOPHEIGHT_IR1_IR2.png')

# %%
fig, axs = pplt.subplots(ncols=2, nrows=2,sharex=0,sharey=0)
#fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)

dsm = ds_sza.where(ds_sza.sza > 100)

colors = ['r','b','g','y']
cc = 0
for i in [100]:
#for i in [103]
    dszmin, dszmax = i, i+1
    ds = dsm.where(dsm.zmaxs > dszmin).where(dsm.zmaxs < dszmax).dropna(dim='time')
    ds= ds.where(ds.top_noise_L1 < 2*1e13).dropna(dim='time')

    #axs[0].scatter(ds.zmaxs, ds.top_noise_L3, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')
    #axs[1].scatter(ds.zmaxs, ds.top_noise_L4, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')

    axs[0].scatter(np.nanmean(ds.L1[:,110:130].values,axis=1)-ds.top_noise_L4, ds.top_noise_L1[:]-ds.top_noise_L4, s=2, color=colors[cc], extend='both')
    axs[1].scatter(np.nanmean(ds.L2[:,110:130].values,axis=1)-ds.top_noise_L4, ds.top_noise_L2[:]-ds.top_noise_L4, s=2, color=colors[cc], extend='both')
    #m=axs[3].scatter(ds.zmaxs, ds.top_noise_L2, c=ds.sza,label={sza},vmin=87,vmax=110, s=1, extend='both')
    #axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

    #axs[0].format(ylabel='O3')
    #axs[1].format(ylabel='O4')

    # plot linear fit
    x = np.nanmean(ds.L1[:,110:130].values,axis=1)-ds.top_noise_L4
    y = ds.top_noise_L1[:]-ds.top_noise_L4
    m, b = np.polyfit(x, y, 1)
    axs[0].plot(x, m*x + b, color=colors[cc], label=f'{i} < Top < {i+1}')
    print(m)

    x = np.nanmean(ds.L2[:,110:130].values,axis=1)-ds.top_noise_L4
    y = ds.top_noise_L2[:]-ds.top_noise_L4
    m, b = np.polyfit(x, y, 1)
    axs[1].plot(x.values, m*x.values + b, color=colors[cc])
    print(m)

    #axs[0].legend(ncol=4,loc='bottom')
    # plot straight line
    #axs[0].format(ylim=[0,0.3*1e14], xlabel='Top height [km]')
    #axs[1].format(ylim=[0,0.3*1e14], xlabel='Top height [km]')
    axs[0].format(xlabel='Signal around 87 km',ylabel='Top signal', title='O1')
    axs[1].format(xlabel='Signal around 87 km',ylabel='Top signal',title='O2')


    # plot fractions in axs[2] and axs[3]
    axs[2].scatter(np.nanmean(ds.L1[:,110:130].values,axis=1)-ds.top_noise_L4.values, (ds.top_noise_L1[:]-ds.top_noise_L4)/(np.mean(ds.L1[:,110:130].values,axis=1)-ds.top_noise_L4), s=2, color=colors[cc], extend='both')
    axs[3].scatter(np.nanmean(ds.L2[:,110:130].values,axis=1)-ds.top_noise_L4.values, (ds.top_noise_L2[:]-ds.top_noise_L4)/(np.mean(ds.L2[:,110:130].values,axis=1)-ds.top_noise_L4), s=2, color=colors[cc], extend='both')

    cc += 1

fig.format(suptitle='top vs 87 km for different top altitudes',fontsize=8)

#fig.colorbar(m,loc='b', label='SZA')
fig.legend(ncol=3,loc='bottom')
#axs[1].legend(ncol=3,loc='bottom')
#axs.format(ylim=[0,0.1*1e14])
axs.format(xlim=[1.e13,1e14])


fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/correlate_top_OSIRIS_ng_TOPHEIGHT.png')

# %%
# %%
fig, axs = pplt.subplots(ncols=2, nrows=1,sharex=0,sharey=0)
#fig.format(title=f'Top signal ({dszmin} < Top < {dszmax} km)',fontsize=8)

dsm = ds_sza.where(ds_sza.sza >100)

colors = ['r','b','g','y']
cc = 0
for i in [100,102]:
    dszmin, dszmax = i, i+1.5
    ds = dsm.where(dsm.zmaxs > dszmin).where(dsm.zmaxs < dszmax)

    #axs[0].scatter(ds.zmaxs, ds.top_noise_L3, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')
    #axs[1].scatter(ds.zmaxs, ds.top_noise_L4, c=ds.sza, label={sza},vmin=87,vmax=110, s=1, extend='both')

    axs[0].scatter(np.max(ds.L1[:,:].values,axis=1), (ds.top_noise_L1[:]-ds.top_noise_L3[:])/(np.mean(ds.L1[:,110:120].values,axis=1)), s=1, color=colors[cc], extend='both')
    axs[1].scatter(np.max(ds.L2[:,:].values,axis=1), (ds.top_noise_L2[:]-ds.top_noise_L3[:])/(np.mean(ds.L2[:,110:120].values,axis=1)), s=1, color=colors[cc], extend='both')
    #m=axs[3].scatter(ds.zmaxs, ds.top_noise_L2, c=ds.sza,label={sza},vmin=87,vmax=110, s=1, extend='both')
    #axs.scatter(ds.time, ds.top_noise_L2[:], color='b',s=3, label='L2')

    axs[0].format(ylabel='O3')
    axs[1].format(ylabel='O4')

    #axs[0].legend(ncol=4,loc='bottom')
    # plot straight line
    #axs[0].format(ylim=[0,0.3*1e14], xlabel='Top height [km]')
    #axs[1].format(ylim=[0,0.3*1e14], xlabel='Top height [km]')
    axs[0].format(xlabel='Signal around 87 km',ylabel='Top signal', title='O1')
    axs[1].format(xlabel='Signal around 87 km',ylabel='Top signal',title='O2')


    cc += 1

fig.format(suptitle='(100 < Top < 101.5 km)',fontsize=8)
axs[0].format(ylim=[0,1])
axs[1].format(ylim=[0,1])
#fig.colorbar(m,loc='b', label='SZA')

#axs.format(xlim=[0,0.5*1e14])

fig.savefig(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/correlate_top_OSIRIS_ng_TOPHEIGHT.png')

# %%
