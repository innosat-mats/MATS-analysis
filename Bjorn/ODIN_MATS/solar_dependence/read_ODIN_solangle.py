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

#%%

def flatten_and_reverse(array, flip):
    """Flatten a nested array and reverse it if flip is True."""
    flat_list = [x for xs in array for x in xs]  # Flatten the array
    return flat_list[::-1] if flip else flat_list

def process_channel(df, channel):
    """Filter, drop duplicates, and set index for the given channel."""
    filtered_df = df[df['channel'] == channel].dropna().reset_index()
    filtered_df = filtered_df.drop_duplicates(subset=['EXPDate']).set_index('EXPDate')
    return filtered_df

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
    #profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e12)
    profile = np.array(image[:, col]*1e12)


    # set heights
    common_heights = np.arange(60,110,0.25)*1000
    profile=np.interp(common_heights,heights,profile)
    return common_heights, profile

#%%

# NLC
NLC=False
DAYGLOW=True
NIGHTGLOW=False

# CHANNEL
chan_str='IR2' # DO NOT CHANGE

OSIRIS_stray=True
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
#OSheights = np.arange(61.25,111.25,0.25) # SHIFT
OSheights = np.arange(60,110,0.25)
common_heights = np.arange(60,110,0.25)*1000

# messy OSIRIS stuff
# Initialize OSIRIS output arrays
n_mjd = len(data['sm']['mjd'])
OSinter = np.zeros((n_mjd, len(OSheights)))
OSinter3 = np.zeros((n_mjd, len(OSheights)))
OSinter4 = np.zeros((n_mjd, len(OSheights)))

# Process each time step
for i in range(n_mjd):
    zt = data['sm']['zt'][i]
    flip = False  # For sorting

    # Flatten zt and check if reversal is needed
    zt = flatten_and_reverse(zt, flip)
    if zt[0] > zt[1]:
        flip = True
        zt = zt[::-1]

    # SELECT CHANNEL FROM ODIN
    Li = data['sm'][f'L{1 if chan_str == "IR1" else 2}'][i]
    Li = flatten_and_reverse(Li, flip)

    # BACKGROUND CHANNELS
    Li3 = flatten_and_reverse(data['sm']['L3'][i], flip)
    Li4 = flatten_and_reverse(data['sm']['L4'][i], flip)

    # Interpolate to OSheights
    OSinter[i, :] = np.interp(OSheights, zt, Li)
    OSinter3[i, :] = np.interp(OSheights, zt, Li3)
    OSinter4[i, :] = np.interp(OSheights, zt, Li4)

    # Apply stray correction if needed
    if OSIRIS_stray:
        OSinter3[i, :] -= OSinter3[i, -4:].mean()
        OSinter4[i, :] -= OSinter4[i, -4:].mean()

# SUBTRACT BACKGROUND
#OSheights = OSheights - 1.25 # SHIFT
ds = xr.Dataset(
    data_vars=dict(
        lat=(["time"], lat),
        lon=(["time"], lon),
        sza=(["time"], sza),
        Li=(["time", "altitude"], OSinter),
        L3=(["time", "altitude"], OSinter3),
        L4=(["time", "altitude"], OSinter4),

    ),
    coords=dict(
        time=timeO,
        altitude=OSheights
    ),
    attrs=dict(description="ODIN"),
)

#%%
# Apply conditions for NIGHTGLOW and DAYGLOW
if NIGHTGLOW:
    ds = ds.where(ds.sza > 98)
if DAYGLOW:
    ds = ds.where(ds.sza < 97)

# Select data within specified date range and drop NaNs
ds = ds.sel(time=slice("2023-02-01", "2023-03-01")).dropna(dim='time')
starttime=datetime(2023,1,1,0,0)
stoptime=datetime(2023,2,1,0,0)
l1b_version="0.6"
download=True
run=True

if download:
    # Download MATS data
    starttime=datetime(2023,1,14,0,0)
    stoptime=datetime(2023,1,17,0,0)
    dmin,dmax = 0, 97
    ccdsel0, ccdsel1 = 1, 5
    filter={'TPsza': [dmin, dmax], 'CCDSEL': [ccdsel0,ccdsel1]}
    dftop_full=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)


else:
    dftop_full=pd.read_pickle(f'/media/waves/AVAGO/data/MATS/MATS_ODIN_DATA/{savestr}/{str(starttime)}_{str(stoptime)}.pkl')

dftop_full=dftop_full.drop(labels=['SID','schedule_description_short'],axis=1)
dftop_full['EXPDate'] = pd.to_datetime(dftop_full['EXPDate'])
dftop_full=dftop_full[0:-1:25]
#%%
# Process IR3 and IR4 channels
ir3 = process_channel(dftop_full, 'IR3')
ir4 = process_channel(dftop_full, 'IR4')
# Select channel for analysis
dftop = process_channel(dftop_full, chan_str)

#%%
endcut=-25

# mean of all MATS measurements done at MATS_times (epoch)
MATS_df=dftop
MATS_df['EXPDate'] = MATS_df.index

#MALi = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
MALi = np.zeros([len(MATS_df), len(common_heights),3]) # if no bgr rem
Ma2 = np.zeros([len(MATS_df), len(common_heights),3]) # if no bgr rem
Ma3 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 
Ma4 = np.zeros([len(MATS_df), len(common_heights)-10+endcut,3]) 


for i in range(0,len(MATS_df)):
    #heights, profile = prepare_profile(MATS_df.iloc[i])

    # nearest time
    ir3_p=ir3.iloc[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')
    ir4_p=ir4.iloc[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')].drop('index')

    # set EXPDate accordingly
    ir3_p['EXPDate'] = ir3.index[ir3.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
    ir4_p['EXPDate'] = ir4.index[ir4.index.get_loc(MATS_df.iloc[i].EXPDate, method='nearest')]
    
    k = 0
    for offs in [0]:
        zsi,MALi[i,:,k],Ma3[i,:,k],Ma4[i,:,k] = prepare_measurment(MATS_df.iloc[i], ir3_p, ir4_p,subtract=False,endcut=endcut,offset=offs)
        k=k+1


 # %%
import math

def solar_zenith_angle_sunset(h, include_refraction=False):
    # Constants
    R_E = 6371  # Radius of the Earth in kilometers
    refraction_correction = 0.5  # Refraction correction in degrees
    
    # Solar zenith angle at sunset without refraction
    zenith_angle = math.degrees(math.acos(R_E / (R_E + h)))
    
    # Apply refraction correction if needed
    if include_refraction:
        zenith_angle -= refraction_correction
    
    return zenith_angle

# Example usage
altitude = 0  # Altitude in kilometers (e.g., 1 km)
sza = solar_zenith_angle_sunset(altitude)
#print(f"Solar zenith angle at sunset for {altitude} km altitude: {90-sza:.2f} degrees")

#%%
fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=4,abc='a.')
ais=[10,15,30,80,110,140,180,190]
#solz=np.linspace(95,70)

j = 0
for ai in ais:
    ai0,ai1=ai,ai+5

    #print(f'AI: {ds.altitude[ai0]} to {ds.altitude[ai1]}')

    x = MATS_df.TPsza.values
    y = np.polyfit(x, np.mean(MALi[:,ai0:ai1,0],axis=1), 3)
    y = np.polyval(y, x)
    xs, ys = zip(*sorted(zip(x, y)))

    x1 = ds.sza
    y1 = np.polyfit(x1, np.mean(ds.Li[:,ai0:ai1],axis=1), 3)
    y1 = np.polyval(y1, x1)
    xs1, ys1 = zip(*sorted(zip(x1, y1)))


    for i in range(0,len(MATS_df)):
        axs[j].scatter(MATS_df.TPsza[i],np.mean(MALi[i,ai0:ai1,0]),color='black',alpha=0.5,s=2)

    #axs[j].plot(xs, ys, color='black', label='MATS')
    #axs[j].plot(xs1, ys1, color='red', label='ODIN')

    #axs[j].plot(xs, ys, color='black', label='MATS')

    for i in range(0,len(ds.time)):
        axs[j].scatter(ds.sza[i],np.mean(ds.Li[i,ai0:ai1]),color='red',alpha=0.5,s=2)

    axs[j].axvline(90+solar_zenith_angle_sunset((ai0+ai1)/2), color='black', linestyle='--',alpha=0.75,linewidth=0.5)
    axs[j].legend()
    axs[j].format(xlabel='SZA',ylabel='Intensity [R]',title=f'alts: {ds.altitude[ai0].values} to {ds.altitude[ai1].values}')

    axs[j].format(ylim=(0, 9e14))

    j = j + 1

fig.suptitle('ODIN vs MATS')
# %%
fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=2,abc='a.')
ais=[20,80,110,140]
j=0
for ai in ais:
    ai0,ai1=ai,ai+5

    x = MATS_df.TPsza.values
    y = np.polyfit(x, np.mean(MALi[:,ai0:ai1,0],axis=1), 3)
    y = np.polyval(y, x)
    xs, ys = zip(*sorted(zip(x, y)))

    x1 = ds.sza
    y1 = np.polyfit(x1, np.mean(ds.Li[:,ai0:ai1],axis=1), 3)
    y1 = np.polyval(y1, x)
    xs1, ys1 = zip(*sorted(zip(x, y1)))

    axs[j].plot(xs1, np.array(ys1)/np.array(ys))
    axs[j].axhline(1, color='black', linestyle='--',alpha=0.75,linewidth=0.5)
    axs[j].axhline(0.8, color='red', linestyle='--',alpha=0.75,linewidth=0.5)

    axs[j].legend()
    axs[j].format(xlabel='SZA',ylabel='OSIRIS / MATS',title=f'alts: {ds.altitude[ai0].values} to {ds.altitude[ai1].values}')

    axs[j].format(ylim=(0.5, 1.5))
    j = j + 1
fig.suptitle('ODIN vs MATS')
# %%
fig, axs = pplt.subplots(figwidth='15cm',ncols=1, nrows=1,abc='a.')

j=0

ais = np.arange(10,190,10)
#ais = np.array([10,20,50,80,110,140,180,190])
#ais=[20,70,120,170]

for ai in ais:
    ai0,ai1=ai,ai+5

    x = MATS_df.TPsza.values
    y = np.polyfit(x, np.mean(MALi[:,ai0:ai1,0],axis=1), 3)
    y = np.polyval(y, x)
    xs, ys = zip(*sorted(zip(x, y)))

    x1 = ds.sza
    y1 = np.polyfit(x1, np.mean(ds.Li[:,ai0:ai1],axis=1), 3)
    y1 = np.polyval(y1, x)
    xs1, ys1 = zip(*sorted(zip(x, y1)))

    alt = (ds.altitude[ai0].values + ds.altitude[ai1].values) / 2

    axs[0].plot(xs1, np.array(ys1)/np.array(ys),color=((ais[-1]-ai)/ais[-1], 0, ai/ais[-1]), label=f'{alt} km')
    axs[0].axhline(1, color='black', linestyle='--',alpha=0.75,linewidth=0.5)
    axs[0].axhline(0.8, color='black', linestyle='--',alpha=0.75,linewidth=0.5)
    j = j + 1

axs[0].legend()
axs[0].format(xlabel='SZA',ylabel='OSIRIS / MATS',title=f'AI: {ds.altitude[ai0].values} to {ds.altitude[ai1].values}')

axs[0].format(ylim=(0.5, 1.5))
    
fig.suptitle('ODIN / MATS')
# %%
import math

def solar_zenith_angle_sunset(h, include_refraction=False):
    # Constants
    R_E = 6371  # Radius of the Earth in kilometers
    refraction_correction = 0.5  # Refraction correction in degrees
    
    # Solar zenith angle at sunset without refraction
    zenith_angle = math.degrees(math.acos(R_E / (R_E + h)))
    
    # Apply refraction correction if needed
    if include_refraction:
        zenith_angle -= refraction_correction
    
    return zenith_angle

# Example usage
altitude = 0  # Altitude in kilometers (e.g., 1 km)
sza = solar_zenith_angle_sunset(altitude)
print(f"Solar zenith angle at sunset for {altitude} km altitude: {90-sza:.2f} degrees")

# %%
