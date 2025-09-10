#%% This is READ_ODIN to line 276, for new data version download simply
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
    #profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e12)
    profile = np.array(image[:, col]*1e12)


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
    lat=np.array(data['lat'])
    lon=np.array(data['lon'])
    time=np.array(data['t'])

    flat_list = []
    for xs in time:
        for x in xs:
            flat_list.append(x)
    timeO=flat_list

    flat_list = []
    for xs in lat:
        for x in xs:
            flat_list.append(x)
    lat=flat_list

    flat_list = []
    for xs in lon:
        for x in xs:
            flat_list.append(x)
    lon=flat_list

    timeO = np.array(timeO)
    lon = np.array(lon)
    lat = np.array(lat)

else:
    lat=np.array(data['sm']['lat90'])
    lon=np.array(data['sm']['lon90'])
    #sza=data['sm']['sza']
    timeO=np.array(data['sm']['time90'])



#    for i in range(0,len(sza)):
#        if type(sza[i]) is list:
#            sza[i] = np.nan 

#    sza = np.array(sza)

    time = np.array(timeO)

jd = timeO + 2400000.5
t = Time(jd, format='jd')
timeO = t.to_datetime()


# SZA at 90
sza = np.zeros(len(timeO))
for i in range(0,len(sza)):
    sza[i]=pysolar.solar.get_altitude(lat[i],
                                      lon[i],
                                      pd.Timestamp(timeO[i],
                                                   tzinfo=timezone.utc).to_pydatetime())
    sza[i] = 90-sza[i]
# xarray
if NLC:
    ds_day_all = xr.Dataset(
        data_vars=dict(
            lat=(["time"], lat),
            lon=(["time"], lon)
        ),
        coords=dict(
            time=timeO,
        ),
        attrs=dict(description="ODIN"),
    )

else:
    ds = xr.Dataset(
        data_vars=dict(
            lat=(["time"], lat),
            lon=(["time"], lon),
            sza=(["time"], sza)
        ),
        coords=dict(
            time=timeO,
        ),
        attrs=dict(description="ODIN"),
    )
    
    if NIGHTGLOW:
        ds_night_all=ds.where(ds.sza > 100)
        ds_day_all=ds.where(ds.sza > 100)
    if DAYGLOW:
        ds_night_all=ds.where(ds.sza > 100,drop=True)
        ds_day_all=ds.where(ds.sza < 90,drop=True)

if NLC:
    ds_day_all=ds_day_all.sel(time = slice("2022-02-01", "2023-04-01")).dropna(dim='time')

else:
    ds_night_all=ds_night_all.sel(time = slice("2023-01-01", "2023-03-01")).dropna(dim='time')
    ds_day_all=ds_day_all.sel(time = slice("2023-01-01", "2023-03-01")).dropna(dim='time')


#%%
## COVERAGE MAP
"""
fig, axs = pplt.subplots(figwidth='15cm',ncols=1, nrows=1,sharex=0,proj='cyl')

axs.scatter(ds_day_all.lon,ds_day_all.lat, c='red', s=1,zorder=4)

axs.format(coast=True,land=True,landcolor='white',landzorder=2,facecolor='cloudy blue',coastcolor='black')
axs.format(latlines=30, lonlines=30,labels=True)
if not NLC:
    axs.scatter(ds_night_all.lon,ds_night_all.lat, c='blue', s=2, zorder=4)
    axs.format(title='OSIRIS MLT coverage (December 2022 - February 2023)')
    lgd=axs.legend([r'SZA $< 90^\circ$',r'SZA $> 100^\circ$'])
    lgd.legendHandles[0]._sizes = [20]
    lgd.legendHandles[1]._sizes = [20]
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/dayglow_nightglow_map.png',format='png')
else:
    axs.format(title='OSIRIS 22/23 NLC coverage')
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/nlc_map.png',format='png')
"""
#%%
"""
### NLC DAYGLOW AIRGLOW COVERAGE
fig, axs = pplt.subplots(figwidth='20cm',ncols=3, nrows=1,sharex=0)



if not NLC:
    #m=axs[0].scatter(ds.time,ds.lat, c='red', s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    #fig.colorbar(m,loc='r',title='longitude')
    axs[0].scatter(ds_night_all.time,ds_night_all.lat, c='blue', s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[0].scatter(ds_day_all.time,ds_day_all.lat, c='red', s=6,alpha=1,cmap='viridis',vmax=360,vmin=0)
    
    axs[1].scatter(ds_night_all.time,ds_night_all.lat, c='blue', s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[2].scatter(ds_day_all.time,ds_day_all.lat, c='red', s=6,alpha=1,cmap='viridis',vmax=360,vmin=0)
    axs[0].format(title='Airglow')
    axs[1].format(title='Nightglow SZA > 98')
    axs[2].format(title='Dayglow (SZA < 90)')
    fig.format(suptitle='OSIRIS MLT airglow coverage')
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/glow_emissions.png',format='png')

else:
    axs[0].format(title='NLC coverage')
    axs[0].scatter(ds_day_all.time,ds_day_all.lat, c=ds_day_all.lon, s=1,alpha=1,cmap='viridis',vmax=360,vmin=0)
    fig.savefig('/home/waves/projects/MATS/MATS-analysis/Bjorn/ODIN_MATS/figures/NLC.png',format='png')
"""
# %% 
# example of a single pass
#starttime=datetime(2023,2,1,0,0)
#stoptime=datetime(2023,3,1,0,0)

for ddd in [1]:
    starttime=datetime(2022,12,20,0,0)
    stoptime=datetime(2023,1,1,0,0)
    l1b_version="0.9"

    download=True
    run=True

    if download:
    # some filter settings - dayglow (ALL SZA)
        if NLC:
            tplat0, tplat1 = -90, -40
            ccdsel0, ccdsel1 = 6, 6
            filter={'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
            dftop_full=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

        elif NIGHTGLOW:
            dmin,dmax = 90, 180
            tplat0, tplat1 = 70, 90
            ccdsel0, ccdsel1 = 1, 5
            filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
            dftop_full=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

        elif DAYGLOW:
            dmin,dmax = 0, 95
            tplat0, tplat1 = -90, -70
            ccdsel0, ccdsel1 = 1, 5
            filter={'TPsza': [dmin, dmax], 'TPlat': [tplat0, tplat1], 'CCDSEL': [ccdsel0,ccdsel1]}
            dftop_full=read_MATS_data(starttime,stoptime,level="1b",version=l1b_version, filter=filter)

        dftop_full['EXPDate'] = pd.to_datetime(dftop_full['EXPDate'])

        dftop_full.to_pickle(f'/media/waves/AVAGO/data/MATS/MATS_ODIN_DATA/{savestr}/v09/{str(starttime)}_{str(stoptime)}.pkl')
    else:
        dftop_full=pd.read_pickle(f'/media/waves/AVAGO/data/MATS/MATS_ODIN_DATA/{savestr}/v09/{str(starttime)}_{str(stoptime)}.pkl')
        dftop_full=dftop_full.drop(labels=['SID','schedule_description_short'],axis=1)
        dftop_full['EXPDate'] = pd.to_datetime(dftop_full['EXPDate'])

# %%
