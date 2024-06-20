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


#%%
result_1d = xr.open_dataset("/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/1D_full/w1_march.nc")
result_1d = xr.open_dataset('/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/1D_full/IR1_MarW1_l1bv0.5_VER.nc')
#%%
ir1band=result_1d
ir1band=ir1band.where((ir1band.latitude < 70) & (ir1band.latitude > -70),drop=True)
#ir1band=ir1band.where((ir1band.longitude < 220) & (ir1band.longitude > 200),drop=True)
#ir1band=ir1band.where((ir1band.latitude < -15) & (ir1band.latitude > -50))
#%% 

fig = pplt.figure(figwidth='12cm')
axs = fig.subplots(ncols=1, nrows=1, proj='cyl')
fig.format(suptitle='Peak emission (Full band) IR1 (~17:30 LT)',abc='a.',abcborder=False)
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

# plot max emission of profile
plotd=ir1band.ver.max(dim='z_r')
#plotd=ir1band.ver[:,10:-10].idxmax('z_r')
vmean=np.mean(plotd)
vstd=np.std(plotd)
vmax,vmin=1.4e5,8e4

m=axs.scatter(ir1band.longitude,ir1band.latitude,c=plotd ,s=0.2,cmap='Thermal',vmin=vmin,vmax=vmax)
#axs.format(latlim=(-10, 50))
fig.colorbar(m, label="1e4 ph/cm3/s",loc='r',length=0.8)

#%% 

fig = pplt.figure(figwidth='12cm')
axs = fig.subplots(ncols=1, nrows=1, proj='cyl')
fig.format(suptitle='Peak height IR1 (~17:30 LT)',abc='a.',abcborder=False)
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

# plot peak height
plotd=ir1band.ver[:,:].idxmax('z_r')
vmean=np.mean(plotd)
vstd=np.std(plotd)
vmax,vmin=88,77

m=axs.scatter(ir1band.longitude,ir1band.latitude,c=plotd ,s=0.2,cmap='Thermal',vmin=vmin,vmax=vmax)
fig.colorbar(m, label="km",loc='r',length=0.8)

#%% 

# plot min emission of profile

fig = pplt.figure(figwidth='12cm')
axs = fig.subplots(ncols=1, nrows=1, proj='cyl')
fig.format(suptitle='Minima (Full band) IR1 (~17:30 LT)',abc='a.',abcborder=False)
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

# plot max emission of profile
plotd=ir1band.ver.min(dim='z_r')
#plotd=ir1band.ver[:,10:-10].idxmax('z_r')
vmean=np.mean(plotd)
vstd=np.std(plotd)
vmax,vmin=3.5e4,1e4

m=axs.scatter(ir1band.longitude,ir1band.latitude,c=plotd ,s=0.2,cmap='Thermal',vmin=vmin,vmax=vmax)
fig.colorbar(m, label="ph/cm3/s",loc='r',length=0.8)


#%% 

fig = pplt.figure(figwidth='20cm')
axs = fig.subplots(ncols=2, nrows=1, proj='cyl')
fig.format(suptitle='IR1 (full band) - Mar 1st-7th')
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

alts=[85,95]

# plot a specific height
for i in range(0,len(alts)):
    plotd=ir1band.sel(z_r=alts[i]).ver
    #plotd=plotd.where(plotd['time.day'] < 25)
    #plotd=ir1band.ver[:,10:-10].idxmax('z_r')
    vmean=np.mean(plotd)
    vstd=np.std(plotd)
    vmax,vmin=vmean+2*vstd,vmean-2*vstd
    #vmax,vmin=1.2e5,5e4
    axs[i].format(title=f'{alts[i]} km')
    m=axs[i].scatter(ir1band.longitude,ir1band.latitude,c=plotd ,s=0.9,cmap='Thermal',vmin=vmin,vmax=vmax)
    axs[i].colorbar(m, label="VER [ph/cm3/s]",loc='r',length=0.8)
    axs[i].format(lonlabels='b', latlabels='l')
fig.save('MarW1_85_95.png',format='png')
#%% 

# plot peak height
#plotd=ir1band.sel(z_r=84).ver/1e6/5*4*np.pi
plotd=ir1band.ver[:,:].idxmax('z_r')
vmean=np.mean(plotd)
vstd=np.std(plotd)
vmax,vmin=87,70
plt.figure(figsize=(8,3))
plt.scatter(ir1band.longitude,ir1band.latitude, c=plotd, s=0.15,vmin=vmin,vmax=vmax)
plt.colorbar()
#%%

plt.figure(figsize=(4,4))
(ir1band.ver[30:-30:1,:]).plot.line(y='z_r',add_legend=False,xlim=([0,2e5]))
#if xa is not None:
#    plt.plot(xa/1e6,retrival_heights, linestyle='--', color='black')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#plt.title(f'({channel}) A-band (full band)')
plt.grid()
plt.tight_layout()
#plt.savefig(f'{channel}VERprofiles_sub.png',format='png')

plt.figure(figsize=(12,4))
(ir1band.ver[30:-30:1,:]).plot.pcolormesh(y='z_r',vmin=0,vmax=2e5,ylim=([60,110]))
plt.title('A-band intensity (full band) photons cm-3 s-1')
plt.tight_layout()
#plt.savefig(f'{channel}.png',format='png')

plt.figure(figsize=(12,4))
#ch.profile[:,:].plot.pcolormesh(y='z',vmin=0,ylim=([60e3,110e3]))
plt.tight_layout()
#plt.savefig(f'{channel}LOS_sub.png',format='png')

#%%

channels=['IR1','IR2']
colors=['red','blue']
i=0
plt.figure(figsize=(4,4))
for channel in channels:

    result_1d = xr.open_dataset(f"{channel}_Mar1_l1bv0.5_VER.nc")
    ir1band=result_1d
    ir1band=ir1band.where((ir1band.latitude < 30) & (ir1band.latitude > -30),drop=True)


    #ir1band_eq=ir1band.where((ir1band.latitude < 20) & (ir1band.latitude > -20),drop=True)

    (ir1band.ver[30:-30:1,:]).plot.line(y='z_r',add_legend=False, color=colors[i], alpha=0.5, linewidth=0.1)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #plt.title(f'({channel}) A-band (full band)')
    plt.grid()
    #plt.legend()

    plt.title('20S to 20N; [1e4 ph/cm3/s]; IR1 red; IR2 blue')
    plt.tight_layout()
    i=i+1
#plt.savefig(f'{channel}VERprofiles_sub.png',format='png')

#%% 

from scipy.interpolate import CubicSpline
overlap=np.array([0.78584363, 0.7656091 , 0.74607243, 0.72726005, 0.70918753,
       0.69185982, 0.67527261, 0.659414  , 0.64426607, 0.62980651,
       0.61600987, 0.60284864, 0.59029418, 0.57831736, 0.56688911,
       0.55598082, 0.54556466, 0.53561371, 0.52610218, 0.51700545,
       0.50830014])
temps=np.array([100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
       230, 240, 250, 260, 270, 280, 290, 300])
tempspl=CubicSpline(overlap[-1::-1],temps[-1::-1])

#%%
plt.figure(figsize=(5,5))
channels=['IR1','IR2']
colors=['red','blue']
i=0

lat0,lat1=-5,5

result_IR1 = xr.open_dataset(f"IR1_Mar1_l1bv0.5_VER.nc")
result_IR2 = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER.nc")

result_IR1=result_IR1.where((result_IR1.latitude < lat1) & (result_IR1.latitude > lat0),drop=True)
result_IR2=result_IR2.where((result_IR2.latitude < lat1) & (result_IR2.latitude > lat0),drop=True)

temps1=tempspl((result_IR1.ver[:,:].values*0.6/result_IR2.ver.values).T)
#plt.plot(temps1,result_IR1.z_r,color='red',linewidth=0.7, alpha=0.5,label='mar')

# another abs
result_IR1 = xr.open_dataset(f"IR1_Mar1_l1bv0.5_VER_sepabs.nc")
result_IR2 = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER_sepabs.nc")

result_IR1=result_IR1.where((result_IR1.latitude < lat1) & (result_IR1.latitude > lat0),drop=True)
result_IR2=result_IR2.where((result_IR2.latitude < lat1) & (result_IR2.latitude > lat0),drop=True)

temps2=tempspl((result_IR1.ver[:,:].values*0.6/result_IR2.ver.values).T)
#plt.plot(temps2,result_IR1.z_r,color='blue',linewidth=0.7, alpha=0.5,label='sep')

plt.plot(temps1,result_IR1.z_r,color='blue',linewidth=0.7, alpha=0.5,label='sep')
plt.plot(temps2,result_IR1.z_r,color='red',linewidth=0.7, alpha=0.5,label='sep')



plt.title(f'lat0: {lat0}, lat1: {lat1}')
#plt.legend()
plt.xlabel('temperature [k]')
#plt.xlim([180,300])
plt.ylabel('[km]')
#"(ir1band.ver[30:-30:1,:]).plot.line(y='z_r',add_legend=False, color=colors[i], alpha=0.5, linewidth=0.1)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#plt.title(f'({channel}) A-band (full band)')
#plt.grid()
#plt.legend()

#plt.title('20S to 20N; [1e4 ph/cm3/s]; IR1 red; IR2 blue')
#plt.tight_layout()
i=i+1
plt.savefig(f'temps_lat0_{lat0}_lat1{lat1}_abscomp.png',format='png')


#%

plt.figure()
plt.plot(result_IR1.longitude,result_IR1.latitude)
plt.xlabel('longitude')
plt.ylabel('latitude')
# %%

#%%
fig, axs = pplt.subplots(figwidth='12cm',ncols=2, nrows=2,abc='a.',sharex=0)
channels=['IR1','IR2']
colors=['red','blue']
i=0

lat0,lat1=-5,5

result_IR1 = xr.open_dataset(f"IR1_Mar1_l1bv0.5_VER.nc")
result_IR1_sep = xr.open_dataset(f"IR1_Mar1_l1bv0.5_VER_sepabs.nc")


result_IR1=result_IR1.where((result_IR1.latitude < lat1) & (result_IR1.latitude > lat0),drop=True)
result_IR1_sep=result_IR1_sep.where((result_IR1_sep.latitude < lat1) & (result_IR1_sep.latitude > lat0),drop=True)

result_IR2_sep = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER_sepabs.nc")
result_IR2 = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER.nc")

result_IR2=result_IR2.where((result_IR2.latitude < lat1) & (result_IR2.latitude > lat0),drop=True)
result_IR2_sep=result_IR2_sep.where((result_IR2_sep.latitude < lat1) & (result_IR2_sep.latitude > lat0),drop=True)


# TEMPS
for i in range(0,len(result_IR1.time)):
    axs[0].plot(result_IR1.isel(time=i).ver,result_IR1.z_r,color='blue',linewidth=0.7, alpha=0.5)
    axs[0].plot(result_IR1_sep.isel(time=i).ver,result_IR1_sep.z_r,color='red',linewidth=0.7, alpha=0.5)


    axs[1].plot(result_IR2.isel(time=i).ver,result_IR2.z_r,color='blue',linewidth=0.7, alpha=0.5)
    axs[1].plot(result_IR2_sep.isel(time=i).ver,result_IR2_sep.z_r,color='red',linewidth=0.7, alpha=0.5)


channels=['IR1','IR2']
colors=['red','blue']
i=0

lat0,lat1=-9,9

result_IR1 = xr.open_dmataset(f"IR1_Mar1_l1bv0.5_VER.nc")
result_IR2 = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER.nc")

result_IR1=result_IR1.where((result_IR1.latitude < lat1) & (result_IR1.latitude > lat0),drop=True)
result_IR2=result_IR2.where((result_IR2.latitude < lat1) & (result_IR2.latitude > lat0),drop=True)

tempsmar=tempspl((result_IR1.ver[:,:].values*0.6/result_IR2.ver.values).T)

# another abs
result_IR1 = xr.open_dataset(f"IR1_Mar1_l1bv0.5_VER_sepabs.nc")
result_IR2 = xr.open_dataset(f"IR2_Mar1_l1bv0.5_VER_sepabs.nc")

result_IR1=result_IR1.where((result_IR1.latitude < lat1) & (result_IR1.latitude > lat0),drop=True)
result_IR2=result_IR2.where((result_IR2.latitude < lat1) & (result_IR2.latitude > lat0),drop=True)

tempssep=tempspl((result_IR1.ver[:,:].values*0.6/result_IR2.ver.values).T)


for i in range(0,len(result_IR1.time)):
    axs[2].plot((tempsmar[:,i]),result_IR1.z_r,color='blue',linewidth=0.7, alpha=0.5)
    axs[2].plot((tempssep[:,i]),result_IR1.z_r,color='red',linewidth=0.7, alpha=0.5)


#axs[0].format(xlim=[5e4,15e4])
axs[2].format(xlim=[150,500], title='retrieved TEMP')


axs[0].format(title='retrieved VER (IR1)')
axs[1].format(title='retrieved VER (IR2)')
# another abs



### ABSORPTION STUFF

axs[3].plot(TMAR,z,color='blue')
axs[3].plot(TSEP,z,color='red')

axs[3].format(ylim=([50,105]),xlim=([150,350]),title='MSIS abs TEMP')

axs[2].format(ylim=([60,105]),xlim=([150,350]))
axs[0].format(ylim=([60,105]))

fig.savefig('changes_absorption.png',format='png')
# %%


import ref_index as RF
from os.path import exists
from numba import jit
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline
from hapi import *
import requests as R
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
#


@jit(nopython=True)
def gauss(x, fwhm):
    sigma = fwhm/2./np.sqrt(2*np.log(2.))
    return (np.exp(-(x*x/2/sigma/sigma))/sigma/np.sqrt(2*np.pi))


Kb = 1.380649e-23  # J/K Boltzmans constant
h = 6.626176e-34  # Plancks constant
AMU = 1.6603145e-27
C = 299792458.  # m/s. Speed of light in vacuo

db_begin('Abanddata')
if not exists('./Abanddata/oxygen2.header'):

    fetch_by_ids('oxygen', [36], 1/776e-7, 1/749e-7)
c2 = 1.4387770  # CM K

#url = "http://129.16.35.2:8080/msis/2019-03-17T12:00:00/0/130/150" # update for MAr
url = "http://129.16.35.2:8080/msis/2019-09-17T17:30:00/0/130/150" # completely different month


atm = R.get(url).json()
N = np.array(atm['n2'])+np.array(atm['o2'])+np.array(atm['o'])
pres = N*Kb*atm['T']
atm['p'] = pres
z = np.array(atm['z'])
o2 = np.array(atm['o2'])/1e6  # to cm-3
TSEP = atm['T']


Kb = 1.380649e-23  # J/K Boltzmans constant
h = 6.626176e-34  # Plancks constant
AMU = 1.6603145e-27
C = 299792458.  # m/s. Speed of light in vacuo

db_begin('Abanddata')
if not exists('./Abanddata/oxygen2.header'):

    fetch_by_ids('oxygen', [36], 1/776e-7, 1/749e-7)
c2 = 1.4387770  # CM K

url = "http://129.16.35.2:8080/msis/2019-03-17T12:00:00/0/130/150" # update for MAr
#url = "http://129.16.35.2:8080/msis/2019-09-17T17:30:00/0/130/150" # completely different month


atm = R.get(url).json()
N = np.array(atm['n2'])+np.array(atm['o2'])+np.array(atm['o'])
pres = N*Kb*atm['T']
atm['p'] = pres
z = np.array(atm['z'])
o2 = np.array(atm['o2'])/1e6  # to cm-3
TMAR = atm['T']


# %%

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

# HOVMÖLLER N AND S (60 to 70)
monthstrings=['w1_march','w2_march']
monthstrings=['w1_february','w2_february']

#monthstrings=['032023']
asc_desc=''

# lats figure 1
lat11, lat12 = (60,70)
# lats figure 2
lat21, lat22 = (-60, -55)

# height
heightind=6

# common scale
common=False

y_bins=np.arange(0,380,20)
days=np.arange(12,28)

remove_zm=False
x_bins = np.arange(-70,70,0.5)
pplt.rc.axesfacecolor = 'gray4'
pplt.rc['grid.color'] ='white'
pplt.rc['grid.alpha'] = 0.5 

fig = pplt.figure(figwidth='8.3cm',figheight='10cm',suptitle='Dayglow (full band) VER (February)')
ax = fig.subplots(ncols=2, nrows=1,xlabel='Longitude [deg]',abc='a.')
height=80

for k in range(1,2):
    loadbool = 1
    meanz=[]
    timez=[]


    CCDmeans=xr.open_dataset(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/1D_full/{monthstrings[0]}.nc')
    CCDmeans=CCDmeans.sel(z_r=height)

    #if remove_zm:
    #    CCDmeans=CCDmeans.groupby_bins("latitude", x_bins) - CCDmeans.groupby_bins("latitude", x_bins).mean(dim='time').rolling(latitude_bins=2).mean()

    for i in range(0,len(days)):

        if (days[i] > 14) and (loadbool == 1):
            CCDmeans=xr.open_dataset(f'/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/1D_full/{monthstrings[1]}.nc')
            CCDmeans=CCDmeans.sel(z_r=height)
            loadbool = 0

        if days[i] != 13:
            y=CCDmeans.where(CCDmeans['time.day'] == days[i], drop=True)

            if (len(y.longitude) > 0):

                timez.append(f"{int(y['time.day'][0])}/{int(y['time.month'][0])}")

                if k == 0:
                    means=y.where((y.latitude<lat12) & (y.latitude>lat11)).groupby_bins("longitude", y_bins).mean(dim='time')
                    y=y.where((y.latitude<lat12) & (y.latitude>lat11))

                if k == 1:

                    means=y.where((y.latitude>lat21) & (y.latitude<lat22)).groupby_bins("longitude", y_bins).mean(dim='time')
                    y=y.where((y.latitude>lat21) & (y.latitude<lat22))

                meanz.append(means.ver.values)

            #x = y.longitude
            #y = y.ver

    #for i in range(0,len(meanz)):
    #    for j in range(0,len(meanz[i])):
    #        if meanz[i][j] < 0:
    #            meanz[i][j] = np.nan
    """
    
    if month=='022023':
        meanz = meanz[7:-1]
        timez=timez[7:-1]

    if common:
        if k == 0:
            mean_cl=np.nanmean(np.array(meanz))
            std_cl=np.nanstd(np.array(meanz))
    else:
    """
    meanz=np.array(meanz)*4*np.pi/1e6/5
    if k == 1:
        

        mean_cl=np.nanmean(np.array(meanz))
        std_cl=np.nanstd(np.array(meanz))


    m=ax[k].contourf(y_bins[1:]-10,np.arange(0,len(meanz)),np.array(meanz),vmin=mean_cl-1.5*std_cl,vmax=mean_cl+1.5*std_cl,cmap='BR',levels=15,extend='both')


    ax[k].format(yticks=pplt.arange(0, len(meanz),3),yticklabels=timez[0:-1:3],grid=True)

ax[0].set_ylabel('Day of measurement [dd/mm]')

ax[0].set_title(f'{lat11}N-{lat12}N')    
ax[1].set_title(f'{-lat22}S-{-lat21}S')
fig.colorbar(m,loc='b',title='VER [ph/s/cm3]')


fig.save('hovemoller.png', format='png')

#%%
import xarray as xr
CCDmeans=xr.open_dataset('/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/1D_full/w1_march.nc')
# %%

# DEBUG IMAGES
from mats_utils.rawdata.read_data import read_MATS_data
starttime=datetime(2023,3,3,0,0)
stoptime=datetime(2023,3,9,23,59) 
dmin,dmax = 0, 95

dftop=read_MATS_data(starttime,stoptime,level="1b",version='0.5', filter={'TPsza': [dmin, dmax], 'TPlat': [-70, 70]})
#%% ISOLATE SMALL REGION
dftop=dftop.where((dftop.TPlat < 22) & (dftop.TPlat > 21)).dropna()

# %%
ascending=True
if ascending:
    midday=DT.time(12, 0, 0)
    dftop['TPlocaltime'] = pd.to_datetime(dftop['TPlocaltime'])
    dftop = dftop[(dftop['TPlocaltime'].dt.time > midday)]
# %%
ir1 = dftop[dftop['channel'] == 'IR1'].dropna().reset_index()#[0:10]
ir3 = dftop[dftop['channel'] == 'IR3'].dropna().reset_index()#[0:10]
ir4 = dftop[dftop['channel'] == 'IR4'].dropna().reset_index()#[0:10]
# %%
fig, axs = pplt.subplots(figwidth='15cm',ncols=2, nrows=2,abc='a.',sharex=0)
fig.format(suptitle='Image means',ylabel='counts')

axs[0].scatter(ir1.EXPDate.values,np.mean(np.stack(ir1.ImageCalibrated)[:,:,:], axis=(1,2)),s=0.4, c=np.stack(ir1.afsTangentH_wgs84),cmap='viridis')
axs[0].format(title='IR1')

axs[2].scatter(ir3.EXPDate.values,np.mean(np.stack(ir3.ImageCalibrated)[:,:,:], axis=(1,2)),s=0.4, c=np.stack(ir3.afsTangentH_wgs84),cmap='viridis')
axs[2].format(title='IR3')

m=axs[3].scatter(ir4.EXPDate.values,np.mean(np.stack(ir4.ImageCalibrated)[:,:,:], axis=(1,2)),s=0.4, c=np.stack(ir4.afsTangentH_wgs84),cmap='viridis')
axs[3].format(title='IR4')

fig.colorbar(m, loc='b',label='afsTangentH_wgs84')
#fig.savefig('limb_radiances.png',format='png')
#plt.scatter(ir3.TPsza,np.mean(np.stack(ir3.ImageCalibrated)[:,:,:], axis=(1,2)), s=0.05)
# %%
fig, axs = pplt.subplots(figwidth='20cm',ncols=2, nrows=2,abc='a.',sharex=0)
fig.format(suptitle='Image means',ylabel='counts')
axs[0].scatter(ir3.EXPDate,np.mean(np.stack(ir3.ImageCalibrated)[:,:,:], axis=(1,2)),s=0.4, c=np.stack(ir3.afsTangentH_wgs84),cmap='viridis')
axs[0].format(title='IR3')
# %%
