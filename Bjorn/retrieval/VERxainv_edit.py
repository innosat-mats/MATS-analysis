#%%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
#import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from numpy.linalg import inv

import sys
sys.path.append('/home/waves/projects/MATS/MATS-analysis/Donal/retrievals')
from oem_functions import linear_oem

# %%
# ir1=xr.load_dataset('/Users/donal/projekt/IR1vertest.nc')
# ir2=xr.load_dataset('/Users/donal/projekt/IR2vertest.nc')
# ir3=xr.load_dataset('/Users/donal/projekt/IR3vertest.nc')
# ir4=xr.load_dataset('/Users/donal/projekt/IR4vertest.nc')
# ir1=xr.load_dataset('/Users/donal/projekt/SIW/IR1Decvertest.nc')
# ir2=xr.load_dataset('/Users/donal/projekt/SIW/IR2Decvertest.nc')
# ir3=xr.load_dataset('/Users/donal/projekt/SIW/IR3Decvertest.nc')
# ir4=xr.load_dataset('/Users/donal/projekt/SIW/IR4Decvertest.nc')
ir1=xr.load_dataset('/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/IR1Feb17vertest_1d.nc')
#ir2=xr.load_dataset('/home/waves/projects/MATS/MATS-analysis/Bjorn/retrieval/IR2Feb17vertest_1d.nc')
#ir3=xr.load_dataset('/Users/donal/projekt/SIW/IR3Mar29vertest.nc')
#ir4=xr.load_dataset('/Users/donal/projekt/SIW/IR4Mar29vertest.nc')
# ir1=xr.load_dataset('/Users/donal/projekt/SIW/IR1Feb11vertest.nc')
# ir2=xr.load_dataset('/Users/donal/projekt/SIW/IR2Feb11vertest.nc')
# ir3=xr.load_dataset('/Users/donal/projekt/SIW/IR3Feb11vertest.nc')
# ir4=xr.load_dataset('/Users/donal/projekt/SIW/IR4Feb11vertest.nc')
# %%
def make_k(ch,i,retrival_heights):
    
    s_140=1700e3
    steps=10 #m steps
    s_steps = np.arange(s_140,s_140 + 2e6,steps) 

    k=np.zeros((len(ch.z),len(retrival_heights)-1))
    d = pd.Timestamp(ch['time'][i].values).replace(tzinfo=sfapi.utc)
    ts = sfapi.load.timescale()
    t=ts.from_datetime(d)
    ecipos=ch['ecipos'][i]
    localR = np.linalg.norm(sfapi.wgs84.latlon(ch.TPlat[i], ch.TPlon[i], elevation_m=0).at(t).position.m)
    for ki,ecivec in enumerate(ch.ecivec[i].values) :
        pos=np.expand_dims(ecipos, axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
        heights=np.linalg.norm(pos,axis=0)-localR
        counts,bins=np.histogram(heights/1000,retrival_heights)
        k[ki,:]=counts*steps
    return k


 # %%
ch=ir1
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
retrival_heights= np.arange(60,111,1)


for n in range(len(ch.time)): 
    profile=ch.profile[n].values
    Se=np.diag(profile)
    #k=make_k(ch,n,retrival_heights)
    k=ch.k[n].values
    xa=np.ones_like(retrival_heights)
    Sa=np.diag(xa) * np.max(profile)
    xa=0*xa
    
    x0, A, Ss, Sm = linear_oem(k, Se, Sa, profile, xa)
    ver.append(x0)
    mr.append(A.sum(axis=1)) #sum over rows 
    error2_retrieval.append(np.diag(Sm))
    error2_smoothing.append(np.diag(Ss))

"""
#### MEAN VER AS A-PRIORI

xa = np.mean(ver, axis=0)
xa[-10:] = 1
xa[0:50] = 1

ch=ir1
mr = []
error2_retrieval = []
error2_smoothing = []
ver = []
limb_fit = []
time_save = []
retrival_heights= np.arange(60,111,1)
Sa=0
k=0
Se=0

for n in range(len(ch.time)): 
    profile=ch.profile[n].values
    Se=np.diag(profile)
    #k=make_k(ch,n,retrival_heights)
    k=ch.k[n].values
    Sa=np.diag(xa) * np.max(profile)

    x0, A, Ss, Sm = linear_oem(k, Se, Sa, profile, xa)
    ver.append(x0)
    mr.append(A.sum(axis=1)) #sum over rows 
    error2_retrieval.append(np.diag(Sm))
    error2_smoothing.append(np.diag(Ss))
"""
# %%
result_1d = xr.Dataset().update({
        'time': (['time'], ch.time.values),
        'z_r': (['z_r',], ch.z_r.values, {'units': 'km'}),
        'ver': (['time','z_r'], ver, {'long name': 'VER', 'units': 'photons cm-3 s-1'}),
        'mr': (['time','z_r'], mr),
        'error2_retrieval': (['time','z_r'], error2_retrieval),
        'error2_smoothing': (['time','z_r'], error2_smoothing),
        #'limb_fit': (['time','pixel'], limb_fit),
        'latitude': (['time',], ch.TPlat.values),
        'longitude': (['time',], ch.TPlon.values),
        'channel': (['time',], ch.channel.values),
        })
ir1band=result_1d
# %%
overlap=np.array([0.78584363, 0.7656091 , 0.74607243, 0.72726005, 0.70918753,
       0.69185982, 0.67527261, 0.659414  , 0.64426607, 0.62980651,
       0.61600987, 0.60284864, 0.59029418, 0.57831736, 0.56688911,
       0.55598082, 0.54556466, 0.53561371, 0.52610218, 0.51700545,
       0.50830014])
temps=np.array([100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
       230, 240, 250, 260, 270, 280, 290, 300])
tempspl=CubicSpline(overlap[-1::-1],temps[-1::-1])
# %%

import matplotlib.pyplot as plt
plt.figure(figsize=(12,6))
ir1.profile.where(ir1.TPsza <98,drop=True).plot(y='z',vmin=0)
# %%
ir1.attrs["title"]="IR1"
# %%
plt.figure(figsize=(12,6))
(ir1.profile[30:230]/1e9*3.57/0.60).plot(y='z',vmin=0)
plt.figure(figsize=(12,6))
(ir1band.ver[30:230]/1e9*3.57/0.60).plot(y='z_r',robust=True,vmin=0)
plt.title('A-band intensity (full band) photons cm-3 s-1')
plt.savefig('Aband_intenstiy_fullband.png',format='png')
# %%
plt.figure(figsize=(12,6))
plt.pcolormesh(ir1band.where((ir1band.z_r > 70) & (ir1band.z_r < 100)).ver[30:230].T/1e9*3.57/0.60,vmin=250,vmax=1500)
#(result_1d.ver[30:50]/1e9*3.57/0.60).plot.line(y='z_r',add_legend=False);
plt.title('A-band intensity (full band) photons cm-3 s-1')
# %%
plt.figure()
ir1.profile[30:230].plot(y='z',vmin=0)
# %%
plt.figure()
ir1.TPlat[30:230].plot()
plt.figure()
ir1.TPlon[30:230].plot()
# %%
plt.figure()
smallband=(ir1band.ver[30:230]/1e9*3.57).values
wideband=(ir2band.ver[30:230]/1e9*8.16).values
lats=ir1band.latitude[30:230]
plt.contourf(lats,retrival_heights,(smallband/wideband).T,np.linspace(0.5,1.5,11) )
plt.clim([0.5,1.5])
plt.colorbar()
# %%
plt.figure()
smallband=(ir1band.ver[30:230]/1e9*3.57).values
wideband=(ir2band.ver[30:230]/1e9*8.16).values
lats=ir1band.latitude[30:230]
plt.contourf(lats,retrival_heights,tempspl(0.7*smallband/wideband).T,np.linspace(100,310,10) )
plt.clim([100,300])
plt.colorbar()
# %%
