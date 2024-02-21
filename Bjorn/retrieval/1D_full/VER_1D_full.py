# %%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
from mats_l2_processing.inverse_model import Sa_inv_tikhonov
#import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline,interp1d
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from numpy.linalg import inv
from fast_histogram import histogramdd
from bisect import bisect_left
import matplotlib.pyplot as plt
import sys
import scipy.sparse as sparse


# %% 

def linear_oem(K, Se, Sa_inv, y, xa):
    # Adapted from Donal (uses Sa_inv instead)
    if len(y.shape) == 1:
        y = y.reshape(len(y),1)
    if len(xa.shape) == 1:
        xa = xa.reshape(len(xa),1)
        
    #if len(y)<len(xa): # m form
    #    G = Sa.dot(K.T).dot(np.linalg.inv(K.dot(Sa).dot(K.T) + Se))
        
    #else: # n form
    Se_inv = np.linalg.inv(Se)
    Sa_inv = Sa_inv.toarray()
    #Sa_inv = np.linalg.inv(Sa)
    G = np.linalg.inv(K.T.dot(Se_inv).dot(K) + Sa_inv).dot(K.T).dot(Se_inv)
#        G= np.linalg.solve(K.T.dot(Se_inv).dot(K) + Sa_inv, (K.T).dot(Se_inv))
    
    x_hat = xa + G.dot(y - K.dot(xa)) 
    A = G.dot(K)
    #Ss = (A - np.identity(len(xa))).dot(Sa).dot((A - np.identity(len(xa))).T) # smoothing error
    Sm = G.dot(Se).dot(G.T) #retrieval noise 
    
    return x_hat.squeeze()

def prepare_profile(ch):
    # TBD comment
    # calculates mean profile and tan heights
    image = image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e15)
    return heights, profile

def generate_ks(df,retrival_heights, s_steps, splinedfactor):
    # TBD comment
    # calculates ks from mean profile
    ts = sfapi.load.timescale()
    profiles = []
    heights = []
    ecipos = []
    ecivecs = []
    ks=[]

    for i in range(len(df)):
        print(i)
        zs, p = prepare_profile(df.iloc[i])

        # remove some pixels top and bottom
        if df.iloc[i]['channel'] == 'IR1':
            zs = zs[10:-25]
            p = p[10:-25]
            nrows = len(zs)

        profiles.append(p)
        heights.append(zs)
        ecipos.append(df['afsGnssStateJ2000'][i][0:3])
        d = df['EXPDate'][i]
        t = ts.from_datetime(d)
        localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat[i], df.TPlon[i], elevation_m=0).at(t).position.m)
        q = df['afsAttitudeState'][i]
        
        quat = R.from_quat(np.roll(q, -1))
        ypixels = np.linspace(0, nrows, 5)
        x, yv = pix_deg(df.iloc[i], int(df['NCOL'][i]/2), ypixels)
        qp = R.from_quat(df['qprime'][i])
        ecivec = np.zeros((3, len(yv)))
        k=np.zeros((nrows,len(retrival_heights)))
        for irow, y in enumerate(yv):
            los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
            ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
        cs_eci=CubicSpline(ypixels,ecivec.T)
        for irow in range(nrows):
            ecivec=cs_eci(irow)
            pos=np.expand_dims(ecipos[-1], axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
            point_height=np.linalg.norm(pos,axis=0)-localR
            
            # absorption weight
            minarg=point_height[:].argmin() # pos of lowest tan height
            target_tangent=zs[irow]/1000
            distances=(s_steps-s_steps[minarg])/1000 #to km for intepolation
            lowertan=bisect_left(tanz,target_tangent)-1
            uppertan=lowertan+1
            lowerfactor=splinedfactor[lowertan](distances)
            upperfactor=splinedfactor[uppertan](distances)
            newfactor=interp1d([tanz[lowertan],tanz[uppertan]],np.array([lowerfactor,upperfactor]).T)
            weight=np.exp(newfactor(target_tangent))

            # calculate k
            counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)),weights=weight)

            k[irow,:]=counts*steps
            ecivecs.append(ecivec)

        ks.append(k)
    ecivecs= np.reshape(ecivecs,(len(df),-1,3))

    return ks, profiles, ecipos, ecivecs, heights

def xainvert(ch,retrival_heights, weight_0, weight_1):
    # TBD comment
    # invert using oem

    ver = []

    for n in range(len(ch.time)): 
        profile=ch.profile[n].values
        Se=np.diag(profile)
        k=ch.k[n].values
        xa=np.ones_like(retrival_heights)
        Sa_inv,_=Sa_inv_tikhonov(np.array([retrival_heights]), weight_0, [weight_1], volume_factors=False, store_terms=False)
        #Sa=np.diag(xa) * np.max(profile) / (1.3*10**11)

        xa=0*xa
        
        x0 = linear_oem(k, Se, Sa_inv, profile, xa)
        ver.append(x0)
        #mr.append(A.sum(axis=1)) #sum over rows 
        #error2_retrieval.append(np.diag(Sm))
        #error2_smoothing.append(np.diag(Ss))

    return ver

#%% 
# load images
channel='IR1'
starttime=datetime(2023,2,17,0,50)
stoptime=datetime(2023,2,17,1,20)

dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.6")
df = dftop[dftop['channel'] == channel].dropna().reset_index()#[0:10]

# absorption weights (TBD: call generation of these)
tanz, splinedfactor=np.load(f"splinedlogfactors{channel}_feb.npy",allow_pickle=True)

# select part of orbit
offset = 10
num_profiles = 250 #use 50 profiles for inversion
df = df.loc[offset:offset+num_profiles]
df = df.reset_index(drop=True)

# %% retrieval settings
retrival_heights= np.arange(60,100,1)
retrival_heights= np.arange(65,100,1)
s_140=1700e3
steps=100 #m steps
s_steps = np.arange(s_140,s_140 + 2e6,steps)

#%% Generate ks // save netcdf
ks, profiles, ecipos, ecivecs, heights = generate_ks(df,retrival_heights,
                                                     s_steps, splinedfactor)
z=np.array(heights).mean(axis=0)

ch = xr.Dataset({
    'time': (['time'], df.EXPDate),
    'channel': (['time'], df.channel),
    'satlat': (['time'], df.satlat, {'long_name': "MATS' Latitude", 'units': 'deg'}),
    'satlon': (['time'], df.satlon),
    'TPlat': (['time'], df.TPlat, {'long_name': "TP Latitude", 'units': 'deg'}),
    'TPlon': (['time'], df.TPlon),
    'TPsza': (['time'], df.TPsza),
    'TPssa': (['time'], df.TPssa),
    'z':  (['z'], z, {'long_name': 'Approx Altitude', 'units': 'm'}),
    'ecipos': (['time', 'xyz'], ecipos),
    'ecivec': (['time', 'z', 'xyz'], ecivecs),
    'profile': (['time', 'z'], profiles, {'long_name': 'LOS  intensity', 'units': 'Photons m-2 nm-1 sr-1 s-1'}),
    'heights': (['time', 'z'], heights),
    'z_r': (['z_r'], retrival_heights),
    'k': (['time','z','z_r'],ks),
})

# save the results
ch.to_netcdf('IR1Feb17vertest_abs_1d.nc')

# %% load
ch=xr.load_dataset('IR1Feb17vertest_abs_1d.nc')

# invert
weight_0 = 0
weight_1 = 0.7e-3
ver = xainvert(ch,retrival_heights, weight_0, weight_1)


result_1d = xr.Dataset().update({
        'time': (['time'], ch.time.values),
        'z_r': (['z_r',], ch.z_r.values, {'units': 'km'}),
        'ver': (['time','z_r'], ver, {'long name': 'VER', 'units': 'photons cm-3 s-1'}),
        'latitude': (['time',], ch.TPlat.values),
        'longitude': (['time',], ch.TPlon.values),
        'channel': (['time',], ch.channel.values),
        })

        
ir1band=result_1d

plt.figure(figsize=(12,3))
(ch.profile[30:230]/1e9*3.57/0.60).plot(y='z',vmin=0)
plt.figure(figsize=(12,2))
(ir1band.ver[30:230:1,5:-5]).plot.pcolormesh(y='z_r')
plt.title('A-band intensity (full band) photons cm-3 s-1')
#plt.savefig('Aband_intenstiy_fullband_feb17.png',format='png')

# %%
import matplotlib.pyplot as plt
plt.figure(figsize=(12,3))
ch.profile.where(ch.TPsza <98,drop=True).plot(y='z',vmin=0)
# %%
ch.attrs["title"]="IR1"
# %%
plt.figure(figsize=(12,3))
(ch.profile[30:230]/1e9*3.57/0.60).plot(y='z',vmin=0)
plt.figure(figsize=(12,3))
(ir1band.ver[30:230:1,:]).plot.pcolormesh(y='z_r')
plt.title('A-band intensity (full band) photons cm-3 s-1')
#plt.savefig('Aband_intenstiy_fullband_feb17.png',format='png')
# %%
plt.figure(figsize=(12,2))
m=plt.pcolormesh(ir1band.latitude[30:230],ir1band.z_r,ir1band.ver[30:230:1,:].T,vmin=1e11,vmax=5e11)
#m=plt.pcolormesh(np.log(ir1band.ver[30:230:1,:].T),vmin=25.5, vmax=27)

plt.colorbar(m)
#(result_1d.ver[30:50]/1e9*3.57/0.60).plot.line(y='z_r',add_legend=False);
plt.title('A-band intensity (full band) photons cm-3 s-1')

# %%
