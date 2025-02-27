# %%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
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

# %%
#dftop = pd.read_pickle('/Users/donal/projekt/SIW/verdec')
splined2dlogfactor=np.load("splined2dlogfactorsIR1.npy",allow_pickle=True).item()
tanz, splinedfactor=np.load("splinedlogfactorsIR1_feb.npy",allow_pickle=True)



#%%
starttime=datetime(2023,2,17,0,50)
stoptime=datetime(2023,2,17,1,20)
#starttime=datetime(2023,3,29,21,0)
#stoptime=datetime(2023,3,29,21,35)
#starttime=datetime(2023,3,18,20,30)
#stoptime=datetime(2023,3,18,21,20)


dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.6")
#%%
# df=df[df['channel']!='NADIR']
df = dftop[dftop['channel'] == 'IR1'].dropna().reset_index()#[0:10]
# %%
#select part of orbit
offset = 10
num_profiles = 250 #use 50 profiles for inversion
df = df.loc[offset:offset+num_profiles]
df = df.reset_index(drop=True)
# %%

def prepare_profile(ch,hotpiximage):
    image = image = np.stack(ch.ImageCalibrated)#-hotpiximage
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[:, col-2:col+2].mean(axis=1)*1e15)
    return heights, profile


# %%
profiles = []
heights = []
ecipos = []
ecivecs = []
ks=[]
retrival_heights= np.arange(60,111,1)
#retrival_heights= np.arange(70,100,1)
s_140=1700e3
steps=100 #m steps
s_steps = np.arange(s_140,s_140 + 2e6,steps)
ts = sfapi.load.timescale()
#plt.figure()
# No need for hotpiximage
#if df.channel[0] == 'IR1': hotpiximage=np.load ('ir1mean.npy')
#elif df.channel[0] == 'IR2': hotpiximage=np.load ('/Users/donal/projekt/SIW/ir2mean.npy')
#elif df.channel[0] == 'IR3': hotpiximage=np.load ('ir3mean.npy')
#elif df.channel[0] == 'IR4': hotpiximage=np.load ('ir4mean.npy')
hotpiximage=None

#plt.figure()

# %%
for i in range(len(df)):
    # for i in range(100):
    print(i)
    zs, p = prepare_profile(df.iloc[i],hotpiximage)
    if df.iloc[i]['channel'] == 'IR1':
        zs = zs[5:-5]
        p = p[5:-5]
        #nrows = len(zs) ### EDIT df['NROW'][i]
    profiles.append(p)
    heights.append(zs)
    ecipos.append(df['afsGnssStateJ2000'][i][0:3])
    d = df['EXPDate'][i]
    t = ts.from_datetime(d)
    localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat[i], df.TPlon[i], elevation_m=0).at(t).position.m)
    q = df['afsAttitudeState'][i]
    

    quat = R.from_quat(np.roll(q, -1))
    ypixels = np.linspace(0, df['NROW'][i], 5)
    x, yv = pix_deg(df.iloc[i], int(df['NCOL'][i]/2), ypixels)
    qp = R.from_quat(df['qprime'][i])
    ecivec = np.zeros((3, len(yv)))
    k=np.zeros((df['NROW'][i],len(retrival_heights)))
    for irow, y in enumerate(yv):
        los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
        ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
    cs_eci=CubicSpline(ypixels,ecivec.T)
    for irow in range(df['NROW'][i]):
        ecivec=cs_eci(irow)
        pos=np.expand_dims(ecipos[-1], axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
        point_height=np.linalg.norm(pos,axis=0)-localR
        
        ### add weight here?
        # splinedlogfactor
        minarg=point_height[:].argmin() # pos of lowest tan height
        target_tangent=zs[irow]/1000
        distances=(s_steps-s_steps[minarg])/1000 #to km for intepolation
        lowertan=bisect_left(tanz,target_tangent)-1
        uppertan=lowertan+1
        lowerfactor=splinedfactor[lowertan](distances)
        upperfactor=splinedfactor[uppertan](distances)
        newfactor=interp1d([tanz[lowertan],tanz[uppertan]],np.array([lowerfactor,upperfactor]).T)
        weight=np.exp(newfactor(target_tangent))

        if irow < 20:
            plt.plot(weight)

        counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)),weights=weight)
        
            

        #counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)))
        k[irow,:]=counts*steps
        ecivecs.append(ecivec)
    ks.append(k)
plt.title('dec')
#plt.savefig('dec_abs.png',format='png')
ecivecs= np.reshape(ecivecs,(len(df),-1,3))
#%%    
z= np.array(heights).mean(axis=0)
inputdata = xr.Dataset({
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
# %%
inputdata.to_netcdf('IR1Feb17vertest_abs_1d.nc')
# %%
