# %%
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline
from skyfield import api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
import xarray as xr
from numpy.linalg import inv

# %%
dftop = pd.read_pickle('/Users/donal/projekt/SIW/verdec')
#%%
#starttime=datetime(2023,3,29,21,0)
#stoptime=datetime(2023,3,29,22,35)
#dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.4")
#%%
# df=df[df['channel']!='NADIR']
df = dftop[dftop['channel'] == 'IR2'].dropna().reset_index()#[0:10]
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
s_140=1700e3
steps=100 #m steps
s_steps = np.arange(s_140,s_140 + 2e6,steps)
ts = sfapi.load.timescale()
if df.channel[0] == 'IR1': hotpiximage=np.load ('ir1mean.npy')
elif df.channel[0] == 'IR2': hotpiximage=np.load ('/Users/donal/projekt/SIW/ir2mean.npy')
elif df.channel[0] == 'IR3': hotpiximage=np.load ('ir3mean.npy')
elif df.channel[0] == 'IR4': hotpiximage=np.load ('ir4mean.npy')
for i in range(len(df)):
    # for i in range(100):
    print(i)
    zs, p = prepare_profile(df.iloc[i],hotpiximage)
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
        counts,bins=np.histogram(point_height/1000,np.hstack((retrival_heights,retrival_heights[-1]+1)))
        k[irow,:]=counts*steps
        ecivecs.append(ecivec)
    ks.append(k)
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
inputdata.to_netcdf('IR2Mar29vertest_1d.nc')
# %%
