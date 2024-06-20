#%%
import numpy as np
import pandas as pd
from numba import njit
from bisect import bisect_left
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos, findheight
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import CubicSpline,interp1d
import scipy.sparse as sp
from skyfield import api as sfapi
from skyfield.framelib import itrs
from os.path import expanduser
import xarray as xr
import pickle
import plotly.graph_objects as go
from fast_histogram import histogramdd

# %%
tanz, splinedfactor=np.load(expanduser("~donal/projekt/SIW/splinedlogfactorsIR1.npy"),allow_pickle=True)
splined2dlogfactor=np.load(expanduser("~donal/projekt/SIW/splined2dlogfactorsIR1.npy"),allow_pickle=True).item()
dftop = pd.read_pickle(expanduser('~donal/projekt/SIW/verdec'))
#%%
# starttime=datetime(2023,3,31,21,0)
# stoptime=datetime(2023,3,31,22,35)
# dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.4")

# with open('verdec2d.pickle', 'wb') as handle:
#     pickle.dump(dftop, handle)#%%
# df=df[df['channel']!='NADIR']
df = dftop[dftop['channel'] == 'IR1'].dropna().reset_index(drop=True)#[0:10]
# %%
# %%


def prepare_profile(ch):
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW']-10)))
    profile = np.array(image[0:-10, col-2:col+2].mean(axis=1)*1e15)
    return heights, profile


@njit(cache=True)
def cart2sph(pos=np.array([[]])):
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    #radius =np.linalg.norm(pos,axis=1)
    longitude = np.arctan2(y, x)
    latitude = np.arcsin(z / radius)
    return radius, longitude,latitude
    #return np.array([radius,longitude,latitude]).T

# %%
profiles = []
heights = []
ecipos = []
ecivecs = []
posecef_sph=[]
retrival_heights= np.arange(30,130,1)

ts = sfapi.load.timescale()

#select part of orbit
offset = 100
num_profiles = 100 #use 50 profiles for inversion
df = df.loc[offset:offset+num_profiles]
df = df.reset_index(drop=True)
k_row = 0

#Generate grid for mid measurement
first = 0
mid = int(num_profiles/2)
last = num_profiles

to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][first])))
posecef_first = to_ecef.apply(df.afsTangentPointECI[first]).astype('float32')
to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][mid])))
posecef_mid = to_ecef.apply(df.afsTangentPointECI[mid]).astype('float32')
to_ecef = R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][last])))
posecef_last = to_ecef.apply(df.afsTangentPointECI[last]).astype('float32')
observation_normal = np.cross(posecef_first,posecef_last)
observation_normal = observation_normal/np.linalg.norm(observation_normal)
posecef_mid_unit = posecef_mid/np.linalg.norm(posecef_mid)
ecef_to_local = R.align_vectors([[1,0,0],[0,1,0]],[posecef_mid_unit,observation_normal])[0]

ecipos.append(df['afsGnssStateJ2000'][mid][0:3])
d = df['EXPDate'][mid]
t = ts.from_datetime(d)
localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat[mid], df.TPlon[mid], elevation_m=0).at(t).position.m)
q = df['afsAttitudeState'][mid]
quat = R.from_quat(np.roll(q, -1))

#calulate tangent geometries for liminted rows and make a spline for it
tanheights, profile = prepare_profile(df.iloc[mid])
ypixels = np.linspace(0, df['NROW'][mid],5)
x, yv = pix_deg(df.loc[mid], int(df['NCOL'][mid]/2), ypixels)
qp = R.from_quat(df['qprime'][mid])
ecivec = np.zeros((3, len(yv)))
for irow, y in enumerate(yv):
    los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
    ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
cs_eci=CubicSpline(ypixels,ecivec.T)
to_ecef=R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][mid])))

#select row to calculate jacobian for
irow = 0
ecivec=cs_eci(irow)
#%%
zs=np.linalg.norm(ecipos[-1])
theta=np.arccos(np.dot(ecipos[-1],ecivec)/zs)
b=2*zs*np.cos(theta)
root=np.sqrt(b**2+4*((120e3+localR)**2 - zs**2))
s_120_1 =(-b-root)/2
s_120_2 =(-b+root)/2

#%%
#s_120_1=findheight(t,ecipos[-1],ecivec,140e3,bracket=(10e5,20e5)).x
#s_120_2=findheight(t,ecipos[-1],ecivec,140e3,bracket=(30e5,60e5)).x
steps=10 #m steps
s_steps = np.arange(s_120_1,s_120_2,steps)
pos=np.expand_dims(ecipos[-1], axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
posecef_i=(to_ecef.apply(pos.T).astype('float32'))
posecef_i = ecef_to_local.apply(posecef_i) #convert to local
posecef_i_sph = cart2sph(posecef_i)   #x: height, y: acrosstrac (angle), z:  along track (angle)
altitude_grid = np.arange(localR+10e3,localR+120e3,1e3)

posecef_i_sph=np.array(posecef_i_sph).T
acrosstrack_grid = np.array([-0.3,0.3])
#acrosstrack_grid = np.linspace(posecef_i_sph[:,1].min(),posecef_i_sph[:,1].max(),1)
alongtrack_grid = np.linspace(posecef_i_sph[:,2].min(),posecef_i_sph[:,2].max(),100)
alongtrack_grid[0] = alongtrack_grid[0]-0.5
alongtrack_grid[-1] = alongtrack_grid[-1]+0.5

#Calculate the weights 
minarg=posecef_i_sph[:,0].argmin()
distances=(s_steps-s_steps[minarg])/1000 #to km for intepolation
target_tangent=tanheights[0]/1000
lowertan=bisect_left(tanz,target_tangent)-1
uppertan=lowertan+1
lowerfactor=splinedfactor[lowertan](distances)
upperfactor=splinedfactor[uppertan](distances)
newfactor=interp1d([tanz[lowertan],tanz[uppertan]],np.array([lowerfactor,upperfactor]).T)
weight=np.exp(newfactor(target_tangent))
#weight=np.exp(splined2dlogfactor(target_tangent,distances)).squeeze()
#hist, edges = np.histogramdd(posecef_i_sph[::1,:],bins=[altitude_grid,acrosstrack_grid,alongtrack_grid]) 
hist = histogramdd(posecef_i_sph[::1,:],bins=[110,1,120],range=[[altitude_grid[0],altitude_grid[-1]],[acrosstrack_grid[0],acrosstrack_grid[-1]],[alongtrack_grid[0],alongtrack_grid[-1]]],weights=weight)
edges=[]
edges.append(np.linspace(altitude_grid[0],altitude_grid[-1],110+1))
edges.append(np.linspace(acrosstrack_grid[0],acrosstrack_grid[-1],1+1))
edges.append(np.linspace(alongtrack_grid[0],alongtrack_grid[-1],120+1))

#%%
fig = go.Figure(data=[go.Scatter3d(x=posecef_i_sph[::1000,0]-localR, y=posecef_i_sph[::1000,1], z=posecef_i_sph[::1000,2],mode='markers')])
fig.show()

#%%
X, Y, Z = np.meshgrid(edges[0][:-1], edges[1][:-1], edges[2][:-1])#,indexing='ij')

fig = go.Figure(data=go.Volume(
    x = X.reshape(-1),
    y = Y.reshape(-1),
    z = Z.reshape(-1),
    value=hist.reshape(-1),
    opacity=0.2, # needs to be small to see through all surfaces
    surface_count=10, # needs to be a large number for good volume rendering
    ))
fig.show()

#%%
k = hist.reshape(-1)
ks = sp.lil_array((df['NROW'].sum()-10*len(df),len(k)))

#calculate jacobian for all measurements

profiles = []
heights = []
ecipos = []
ecivecs = []
posecef_sph=[]

for i in range(len(df)):
    print(i)

    zs, p = prepare_profile(df.iloc[i])
    profiles.append(p)
    heights.append(zs)

    ecipos.append(df.iloc[i]['afsGnssStateJ2000'][0:3])
    d = df.iloc[i]['EXPDate']
    t = ts.from_datetime(d)
    localR = np.linalg.norm(sfapi.wgs84.latlon(df.TPlat.iloc[i], df.TPlon.iloc[i], elevation_m=0).at(t).position.m)
    q = df['afsAttitudeState'].iloc[i]
    quat = R.from_quat(np.roll(q, -1))

    #calulate tangent geometries for center column limited rows and make a spline for it
    ypixels = np.linspace(0, df['NROW'].iloc[i], 5)
    x, yv = pix_deg(df.iloc[i], int(df['NCOL'].iloc[i]/2), ypixels)
    qp = R.from_quat(df['qprime'].iloc[i])
    ecivec = np.zeros((3, len(yv)))
    k=np.zeros((df['NROW'].iloc[i]-10,len(retrival_heights)))
    for irow, y in enumerate(yv):
        los = R.from_euler('xyz', [0, y, x], degrees=True).apply([1, 0, 0])
        ecivec[:, irow] = np.array(quat.apply(qp.apply(los)))
    cs_eci=CubicSpline(ypixels,ecivec.T)
    
    to_ecef=R.from_matrix(itrs.rotation_at(ts.from_datetime(df['EXPDate'][i])))
    satecefpos=to_ecef.apply(ecipos[-1])
    
    #calculate jacobian for each row
    for irow in range(df['NROW'].iloc[i]-10):
        ecivec=cs_eci(irow)
        ecefvec=to_ecef.apply(ecivec)
        zs=np.linalg.norm(ecipos[-1])
        theta=np.arccos(np.dot(ecipos[-1],ecivec)/zs)
        b=2*zs*np.cos(theta)
        root=np.sqrt(b**2+4*((120e3+localR)**2 - zs**2))
        s_120_1 =(-b-root)/2
        s_120_2 =(-b+root)/2
        #s_120_1=findheight(t,ecipos[-1],ecivec,120e3,bracket=(1e3,20e5)).x
        #s_120_2=findheight(t,ecipos[-1],ecivec,120e3,bracket=(30e5,40e5)).x
        #if ((s_120_2 - s_120_1) < 1e3) : print('wow',s_120_1,s_120_2)
        s_steps = np.arange(s_120_1,s_120_2,steps)
        #pos=np.expand_dims(ecipos[-1], axis=0).T+s_steps*np.expand_dims(ecivec, axis=0).T
        #posecef_i=(to_ecef.apply(pos.T).astype('float32'))
        posecef_i=(np.expand_dims(satecefpos, axis=0).T+s_steps*np.expand_dims(ecefvec, axis=0).T).astype('float32')
        posecef_i = ecef_to_local.apply(posecef_i.T) #convert to local (for middle alongtrack measurement)
        posecef_i_sph = cart2sph(posecef_i)   
        posecef_i_sph=np.array(posecef_i_sph) .T   
        #hist, _ = np.histogramdd(posecef_i_sph[::1,:],edges)
        #Calculate the weights 
        minarg=posecef_i_sph[:,0].argmin()
        distances=(s_steps-s_steps[minarg])/1000 #to km for intepolation
        target_tangent=tanheights[irow]/1000
        lowertan=bisect_left(tanz,target_tangent)-1
        uppertan=lowertan+1
        lowerfactor=splinedfactor[lowertan](distances)
        upperfactor=splinedfactor[uppertan](distances)
        newfactor=interp1d([tanz[lowertan],tanz[uppertan]],np.array([lowerfactor,upperfactor]).T)
        weight=np.exp(newfactor(target_tangent))
        #weight=np.exp(splined2dlogfactor(target_tangent,distances)).squeeze()
        hist = histogramdd(posecef_i_sph[::1,:],range=[[edges[0][0],edges[0][-1]],[edges[1][0],edges[1][-1]],[edges[2][0],edges[2][-1]]],bins=[110,1,120],weights=weight)
        k = hist.reshape(-1)
        ks[k_row,:] = k
        k_row = k_row+1
    
        
    
# ecivecs= np.reshape(ecivecs,(len(df),-1,3))
#posecef_sph= np.reshape(posecef_sph,(-1,3))

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
    'ecipos': (['time', 'xyz'], ecipos),   
    'z':  (['z'], z, {'long_name': 'Approx Altitude', 'units': 'm'}),
    'profile': (['time', 'z'], profiles, {'long_name': 'LOS  intensity', 'units': 'Photons m-2 nm-1 sr-1 s-1'}),
    'heights': (['time', 'z'],  heights),
    'ret_grid_z': (['z_r'], edges[0]),
    'ret_grid_lon': (['lon_r'], edges[1]),
    'ret_grid_lat': (['lat_r'], edges[2]),
})

inputdata.to_netcdf('IR1mars31vertest_100-200.nc')
# %%
filename = "jacobian1_100-200.pkl"
with open(filename, "wb") as file:
    pickle.dump((edges, ks, ecef_to_local), file)

# %%
