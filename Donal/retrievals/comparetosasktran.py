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

import sasktran as sk
from sasktran.geometry import VerticalImage
import mjdtools

# %%
#tanz, splinedfactor=np.load(expanduser("~donal/projekt/SIW/splinedlogfactorsIR1.npy"),allow_pickle=True)
#splined2dlogfactor=np.load(expanduser("~donal/projekt/SIW/splined2dlogfactorsIR1.npy"),allow_pickle=True).item()
dftop = pd.read_pickle(expanduser('~donal/projekt/SIW/verdec'))
#%%
starttime=datetime(2023,3,23,21,0)
stoptime=datetime(2023,3,23,22,35)
#starttime=datetime(2023,5,5,21,0)
#stoptime=datetime(2023,5,5,22,35)
starttime=datetime(2022,12,25,21,0)
stoptime=datetime(2022,12,25,22,35)
dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.5")
#%%
# with open('verdec2d.pickle', 'wb') as handle:
#     pickle.dump(dftop, handle)#%%
# df=df[df['channel']!='NADIR']
ir1 = dftop[dftop['channel'] == 'IR1'].dropna().reset_index(drop=True)#[0:10]
ir2 = dftop[dftop['channel'] == 'IR2'].dropna().reset_index(drop=True)#[0:10]
ir3 = dftop[dftop['channel'] == 'IR3'].dropna().reset_index(drop=True)#[0:10]
ir4 = dftop[dftop['channel'] == 'IR4'].dropna().reset_index(drop=True)#[0:10]
uv1 = dftop[dftop['channel'] == 'UV1'].dropna().reset_index(drop=True)#[0:10]
uv2 = dftop[dftop['channel'] == 'UV2'].dropna().reset_index(drop=True)#[0:10]

# %%

# %%


def prepare_profile(ch):
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[0:-10, col-2:col+2].mean(axis=1)*1e13)
    profile = profile*1000/ch.TEXPMS #until files fixed
    common_heights = np.arange(60000,110250,250)
    profile=np.interp(common_heights,heights[0:-10],profile)
    return common_heights, profile, heights


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
offset = 0
num_profiles = 900 #use 50 profiles for inversion
ir1 = ir1.loc[offset:offset+num_profiles].reset_index(drop=True)
ir2 = ir2.loc[offset:offset+num_profiles].reset_index(drop=True)
ir3 = ir3.loc[offset:offset+num_profiles].reset_index(drop=True)
ir4 = ir4.loc[offset:offset+num_profiles].reset_index(drop=True)
#offset=offset-414
#uv1 = uv1.loc[offset:offset+num_profiles].reset_index(drop=True)
#uv2 = uv2.loc[offset:offset+num_profiles].reset_index(drop=True)
#%%
tanalts_km = np.arange(60, 110.25, 0.25)
szas=[]
ssas=[]
azis=[]

scales=[]
means=[]
intens=[]
# First recreate our geometry and atmosphere classes
geometry = VerticalImage()
for i in range(0,num_profiles,10):
    z1,p1,orig_heights=prepare_profile(ir1.iloc[i])
    z2,p2,_=prepare_profile(ir2.iloc[i])
    z3,p3,_=prepare_profile(ir3.iloc[i])
    z4,p4,_=prepare_profile(ir4.iloc[i])
    #z1,puv1=prepare_profile(uv1.iloc[i])
    #z2,puv2=prepare_profile(uv2.iloc[i])
    plt.figure()
    #plt.semilogx(p1,z1,p2,z2,p3-p3[-1]/1.1,z3,p4-p4[-1]/1.1,z4) #,puv1,z1,puv2,z2)
    plt.semilogx(p1,z1,p2,z2,p3,z3,p4,z4) #,puv1,z1,puv2,z2)
    plt.legend(['IR1','IR2','IR3','IR4','UV1','UV2',])
    plt.title('TP sza {:5.2f} TP ssa {:5.2f}  TP lat {:5.2f} TP lon {:5.2f}'.format(ir3.TPsza[i],ir3.TPssa[i],ir3.TPlat[i],ir3.TPlon[i]))
    #plt.xlim([1e15, 3e17])
    # plt.figure()
    # plt.plot(p1-(p3+p4)/2+(p3[-1]+p4[-1])*1.1/2,z1,p2-(p3+p4)/2+(p3[-1]+p4[-1])*1.1/2,z2)
    # plt.legend(['IR1','IR2','IR3','IR4',])

    mjd=mjdtools.utc2mjd(ir3.EXPDate[i].year,ir3.EXPDate[i].month,ir3.EXPDate[i].day,ir3.EXPDate[i].hour,ir3.EXPDate[i].minute,ir3.EXPDate[i].second)
    geometry.from_sza_ssa(sza=ir3.TPsza[i], ssa=ir3.TPssa[i], lat=ir3.TPlat[i], lon=ir3.TPlon[i], tanalts_km=tanalts_km, mjd=mjd, locallook=0,
                        satalt_km=600, refalt_km=20)
    szas.append(ir3.TPsza[i])
    ssas.append(ir3.TPssa[i])
    azis.append(ir3.nadir_az[i])
    atmosphere = sk.Atmosphere()
    atmosphere.atmospheric_state=sk.MSIS90(max_altitude_km=140)


    atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
    atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90(max_altitude_km=140))

    # And now make the engine
    engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
    engine.top_of_atmosphere_altitude = 111000
    # Choose some wavelengths to do the calculation at
    engine.wavelengths = [754, 772]

    # And do the calculation
    radiance = engine.calculate_radiance()
    spectrum=sk.SolarSpectrum()
    irr=spectrum.irradiance(engine.wavelengths)*1e4 #to m**2
    #fig, ax = plt.subplots()
    plt.semilogx (irr*radiance.T, tanalts_km*1000)
    plt.legend(['IR1','IR2','IR3','IR4','s3','s4'])
    plt.show()
    plt.figure()
    meas=np.array([p3,p4]).T
    model=  irr*radiance.T  
    scale=meas[10,:]/model[10,:]
    scales.append(scale)
    plt.plot(meas-scale*model,tanalts_km*1000)
    difference=meas-scale*model
    means.append(difference[40:160,:].mean(axis=0))
    intens.append((np.array([p1,p2]).T)[40:160,:].sum(axis=1))
    plt.legend(['s3','s4'])
    plt.title('scaling factors = {:5.2f},{:5.2f}'.format(*scale))
    plt.show()
# %%
plt.figure()
plt.plot(szas,means)
plt.title('Mean Straylight')
plt.xlabel('SZA')
plt.legend(['s3','s4'])
plt.figure()
plt.plot(szas,scales)
plt.title('Scale')
plt.xlabel('SZA')
plt.legend(['s3','s4'])

plt.figure()
plt.plot(ssas,means)
plt.title('Mean Straylight')
plt.xlabel('SSA')
plt.legend(['s3','s4'])
plt.figure()
plt.plot(ssas,scales)
plt.title('Scale')
plt.xlabel('SSA')
plt.legend(['s3','s4'])

plt.figure()
plt.plot(azis,means)
plt.title('Mean Straylight')
plt.xlabel('Azi')
plt.legend(['s3','s4'])
plt.figure()
plt.plot(azis,scales)
plt.title('Scale')
plt.xlabel('Azi')
plt.legend(['s3','s4'])

    # %%
plt.figure()
plt.plot(szas,scales)
plt.title('Scale')
plt.xlabel('SZA')
plt.legend(['s3','s4'])
plt.ylim([0,2])
plt.figure()
plt.plot(ssas,scales)
plt.title('Scale')
plt.xlabel('SSA')
plt.legend(['s3','s4'])
plt.ylim([0,2])
plt.figure()
plt.plot(azis,scales)
plt.title('Scale')
plt.xlabel('Azi')
plt.legend(['s3','s4'])
plt.ylim([0,2])
# %%
