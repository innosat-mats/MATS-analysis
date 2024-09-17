#%%
import numpy as np
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos, findheight
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.io import loadmat
from skyfield import api as sfapi
from skyfield.framelib import itrs
from scipy.spatial.transform import Rotation as R
from numpy.linalg import norm

import sasktran as sk
from sasktran.geometry import VerticalImage
import mjdtools
# %%
today=datetime.now()
mjd=mjdtools.utc2mjd(today.year,today.month,today.day,today.hour,today.minute,today.second)
geometry = VerticalImage()
tanalts=np.arange(50,110)
geometry.from_sza_ssa(sza=80, ssa=90, lat= -85, lon= 10, tanalts_km=tanalts, mjd=mjd, locallook=0,
                        satalt_km=600, refalt_km=20)

atmosphere = sk.Atmosphere()
atmosphere.atmospheric_state=sk.MSIS90(max_altitude_km=140)
atmosphere.brdf=1

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
modelalb1=  irr*radiance.T 
atmosphere.brdf=0
radiance = engine.calculate_radiance()
spectrum=sk.SolarSpectrum()
irr=spectrum.irradiance(engine.wavelengths)*1e4 #to m**2
modelalb0=  irr*radiance.T 
plt.figure()
plt.semilogx(modelalb0,tanalts)
plt.title('Albedo = 0')
plt.legend(['s3','s4'])
plt.figure()
plt.semilogx(modelalb1,tanalts)
plt.title('Albedo = 1')
plt.legend(['s3','s4'])
plt.figure()
plt.plot(modelalb1/modelalb0,tanalts)
plt.legend(['s3','s4'])
plt.title('Ratio Albedo = 1 to Albedo = 0')
#fig, ax = plt.subplots()
# %%
data=loadmat('./sm_ABand_2022_2023_v4.mat',squeeze_me=True)
lat=np.array(data['sm']['lat90'])
lon=np.array(data['sm']['lon90'])
mjd=data['sm']['time90']
dates=[mjdtools.mjd2utc(m) for m in mjd]
ts = sfapi.load.timescale()
times=[ts.from_datetime(d) for d in dates]
# %%
for ch in ['L1','L2','L3','L4'] :
    locals()[ch]=data['sm'][ch]
zt=data['sm']['tangent_altitude']    
lookvec=data['sm']['look_vector']
observer=data['sm']['observer_position']
to_ecef = R.from_matrix(itrs.rotation_at(times[0]))

# %%
eph=sfapi.load('de421.bsp')

sun=eph['sun']
earth=eph['earth']
sza=[]
ssa=[]
for i in range(len(dates)):
    sza.append(90-((earth+sfapi.wgs84.latlon(lat[i],lon[i],elevation_m=90000)).at(times[i]).observe(sun).apparent().altaz()[0].degrees) )
    to_ecef = R.from_matrix(itrs.rotation_at(times[i]))
    sunpos=sun.at(times[i]).position.km
    closest_index = np.abs(zt[i] - 90).argmin()
    lveci=to_ecef.inv().apply(lookvec[i][closest_index,:])    
    ssa.append(np.arccos(np.dot(sunpos,lveci)/(norm(sunpos)*norm(lveci)))*180/np.pi)
# %%

def makeSKprofiles(mjd,lat=-85.0,lon=10.0,ssa=90.0,sza=90.0,tanalts=np.arange(50,110),observer=[0,0,0],lookvector=[0,1,0]):
    geometry = VerticalImage()
    los=[sk.LineOfSight(mjd=mjd,observer=observer[iz],look_vector=lookvector[iz]) for iz in range(len(tanalts))]
    #geometry.from_sza_ssa(sza=sza, ssa=ssa, lat= lat, lon= lon , tanalts_km=tanalts, mjd=mjd, locallook=0,
    #                        satalt_km=600, refalt_km=90)
    geometry.lines_of_sight=los
    atmosphere = sk.Atmosphere()
    atmosphere.atmospheric_state=sk.MSIS90(max_altitude_km=140)
    atmosphere.brdf=1

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
    modelalb1=  irr*radiance.T 
    atmosphere.brdf=0
    radiance = engine.calculate_radiance()
    spectrum=sk.SolarSpectrum()
    irr=spectrum.irradiance(engine.wavelengths)*1e4 #to m**2
    modelalb0=  irr*radiance.T 
    return(modelalb0,modelalb1)
# %%
i=2866
print(dates[i],ssa[i],sza[i])    
malb0,malb1=makeSKprofiles(mjd[i],lat=lat[i],lon=lon[i],ssa=ssa[i],sza=sza[i],tanalts=zt[i],observer=observer[i],lookvector=lookvec[i])
# %%
plt.figure()
#plt.semilogx(malb0,zt[i])
plt.semilogx(malb1,zt[i])
plt.semilogx(L3[i]-L3[i][-4],zt[i],'.')
plt.semilogx(L4[i]-L4[i][-4],zt[i],'.')
plt.ylim(40,110)
plt.grid()
plt.legend(['s3','s4','L3','L4'])


# %%
target_date = datetime(2023, 2, 15, 13,19 )
closest_index = np.argmin(np.abs(np.array(dates) - target_date))
