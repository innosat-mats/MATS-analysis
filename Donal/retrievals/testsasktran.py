#%%
import numpy as np
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos, findheight,findtangent
from mats_l1_processing.pointing import pix_deg
import matplotlib.pylab as plt
from scipy.io import loadmat
from scipy.interpolate import interp1d
from skyfield import api as sfapi
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.units import Distance
from skyfield.positionlib import Geocentric, ICRF
from scipy.spatial.transform import Rotation as R
from numpy.linalg import norm

import sasktran as sk
from sasktran.geometry import VerticalImage
import mjdtools
%matplotlib widget

 # %%
def prepare_profile(ch):
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    profile = np.array(image[0:-10, col-2:col+2].mean(axis=1)*1e12)
    #profile = profile*1000/ch.TEXPMS #until files fixed
    common_heights = np.arange(60000,110250,250)
    profile=np.interp(common_heights,heights[0:-10],profile)
    return common_heights, profile, heights

# %%
data=loadmat('./sm_ABand_2022_2023_v4.mat',squeeze_me=True)
lat=np.array(data['sm']['lat90'])
lon=np.array(data['sm']['lon90'])
mjd=data['sm']['time90']
dates=[mjdtools.mjd2utc(m) for m in mjd]
dates = [d.replace(tzinfo=timezone.utc) for d in dates]
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

def makeSKprofiles(mjd,tanalts=np.arange(50,110),observer=[0,0,0],lookvector=[0,1,0]):
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
    engine.num_diffuse_profiles = 5

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
malb0,malb1=makeSKprofiles(mjd[i],tanalts=zt[i],observer=observer[i],lookvector=lookvec[i])
# %%
plt.figure()
plt.semilogx(malb0,zt[i])
plt.semilogx(malb1,zt[i])
plt.semilogx(L3[i]-L3[i][-2],zt[i],'.')
plt.semilogx(L4[i]-L4[i][-2],zt[i],'.')
#plt.semilogx(L3[i],zt[i],'.')
#plt.semilogx(L4[i],zt[i],'.')
plt.ylim(20,110)
plt.grid()
plt.legend(['s3-0','s4-0','s3-1','s4-1','L3','L4'])

plt.title("{} sza={:5.2f} ssa={:5.2f}".format(dates[i],sza[i],ssa[i]))  

plt.figure()
#plt.semilogx(malb0,zt[i])
#plt.semilogx(malb1,zt[i])
plt.plot((L3[i]-L3[i][-4])/malb0[:,0],zt[i],'.')
plt.plot((L4[i]-L4[i][-4])/malb0[:,1],zt[i],'.')
#plt.plot(L3[i]/malb0[:,0],zt[i])
#plt.plot(L4[i]/malb0[:,1],zt[i])
plt.legend(['L3/s3','L4/s4'])
plt.ylim(20,110)
plt.xlim([0,5])
plt.grid()

plt.title("{} sza={:5.2f} ssa={:5.2f}".format(dates[i],sza[i],ssa[i]))  


# %%
target_date = datetime(2023,2, 15, 13,26,21 , tzinfo=timezone.utc)
closest_index = np.argmin(np.abs(np.array(dates) - target_date))
print (closest_index)

# %%
plt.figure()
plt.hist(sza,bins=range(85,95,1))
count,classes= np.histogram(sza,bins=range(85,95,1))
# %%
#find the indices of the sza for each class given by classes
bin_indices = np.digitize(sza, classes)
indices_by_bin = {i: np.where(bin_indices == i)[0] for i in range(1, len(classes))}

# %%
m0iL3=[]
m0iL4=[]
m1iL3=[]
m1iL4=[]
L3i=[]
L4i=[]
ztcommon=np.arange(20,110,2)
#for i in range(1, len(classes)):
for i in [7]:
    for index in indices_by_bin[i]:
        malb0,malb1=makeSKprofiles(mjd[index],tanalts=zt[index],observer=observer[index],lookvector=lookvec[index])
        if zt[index][0] > zt[index][-1] : dir = -1
        else : dir = 1      
        f_L3 = interp1d(zt[index][::dir], L3[index][::dir], bounds_error=False, fill_value="extrapolate")
        f_L4 = interp1d(zt[index][::dir], L4[index][::dir], bounds_error=False, fill_value="extrapolate")
        malb0, malb1 = makeSKprofiles(mjd[index], tanalts=zt[index], observer=observer[index], lookvector=lookvec[index])
        f_malb0L3 = interp1d(zt[index][::dir], malb0[::dir, 0], bounds_error=False, fill_value="extrapolate")
        f_malb0L4 = interp1d(zt[index][::dir], malb0[::dir, 1], bounds_error=False, fill_value="extrapolate")

        f_malb1L3 = interp1d(zt[index][::dir], malb1[::dir, 0], bounds_error=False, fill_value="extrapolate")
        f_malb1L4 = interp1d(zt[index][::dir], malb1[::dir, 1], bounds_error=False, fill_value="extrapolate")

        L3i.append(f_L3(ztcommon))
        L4i.append(f_L4(ztcommon))
        m0iL3.append(f_malb0L3(ztcommon))
        m0iL4.append(f_malb0L4(ztcommon))
        m1iL3.append(f_malb1L3(ztcommon))
        m1iL4.append(f_malb1L4(ztcommon))
     



# %%
plt.figure()
for i,index in enumerate(indices_by_bin[1]):
    #plt.plot((L3[index]-L3[index][-4])/m0[i][:,0],zt[index],'.')
    #plt.plot((L4[index]-L4[index][-4])/m0[i][:,1],zt[index],'.')
    
    plt.plot((L3i[i])/m0iL3[i],ztcommon,'.')
    plt.plot((L4i[i])/m0iL4[i],ztcommon,'.')
    plt.legend(['L3/s3','L4/s4'])
    plt.ylim(20,110)
    plt.xlim([0,5])
    plt.grid()
# %%
plt.figure()


# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L3i) / np.asarray(m0iL3), axis=0)
std_values = np.nanstd(np.asarray(L3i) / np.asarray(m0iL3), axis=0)

# Plot the mean
plt.plot(mean_values, ztcommon, label='Mean Albedo = 0')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')

# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L3i) / np.asarray(m1iL3), axis=0)
std_values = np.nanstd(np.asarray(L3i) / np.asarray(m1iL3), axis=0)

# Plot the mean
plt.plot(mean_values, ztcommon, label='Mean Albedo = 1')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')

plt.legend()
plt.ylim(20,110)
plt.xlim([0,5])
plt.title ('Ratio between L3 and saskran s3')
plt.grid()
plt.figure()


# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L4i) / np.asarray(m0iL4), axis=0)
std_values = np.nanstd(np.asarray(L4i) / np.asarray(m0iL4), axis=0)

# Plot the mean
plt.plot(mean_values, ztcommon, label='Mean Albedo = 0')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')

# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L4i) / np.asarray(m1iL4), axis=0)
std_values = np.nanstd(np.asarray(L4i) / np.asarray(m1iL4), axis=0)

# Plot the mean
plt.plot(mean_values, ztcommon, label='Mean Albedo = 1')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')


plt.legend()
plt.title ('Ratio between L4 and saskran s4')
plt.ylim(20,110)
plt.xlim([0,5])
plt.grid()
# %%
plt.figure()


# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L3i) - np.asarray(m0iL3), axis=0)
std_values = np.nanstd(np.asarray(L3i) - np.asarray(m0iL3), axis=0)

# Plot the mean
plt.semilogx(mean_values, ztcommon, label='Mean Albedo = 0')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')

#Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L3i) - np.asarray(m1iL3), axis=0)
std_values = np.nanstd(np.asarray(L3i) - np.asarray(m1iL3), axis=0)

# Plot the mean
plt.semilogx(mean_values, ztcommon, label='Mean Albedo = 1')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')


plt.legend()
plt.title ('Difference between L3 and saskran s3')
plt.ylim(20,110)
#plt.xlim([0,5])
plt.grid()
plt.figure()

i=1
# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L4i) - np.asarray(m0iL4), axis=0)
std_values = np.nanstd(np.asarray(L4i) - np.asarray(m0iL4), axis=0)

# Plot the mean
plt.semilogx(mean_values, ztcommon, label='Mean Albedo = 0')

# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')

# Calculate the mean and standard deviation
mean_values = np.nanmean(np.asarray(L4i) - np.asarray(m1iL4), axis=0)
std_values = np.nanstd(np.asarray(L4i) - np.asarray(m1iL4), axis=0)

# Plot the mean
plt.semilogx(mean_values, ztcommon, label='Mean Albedo = 1')


# Add a shaded area for ± one standard deviation
plt.fill_betweenx(ztcommon, mean_values - std_values, mean_values + std_values, alpha=0.3, label='±1 Std Dev')
plt.legend()
plt.title ('Difference between L4 and saskran s4')

plt.ylim(20,110)
#plt.xlim([0,5])
plt.grid()
# %%
sel=np.where(np.asarray(sza)>100)[0]
tops=[np.mean(L3[i][-1:-4:-1]) for i in sel]
middles=[np.mean(L3[i][np.where(np.asarray(zt)[i]>90) and np.where(np.asarray(zt)[i]<95)]) for i in sel]
plt.figure()
plt.plot(np.asarray(dates)[sel],tops,'.')
plt.ylim([-2e13,2e13])
plt.figure()
plt.plot(tops,middles,'.')
plt.xlim(-1e13,1e13)
plt.ylim(-1e13,1e13)

# %%
sel=np.where(np.asarray(sza)>100)[0]
tops=[np.mean(L1[i][-1:-4:-1]) for i in sel]
middles=[np.mean(L1[i][np.where(np.asarray(zt)[i]>90) and np.where(np.asarray(zt)[i]<95)]) for i in sel]
plt.figure()
plt.plot(np.asarray(dates)[sel],tops,'.')
plt.ylim([-2e13,2e13])
plt.figure()
plt.plot(middles,tops,'.')
plt.ylabel('top')
plt.xlabel('middle')
plt.title('L1')
plt.xlim(-1e13,10e13)
plt.ylim(-1e13,5e13)
# %%
sel=np.where(np.asarray(sza)>100)[0]
tops1=[np.mean(L2[i][(np.asarray(zt)[i]>=100) & (np.asarray(zt)[i]<=102)])for i in sel]
middles1=[np.mean(L2[i][(np.asarray(zt)[i]>87) & (np.asarray(zt)[i]<92)]) for i in sel]
plt.figure()
plt.plot(np.asarray(dates)[sel],tops1,'.')
plt.ylim([-2e13,2e13])
plt.figure()
plt.plot(middles1,tops1,'.',alpha=0.5)
#nan filter middles and tops 
middles1 = np.array(middles1)
tops1 = np.array(tops1)
middles = middles1[~np.isnan(middles1) & ~np.isnan(tops1)]
tops = tops1[~np.isnan(tops1) & ~np.isnan(middles1)]

# add straight line fit to tops versus middles and plot
from scipy.stats import linregress  
slope, intercept, r_value, p_value, std_err = linregress(middles,tops) 
plt.plot(middles, slope*np.asarray(middles) + intercept, 'r')
plt.ylabel('top')
plt.xlabel('middle')
plt.title('L2')
plt.xlim(0e13,7e13)
plt.ylim(0e13,6e13)
# %%
sel=np.where(np.asarray(sza)>100)[0]
tops1=[np.mean(L1[i][(np.asarray(zt)[i]>=100) & (np.asarray(zt)[i]<=102)]) for i in sel]
middles1=[np.mean(L1[i][(np.asarray(zt)[i]>=87) & (np.asarray(zt)[i]<=92)]) for i in sel]
plt.figure()
plt.plot(np.asarray(dates)[sel],tops1,'.')
plt.ylim([-2e13,2e13])
plt.figure()
plt.plot(middles1,tops1,'.',alpha=0.5)
#nan filter middles and tops 
middles1 = np.array(middles1)
tops1 = np.array(tops1)
middles = middles1[~np.isnan(middles1) & ~np.isnan(tops1)]
tops = tops1[~np.isnan(tops1) & ~np.isnan(middles1)]

# add straight line fit to tops versus middles and plot
from scipy.stats import linregress  
slope, intercept, r_value, p_value, std_err = linregress(middles,tops) 
plt.plot(middles, slope*np.asarray(middles) + intercept, 'r')
plt.ylabel('top')
plt.xlabel('middle')
plt.title('L1')
plt.xlim(0e13,7e13)
plt.ylim(0e13,5e13)
# %%
OSsza = []
OSTPpos = []

for j in range(len(observer[i])):
    to_ecef = R.from_matrix(itrs.rotation_at(times[i]))
    sunpos = sun.at(times[i]).position.km

    lveci = to_ecef.inv().apply(lookvec[i][j, :])
    ecipos = to_ecef.inv().apply(observer[i][j, :])
    res = findtangent(times[i], ecipos, lveci)
    TPposeci = ecipos + res.x * lveci
    TPpos = wgs84.geographic_position_of(
        Geocentric(position_au=Distance(m=TPposeci).au, t=times[i])
    )

    OSTPpos.append([
        TPpos.latitude.degrees,
        TPpos.longitude.degrees,
        TPpos.elevation.m
    ])
    OSsza.append(
        90 - ((earth + TPpos).at(times[i]).observe(sun).apparent().altaz()[0].degrees)
    )
# %%

