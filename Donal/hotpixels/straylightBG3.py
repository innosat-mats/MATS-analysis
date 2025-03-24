#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as DT
from scipy.signal import savgol_filter
import hpfunctions
import xarray as xr
import skyfield.api as sfapi
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance
from pandas import DataFrame, Timestamp, to_datetime 
import sasktran as sk
from sasktran.geometry import VerticalImage
%matplotlib widget

#%%
chan='IR3'
ch = xr.open_zarr(f"s3://test-release-v1.0.1/mats-level-1b-limb-cropd-{chan}.zarr", storage_options={"anon":True})
#%%
def makeSKprofiles(mjd,imageshape=[63,9],observer=[0,0,0],lookvector=[0,1,0]):
    geometry = VerticalImage()
    XX,YY=np.meshgrid(np.arange(imageshape[0]),np.arange(imageshape[1]),sparse=True)
    los=[sk.LineOfSight(mjd=mjd,observer=observer,look_vector=lookvector[ix,iy,:]) for ix in range(imageshape[0]) for iy in range(imageshape[1])]   
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
    #engine.num_diffuse_profiles = 5

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

def do_geoloc(ccditem):
    geocoords=[ccditem.geo_coord_x.to_numpy(), ccditem.geo_coord_y.to_numpy()]
    geodata=ccditem.geoloc.to_numpy()
    from scipy.interpolate import RegularGridInterpolator

    interpolator = RegularGridInterpolator(geocoords, geodata, method='cubic',bounds_error=False,fill_value=None)   
    fullxgrid = np.arange(ccditem["NCOL"] + 1)
    fullygrid = np.arange(ccditem["NROW"])
    XX, YY = np.meshgrid(fullxgrid, fullygrid, sparse=True)
    geoloc=interpolator((XX,YY))
    return geoloc
# %%
plt.figure()
plt.plot(ch.TPsza)
days=ch.sel(time=slice("2023-04-23T00:00","2023-04-29T23:59"),drop=True)
#%%
bright =days.where(days.TPsza < 90,drop=True)
dark = days.where(days.TPsza > 97,drop=True)
dark=dark.isel(time=range(0,len(dark.time),10))
bright=bright.isel(time=range(0,len(bright.time),10))
# %%
dark.ImageCalibrated[0].shape
mapdate,hp=hpfunctions.gethpm(dark.time[0],chan)
# %%
plt.figure(figsize=(6, 4))
dark.ImageCalibrated.mean("time").plot(robust=True)
plt.figure(figsize=(6, 4))
bright.ImageCalibrated.mean("time").plot(robust=True)
plt.figure(figsize=(6, 4))
plt.imshow(hp,origin='lower', aspect='auto')
plt.colorbar()
# %%
# %%
# %%
results0 = []
results1 = []

for imageno in range(len(bright.time)):
    d = to_datetime(bright.sel(time=bright.time[imageno])["time"].values).tz_localize("UTC").to_pydatetime()
    ts = sfapi.load.timescale()
    t = ts.from_datetime(d)
    ecipos = bright.sel(time=bright.time[imageno]).afsGnssStateJ2000[:3]
    ecefpos = ICRF(Distance(m=ecipos.T).au, t=t, center=399).frame_xyz(itrs).m.T
    geoloc = do_geoloc(bright.isel(time=imageno))
    ecefpos_repeated = np.tile(ecefpos, [63, 9, 1])
    los = geoloc[:, :, 1:] - ecefpos_repeated
    skprofs0, skprofs1 = makeSKprofiles(t.tt - 2400000.5, observer=ecefpos, lookvector=los)
    skimage1 = skprofs1[:, 0].reshape(63, -1)
    skimage0 = skprofs0[:, 0].reshape(63, -1)
    diff_image1 = bright.ImageCalibrated.isel(time=imageno) - skimage1  
    diff_image0 = bright.ImageCalibrated.isel(time=imageno) - skimage0
    results0.append(diff_image0)
    results1.append(diff_image1)

bright['ImageCalibratedMinusSKImage0'] = xr.concat(results0, dim='time')
bright['ImageCalibratedMinusSKImage1'] = xr.concat(results1, dim='time')
# %%
st=85
plt.figure()
bright.ImageCalibratedMinusSKImage0.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time").plot.line(hue="im_col")
plt.gca().set_prop_cycle(None)
bright.ImageCalibratedMinusSKImage1.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time").plot.line(hue="im_col",linestyle='dashed')
# %%
st2=85
plt.figure()
bright.ImageCalibratedMinusSKImage0.where((bright.TPsza>st2) & (bright.TPsza<st2+1)).mean("time").plot(robust=True) 
plt.figure()
bright.ImageCalibratedMinusSKImage1.where((bright.TPsza>st2) & (bright.TPsza<st2+1)).mean("time").plot(robust=True) 
# %%
plt.figure()
dark.ImageCalibrated.where((dark.TPsza>106) & (dark.TPsza<107)).mean("time").plot(robust=True) 
plt.figure()
dark.ImageCalibrated.where((dark.TPsza>109) & (dark.TPsza<110)).mean("time").plot(robust=True) 
# %%
plt.figure()
(dark.ImageCalibrated.where((dark.TPsza>109) & (dark.TPsza<110)).mean("time")
 - dark.ImageCalibrated.where((dark.TPsza>106) & (dark.TPsza<107)).mean("time")).plot(robust=True)

# %%78ut
bright=xr.open_dataset('bright')
# %%
plt.figure()
dark.ImageCalibrated.where((dark.TPssa>100) & (dark.TPssa<101)).mean("time").plot(robust=True) 
plt.figure()
dark.ImageCalibrated.where((dark.TPssa>109) & (dark.TPssa<110)).mean("time").plot(robust=True) 
# %%
import matplotlib.animation as animation
ds= dark.where((dark.TPssa>105) & (dark.TPssa<106),drop=True)   
fig, ax = plt.subplots()

def update(frame):
    ax.clear()
    ds.ImageCalibrated.isel(time=frame).plot(ax=ax, robust=False, add_colorbar=False)
    ax.set_title(f"Frame {frame}")

ani = animation.FuncAnimation(fig, update, frames=len(ds.time), repeat=False)
plt.show()
# %%
st=109
ds= dark.where((dark.TPsza>st) & (dark.TPsza<st+1),drop=True) 
for imageno in [10]:
    d = to_datetime(ds.sel(time=ds.time[imageno])["time"].values).tz_localize("UTC").to_pydatetime()
    ts = sfapi.load.timescale()
    t = ts.from_datetime(d)
    ecipos = ds.sel(time=ds.time[imageno]).afsGnssStateJ2000[:3]
    ecefpos = ICRF(Distance(m=ecipos.T).au, t=t, center=399).frame_xyz(itrs).m.T
    geoloc = do_geoloc(ds.isel(time=imageno))
    ecefpos_repeated = np.tile(ecefpos, [63, 9, 1])
    los = geoloc[:, :, 1:] - ecefpos_repeated
    skprofs0, skprofs1 = makeSKprofiles(t.tt - 2400000.5, observer=ecefpos, lookvector=los)
    skimage1 = skprofs1[:, 0].reshape(63, -1)
    skimage0 = skprofs0[:, 0].reshape(63, -1)
    diff_image1 = ds.ImageCalibrated.isel(time=imageno) - skimage1  
    diff_image0 = ds.ImageCalibrated.isel(time=imageno) - skimage0
    results0.append(diff_image0)
    results1.append(diff_image1)
# %%
plt.figure()
ds.ImageCalibrated.isel(time=imageno).plot()
plt.figure()
diff_image0.plot()
plt.figure()
diff_image1.plot()
plt.figure()
diff_image0.plot.line(hue="im_col")
plt.figure()
diff_image1.plot.line(hue="im_col")
# %%
darkbg=dark.ImageCalibrated.where((dark.TPssa>100) & (dark.TPssa<110)).mean("time")
# %%
st=70
plt.figure()
(bright.ImageCalibratedMinusSKImage0.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time")-darkbg).plot.line(hue="im_col")
plt.gca().set_prop_cycle(None)
(bright.ImageCalibratedMinusSKImage1.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time")-darkbg).plot.line(hue="im_col",linestyle='dashed')
# %%
st2=70
plt.figure()
(bright.ImageCalibratedMinusSKImage0.where((bright.TPsza>st2) & (bright.TPsza<st2+1)).mean("time")-darkbg).plot(robust=True) 
plt.figure()
(bright.ImageCalibratedMinusSKImage1.where((bright.TPsza>st2) & (bright.TPsza<st2+1)).mean("time")-darkbg).plot(robust=True) 
# %%
from scipy.stats import linregress
st = 89
bright_diff0 = bright.ImageCalibratedMinusSKImage0.where((bright.TPsza > st) & (bright.TPsza < st + 1)).mean("time") - darkbg
bright_diff1 = bright.ImageCalibratedMinusSKImage1.where((bright.TPsza > st) & (bright.TPsza < st + 1)).mean("time") - darkbg

im_col_range = range(len(bright_diff0.im_col))
slopes = []
intercepts = []

for col in im_col_range:
    y0 = bright_diff0.sel(im_col=col).values[20:40]
    y1 = bright_diff1.sel(im_col=col).values[20:40]
    x = np.arange(20, 40)
    
    slope0, intercept0, _, _, _ = linregress(x, (y0 + y1) / 2)
    
    slopes.append(slope0)
    intercepts.append(intercept0)

print("Slopes for ImageCalibratedMinusSKImage:", slopes)
print("Intercepts for ImageCalibratedMinusSKImage:", intercepts)

# %%

plt.figure()

(bright.ImageCalibratedMinusSKImage0.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time")-darkbg).plot.line(hue="im_col")
plt.gca().set_prop_cycle(None)
(bright.ImageCalibratedMinusSKImage1.where((bright.TPsza>st) & (bright.TPsza<st+1)).mean("time")-darkbg).plot.line(hue="im_col",linestyle='dashed')
plt.gca().set_prop_cycle(None)
for col in im_col_range:
    x = np.arange(len(bright_diff0.im_row))
    fitted_line = slopes[col] * x + intercepts[col]
    plt.gca().plot(x, fitted_line)
plt.show()
# %%
# %%
coeffs = {}

for st in range(70, 90):
    bright_diff0 = bright.ImageCalibratedMinusSKImage0.where((bright.TPsza > st) & (bright.TPsza < st + 1)).mean("time") - darkbg
    bright_diff1 = bright.ImageCalibratedMinusSKImage1.where((bright.TPsza > st) & (bright.TPsza < st + 1)).mean("time") - darkbg

    im_col_range = range(len(bright_diff0.im_col))
    slopes = []
    intercepts = []

    for col in im_col_range:
        y0 = bright_diff0.sel(im_col=col).values[20:40]
        y1 = bright_diff1.sel(im_col=col).values[20:40]
        x = np.arange(20, 40)
        
        slope0, intercept0, _, _, _ = linregress(x, (y0 + y1) / 2)
        
        slopes.append(slope0)
        intercepts.append(intercept0)
    
    coeffs[st] = {'slopes': slopes, 'intercepts': intercepts}

print("Coefficients for each st value:", coeffs)
# %%
coeffs_xr = xr.Dataset(
    {
        'slopes': (('st', 'im_col'), np.array([coeffs[st]['slopes'] for st in coeffs])),
        'intercepts': (('st', 'im_col'), np.array([coeffs[st]['intercepts'] for st in coeffs]))
    },
    coords={
        'st': list(coeffs.keys()),
        'im_col': im_col_range
    }
)

print(coeffs_xr)
# %%
plt.figure()
coeffs_xr.intercepts.plot()
plt.figure()
coeffs_xr.slopes.plot()
# %%
