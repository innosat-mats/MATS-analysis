# %%
from mats_utils.geolocation.coordinates import findtangent
import scipy
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation as R
from skyfield.api import load
from mats_l1_processing.pointing import pix_deg
from skyfield.units import Distance
from skyfield.toposlib import wgs84
from skyfield.positionlib import Geocentric
from datetime import timedelta
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from Keogram import makeStripMatrix
from aacgmv2 import get_aacgm_coord
#Shepherd, S. G. (2014), Altitude‐adjusted corrected geomagnetic coordinates: Definition and functional approximations, 
#Journal of Geophysical Research: Space Physics, 119, 7501–7521, doi:10.1002/2014JA020264.
# https://aacgmv2.readthedocs.io/en/latest/usage.html#convert-geographic-magnetic-coordinates

# %% Function for tangent point geodetic coordinates
def TPpos(ccditem):
    """Function giving the GPS TP in geodetic coordinates
    Arguments:
        ccditem or dataframe that has the 'afsTangentPointECI' attribute
    Returns:
        TPlat: latitude of TP (degrees)
        TPlon: longitude of TP (degrees)
        TPalt: Geodetic altitude (km)
    """
    eci= ccditem['afsTangentPointECI']
    d = ccditem['EXPDate']
    ts= load.timescale()
    t = ts.from_datetime(d)
    pos = Geocentric(position_au=Distance(
        m=eci).au, t=t)
    position = wgs84.geographic_position_of(pos)
    TPalt = position.elevation.km
    TPlat = position.latitude.degrees
    TPlon = position.longitude.degrees
    return TPlat,TPlon,TPalt

def save_TPMLT(items, filename):
    """Saves the magnetic coordinates of TP separately
        """
    MLT = []
    Mlat = []
    for k, item in items.iterrows():
        TPlat,TPlon,TPalt = TPpos(item)
        mlat, mlon, mlt = get_aacgm_coord(TPlat,TPlon,TPalt,item.EXPDate, method='ALLOWTRACE')
        MLT.append(mlt)
        Mlat.append(mlat)
    scipy.io.savemat(filename +'MLT.mat',{filename +'MLT': MLT, 'label':'MLT'}) #saves to matlabfile
    scipy.io.savemat(filename +'MLat.mat',{filename +'MLat': Mlat, 'label':'MLat'}) #saves to matlabfile
    return

def set_strip_spec(strip,item):
    """Sets the properties of the non aurora strip objects"""
    TPlat,TPlon,TPalt = TPpos(item)
    mlat, mlon, mlt = get_aacgm_coord(TPlat,TPlon,TPalt,item.EXPDate, method='ALLOWTRACE')
    strip.maxrow = 107  #tried different values, row 107 gave similar values using colpos as TPpos
    strip.maxalt = TPalt
    strip.MagLT = mlt
    strip.Maglat = mlat
    strip.Maglon = mlon
    return

# %% Aurora analysis functions 
def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    """Returns the geodetic coordinates of a pixel
    Arguments:
        dataframe (raw MATS-data), x = column index
    Returns: TPgeo
        TPgeo[iy,0] = lat.degrees
        TPgeo[iy,1] = lon.degrees 
        TPgeo[iy,2] = alt.km
    """
    if nheights == None:
        nheights = ccditem['NROW']
    d = ccditem['EXPDate']
    ts = load.timescale()  #loads earth rotation data
    t = ts.from_datetime(d)
    #gets the cameras position and attitude data
    ecipos = ccditem['afsGnssStateJ2000'][0:3] 
    #uses the J2000 equinox epoch, for the global navigation satellite system GNSS
    # afs = Atomic Frequency Standards
    q = ccditem['afsAttitudeState']  #Written in Euler parameters
    quat = R.from_quat(np.roll(q, -1))  #the last element is put first
    qprime = R.from_quat(ccditem['qprime'])
    ypixels = np.linspace(0, ccditem['NROW'], nheights)
    TPpos = np.zeros((len(ypixels), 3))
    TPgeo = np.zeros((len(ypixels), 3))
    xdeg, ydeg = pix_deg(ccditem, x, ypixels) #Function to get the x and y angle from a pixel relative to the center of the CCD
    for iy, y in enumerate(ydeg):
        
        los = R.from_euler('XYZ', [0, y, xdeg], degrees=True).apply([1, 0, 0])
        ecivec = quat.apply(qprime.apply(los))
        res = findtangent(t, ecipos, ecivec)
        TPpos[iy, :] = ecipos+res.x*ecivec
        posGC = Geocentric(position_au=Distance(
        m=TPpos[iy,:]).au, t=t)
        #lat,lon,rad = pos.frame_latlon(itrs) old setting for geocentric position
        position = wgs84.geographic_position_of(posGC)  #gives geodetic coordinates from geocentric position
        alt = position.elevation
        lat = position.latitude
        lon = position.longitude
        TPgeo[iy,0] = lat.degrees
        TPgeo[iy,1] = lon.degrees 
        TPgeo[iy,2] = alt.km
    if splineTPgeo:
        return CubicSpline(ypixels, TPgeo)
    else:
        return TPgeo 

def set_aurora_spec(strip,item,row):
    "Sets the properties in the strip objects, including position and intensity"
    centercol = 22
    TPgeo = col_pos(item,centercol)
    [lat,lon,altitude] = TPgeo[row,:]
    mlat, mlon, mlt = get_aacgm_coord(lat, lon, altitude, strip.time, method='ALLOWTRACE')
    strip.maxrow = row
    strip.maxlat = lat #peak point
    strip.maxlon = lon #peak point
    strip.maxalt = altitude #peak point
    strip.MagLT = mlt
    strip.Maglat = mlat
    strip.Maglon = mlon
    #intensity integration of the peak strip image, if I want only a pixel value, ccd_strip.item(row)
    intensity = IntensityPeak(strip)
    strip.totI = intensity
    return

def IntensityPeak(aurorastrip):
    "Integrates the intensity of a strips full image"
    collow = 14
    coltop = 30
    
    rowlow = aurorastrip.maxrow - 4
    rowtop = aurorastrip.maxrow + 4
    if rowtop >= len(aurorastrip.image):
        im_part = aurorastrip.image[rowlow:len(aurorastrip.image),collow:coltop]
    else:
        im_part = aurorastrip.image[rowlow:rowtop,collow:coltop]
    im_sum = np.sum(im_part) #sum the surrounding part of the peak
    return im_sum

def get_all_altitudes(strips):
    """Arguments: List of strips (class CenterStrip)
    Returns: list with all altitudes from given list"""
    allaltitudes = []
    for index, strip in enumerate(strips):
        altitude = strip.maxalt
        allaltitudes.append(altitude)
    return allaltitudes

def get_aurora_max(aurorastrips,filedate):
    """Function gets the peak point from each aurora cluster
    Arguments: list of aurora strips (class CenterStrip) and name of file Ex: 3WFeb
    Returns: List of peak points, and list of NH peaks and SH peaks"""
    peak_strips = [] # list of the peak point strips for all events.
    peak_stripsNH = []
    peak_stripsSH = []
    allaltitudes = get_all_altitudes(aurorastrips)
    n = 0                    
    for i in range(n,len(aurorastrips)-1):
        strip = aurorastrips[i]
        nextstrip = aurorastrips[i+1]
        deltat = nextstrip.time-strip.time
        
        if deltat < timedelta(minutes=4) and i != len(aurorastrips)-2: #Belongs to same cluster
            continue
        else:
            #New full cluster, check the peak for that cluster
            event = allaltitudes[n:i+1]
            if len(event) < 4:
                continue
            #finds the point of aurora cluster with highest altitude
            peak = max(aurorastrips[n:i+1],key=lambda x: x.maxalt)
            #ind = np.argmax(event)
            #pixelI = aurorastrips[n+ind].strip
            #rad = aurorastrips[n+ind].maxrow
            #print(rad, pixelI[rad],aurorastrips[n+ind].time)
            if peak.maxlat > 0 :
                peak_stripsNH.append(peak)
            else:
                peak_stripsSH.append(peak)
            peak_strips.append(peak)
            n = i+1
    save_strips(peak_strips,filedate +'peaks.mat',filedate +'peaks')
    save_strips(peak_stripsNH,filedate +'peaksNH.mat',filedate +'peaksNH')
    save_strips(peak_stripsSH,filedate +'peaksSH.mat',filedate +'peaksSH')

    return

def save_strips(strips,filedate,structname):
    "Creates a panda object of the strip list and saves it to matfile"
    fullIMG = [] #newly added, has not been re-run with all data
    maxalt = [] #altitude of max intensity point
    maxrow = [] #row of max intensity point
    maxlat = [] #geodetic latitude of max intensity point
    maxlon = [] #geodetic longitude of max intensity point
    times = []
    latitudes = []   #tangent point latitudes
    intensities = []  #full image integrated intensities
    #magnetic aacgm coordinates of the max intensity point
    Mlat = []
    Mlon = []
    MagLT = []

    for strip in strips:
        fullIMG.append(strip.image)
        timestamp = strip.time
        maxrow.append(strip.maxrow)
        maxalt.append(strip.maxalt)
        maxlat.append(strip.maxlat)
        maxlon.append(strip.maxlon)
        times.append(timestamp.strftime("%d/%m %H:%M:%S"))
        intensities.append(strip.totI)
        latitudes.append(strip.latitude)
        Mlat.append(strip.Maglat)
        Mlon.append(strip.Maglon)
        MagLT.append(strip.MagLT)
    pandastrips = pd.DataFrame({'image':fullIMG,'row': maxrow ,'alt': maxalt, 'maxlat': maxlat, 'maxlon': maxlon, 'maxI': intensities, 'lat' : latitudes ,'time' : times, 'MLT' : MagLT, 'Mlat' : Mlat, 'Mlon' : Mlon})
    pandastrips.to_pickle('MatsData/'+ structname)  
    scipy.io.savemat(filedate, {structname: pandastrips.to_dict('list')})
    return

# %%
