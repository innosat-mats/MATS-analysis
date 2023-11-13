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
# %% Other functions and old aurora finder
def IntensityEvent(aurorastrips):
    "Adds all aurora images together from from each orbit. Saves in list"
    airglowlim = 160
    orbit_peak = []
    start_time = aurorastrips[0].time
    sumorbit = 0

    for strip in aurorastrips:
        
        if strip.time-start_time < timedelta(minutes=7):
            im_part = strip.image[airglowlim:,:]
            im_sum = np.sum(im_part) #sum the top part of image
            sumorbit = sumorbit + im_sum

            if strip == aurorastrips[-1]:
                orbit_peak.append(sumorbit)

        else:
            orbit_peak.append(sumorbit)
            sumorbit = 0
            im_part = strip.image[airglowlim:,:]
            im_sum = np.sum(im_part) #sum the top part of image
            sumorbit = im_sum 
            start_time = strip.time

    #orbit_mean = np.sum(orbit_peak)/len(orbit_peak)
    
    return orbit_peak

def KeogramAltitudes(items,channel):
    """altitudes for all strips and keogram for it, above row 160"""
    centercol = 22
    dates = []
    airglowlim = 160
    time_strings = []  #as strings
    altitudes= []
    for n, item in items.iterrows():
        ccdimage = item['ImageCalibrated']
        ccd_strip= ccdimage[:,centercol]  #creates strip object #full image [0:200,0:45]
        #only some of the points above certain limit
        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
        #calculates altitude of the max value of each column
        TPgeo = col_pos(item,centercol)
        [lat,lon,altitude] = TPgeo[row,:]
        #print(row,altitude)
        date = item.EXPDate
        dates.append(date)
        time_strings.append(date.strftime("%d/%m %H:%M:%S"))
        altitudes.append(altitude)
        #print(row,altitude,date)

    matrix, stripslist = makeStripMatrix(items,channel)
    scipy.io.savemat('keogram_mat',{'keogram_mat': matrix, 'label':'values'}) #saves to matlabfile
    scipy.io.savemat('alt15feb',{'alt15feb': altitudes, 'label':'altitudes'}) #saves to matlabfile
    scipy.io.savemat('time15feb',{'time15feb': time_strings, 'label':'times'}) #saves to matlabfile

    plt.plot(time_strings,altitudes, '.')
    plt.xticks(time_strings[::int(len(dates)/13)], rotation=30)
    plt.title('15 February')
    plt.ylabel('Altitude (km)')
    plt.grid()
"""
def aurora_strips(items, numdays, Tperiod):
    #returns the aurora strips and the peak maximums, save in list
    n = 0
    orb = 0
    centercol = 22
    airglowlim = 160
    auroramean = 50
    aurorastrips = []
    aurorastripsNH = []
    aurorastripsSH = []
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        #this for loop goes through the images starting from the end of previous orbit
        for i in range(n, len(items)-1):
            startday = items.iloc[orb].EXPDate

            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6 and i < len(items)-2:
                continue
            else:  #if this is True, next image will belong to next orbit.                         
                if items.iloc[i].TPlat > 0: #north hemisphere
                    auroraintensity = 55
                    #creates orbit from index n to i
                    NH = items.iloc[orb:i+1]
                    if len(NH) == 0 : #if empty, go to next hemisphere
                        continue
                    # This for loop goes through the images belonging to NH
                    for k, ccd in NH.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim

                        #if row >= airglowlim + 5 and ccd_strip.item[row] > auroraintensity:
                        top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check that aurora is there as well
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        if ccd_strip.item(top_max) > auroraintensity:
                            if top_mean > auroramean:
                                
                                #create strip objects of the aurora columns
                                new_strip = CenterStrip(ccd)
                                new_strip.makeVerticalStrip()
                                ccd_strip = new_strip.strip
                                #sets the position coordinates of the max intensity point of strips with aurora
                                set_aurora_spec(new_strip,ccd,row)

                                #list of aurora strip objects
                                aurorastripsNH.append(new_strip)
                                aurorastrips.append(new_strip)       
                elif items.iloc[i].TPlat < 0: #south hemisphere
                    auroraintensity = 70
                    SH = items.iloc[orb:i+1]
                    if len(SH) == 0 :
                        continue
                    for k, ccd in SH.iterrows():
                        st1 = time.time()
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
                        top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        
                        #Check to not include SAA strips in list of aurora strips
                        if SH.iloc[k].TPlat >= -60 and SH.iloc[k].TPlon >= 0 and SH.iloc[k].TPlon <= 40 or (SH.iloc[k].TPlat >= -60 and SH.iloc[k].TPlon >= 300):
                            pass
                        elif ccd_strip.item(top_max) > auroraintensity:
                            if top_mean > auroramean:
                                #create strip objects of the aurora columns only
                                new_strip = CenterStrip(ccd)
                                new_strip.makeVerticalStrip()

                                #sets the position coordinates of the max intensity point of strips with aurora
                                set_aurora_spec(new_strip,ccd,row)
                                
                                #list of aurora strip objects
                                aurorastripsSH.append(new_strip)
                                aurorastrips.append(new_strip) 
                        print(time.time()-st1)                

                orb = i+1 
                nextday_startdate = items.iloc[orb].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if startday.day != nextday_startdate.day:
                    n = orb
                    #then we want to quit this for loop and start a new day
                    print("new day", nextday_startdate)
                    
                    break
    save_strips(aurorastripsNH,'aurorastripsNH.mat','aurorastripsNH')
    save_strips(aurorastripsSH,'aurorastripsSH.mat','aurorastripsSH')
    save_strips(aurorastrips,'aurorastrips.mat','aurorastrips')

    return aurorastrips
"""
# %% Function for tangent point geodetic coordinates
def TPpos(ccditem):
    """Function giving the GPS TP in geodetic coordinates
    Arguments:
        ccditem or dataframe with the 'afsTangentPointECI'
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
    """Saves the magnetic coordinates of TP separately"""
    MLT = []
    Mlat = []
    for k, ccd in items.iterrows():
        TPlat,TPlon,TPalt = TPpos(ccd)
        mlat, mlon, mlt = get_aacgm_coord(TPlat,TPlon,TPalt,ccd.EXPDate, method='ALLOWTRACE')
        MLT.append(mlt)
        Mlat.append(mlat)
    scipy.io.savemat(filename +'MLT.mat',{filename +'MLT': MLT, 'label':'MLT'}) #saves to matlabfile
    scipy.io.savemat(filename +'MLat.mat',{filename +'MLat': Mlat, 'label':'MLat'}) #saves to matlabfile
    return

def set_strip_spec(strip,ccd):
    """Sets the properties of the non aurora strip objects"""
    TPlat,TPlon,TPalt = TPpos(ccd)
    mlat, mlon, mlt = get_aacgm_coord(TPlat,TPlon,TPalt,ccd.EXPDate, method='ALLOWTRACE')
    strip.maxrow = 107  #tried different values, row 107 gave similar values using colpos as TPpos
    strip.maxalt = TPalt
    strip.MagLT = mlt
    strip.Maglat = mlat
    strip.Maglon = mlon
    return

# %% Aurora analysis functions 
def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    """Returns the geodetic coordinates of a pixel, lat, lon and altitude"""
    if nheights == None:
        nheights = ccditem['NROW']
    d = ccditem['EXPDate']
    ts = load.timescale()  #loads earth rotation data
    t = ts.from_datetime(d)
    #gets the cameras position and attitude data
    ecipos = ccditem['afsGnssStateJ2000'][0:3] #uses the J2000 equinox epoch,for the global navigation satellite system
    q = ccditem['afsAttitudeState']  #Written in Euler parameters
    quat = R.from_quat(np.roll(q, -1))  #the last element is put first
    qprime = R.from_quat(ccditem['qprime'])
    ypixels = np.linspace(0, ccditem['NROW'], nheights)
    TPpos = np.zeros((len(ypixels), 3))
    TPgeo = np.zeros((len(ypixels), 3))
    xdeg, ydeg = pix_deg(ccditem, x, ypixels)
    for iy, y in enumerate(ydeg):
        
        los = R.from_euler('XYZ', [0, y, xdeg], degrees=True).apply([1, 0, 0])
        ecivec = quat.apply(qprime.apply(los))
        res = findtangent(t, ecipos, ecivec)
        TPpos[iy, :] = ecipos+res.x*ecivec
        posGC = Geocentric(position_au=Distance(
        m=TPpos[iy,:]).au, t=t)
        #lat,lon,rad = pos.frame_latlon(itrs) old setting for geocentric position
        position = wgs84.geographic_position_of(posGC)  #gives geodetic coordinates
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

def set_aurora_spec(strip,ccd,row):
    "Sets the properties in the strip objects, including position and intensity"
    centercol = 22
    TPgeo = col_pos(ccd,centercol)
    [lat,lon,altitude] = TPgeo[row,:]
    mlat, mlon, mlt = get_aacgm_coord(lat, lon, altitude, strip.time, method='ALLOWTRACE')
    
    strip.maxrow = row
    strip.maxlat = lat
    strip.maxlon = lon
    strip.maxalt = altitude
    strip.MagLT = mlt
    strip.Maglat = mlat
    strip.Maglon = mlon
    #intensity integration of the peak strip image, if I want only a pixel value, ccd_strip.item(row)
    intensity = IntensityPeak(strip)
    strip.maxI = intensity
    return

def IntensityPeak(aurorastrip):
    "Integrates the intensity of a strips original image"
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
        
        if deltat < timedelta(minutes=4): #Belongs to same cluster
            continue
        else:
            #New cluster, check the peak for the previous cluster
            event = allaltitudes[n:i+1]
            print(n,i,len(event)) 
            if len(event) < 4:
                continue
            
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

def save_strips(strips,filename,structname):
    "Creates a panda object of the strip and saves it to matfile"
    maxalt = []
    maxrow = []
    maxlat = []
    maxlon = []
    times = []
    latitudes = []
    intensities = []
    #magnetic aacgm coordinates
    Mlat = []
    Mlon = []
    MagLT = []

    for strip in strips:
        timestamp = strip.time
        maxrow.append(strip.maxrow)
        maxalt.append(strip.maxalt)
        maxlat.append(strip.maxlat)
        maxlon.append(strip.maxlon)
        times.append(timestamp.strftime("%d/%m %H:%M:%S"))
        intensities.append(strip.maxI)
        latitudes.append(strip.latitude)
        Mlat.append(strip.Maglat)
        Mlon.append(strip.Maglon)
        MagLT.append(strip.MagLT)
    pandastrips = pd.DataFrame({'row': maxrow ,'alt': maxalt, 'maxlat': maxlat, 'maxlon': maxlon, 'maxI': intensities, 'lat' : latitudes ,'time' : times, 'MLT' : MagLT, 'Mlat' : Mlat, 'Mlon' : Mlon})  
    pandastrips.to_pickle(structname)
    scipy.io.savemat(filename, {structname: pandastrips.to_dict('list')})

    return

# %%
