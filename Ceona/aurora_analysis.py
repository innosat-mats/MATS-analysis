# %%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import findtangent
import scipy
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation as R
from skyfield.api import load
from mats_l1_processing.pointing import pix_deg
from skyfield.units import Distance
from skyfield.framelib import itrs
from skyfield.toposlib import wgs84
from skyfield.positionlib import Geocentric, ICRF
import datetime as DT
from datetime import datetime, timedelta, timezone
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from Keogram import makeStripMatrix, CenterStrip
import time

# %%
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
    "altitudes for all strips and keogram for it, above row 160"
    centercol = 22
    dates = []
    time_strings = []  #as strings
    altitudes= []
    for n, item in items.iterrows():
        ccdimage = item['ImageCalibrated']
        ccd_strip= ccdimage[:,centercol]  #creates strip object #full image [0:200,0:45]
        #only some of the points above certain limit
        row = np.argmax(ccd_strip[160:]) + 160 
        "calculates altitude of the max value"
        TPgeo = col_pos(item,centercol)
        [lat,lon,altitude] = TPgeo[row,:]
        #print(row,altitude)
        date = item.EXPDate
        dates.append(date)
        time_strings.append(date.strftime("%d/%m %H:%M:%S"))
        altitudes.append(altitude)
        #print(row,altitude,date)

    matrix = makeStripMatrix(items,channel)
    scipy.io.savemat('keogram_mat',{'keogram_mat': matrix, 'label':'values'}) #saves to matlabfile
    scipy.io.savemat('alt15feb',{'alt15feb': altitudes, 'label':'altitudes'}) #saves to matlabfile
    scipy.io.savemat('time15feb',{'time15feb': time_strings, 'label':'times'}) #saves to matlabfile

    plt.plot(time_strings,altitudes, '.')
    plt.xticks(time_strings[::int(len(dates)/13)], rotation=30)
    plt.title('15 February')
    plt.ylabel('Altitude (km)')
    plt.grid()

# %% Used to get the geodetic coordinates of a pixel
def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    "Returns the geodetic coordinates of a pixel, lat, lon and altitude"
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

def set_aurora_spec(strip,ccd,row,centercol):
    "Sets the properties in the strip objects, including position and intensity"
    TPgeo = col_pos(ccd,centercol)
    [lat,lon,altitude] = TPgeo[row,:]
    strip.maxrow = row
    strip.maxlat = lat
    strip.maxlon = lon
    strip.maxalt = altitude
    intensity = IntensityPeak(strip)
    #intensity integration of the peak strip image, if I want only a pixel value, ccd_strip.item(row)
    strip.maxI = intensity
    return [lat,lon,altitude,intensity]

def IntensityPeak(aurorastrip):
    "Integrate the intensity of a strips original image"
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

# %%  
def aurora_strips(items, numdays, Tperiod):
    "returns the aurora strips and the peak maximums, save in list"
    n = 0
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
            orbit_startdate = items.iloc[n].EXPDate
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6:
                continue
            else:  #if this is True, next image will belong to next orbit.                         
                
                if items.iloc[i].TPlat > 0: #north hemisphere
                    auroraintensity = 55
                    #creates orbit from index n to i
                    NH = items.iloc[n:i]
                    if len(NH) == 0 : #if empty, go to next hemisphere
                        continue

                    # This for loop goes through the images belonging to NH
                    for n, ccd in NH.iterrows():
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
                                set_aurora_spec(new_strip,ccd,row,centercol)

                                #list of aurora strip objects
                                aurorastripsNH.append(new_strip)
                                aurorastrips.append(new_strip)       

                if items.iloc[i].TPlat < 0: #south hemisphere
                    if items.iloc[i].TPlat > -50 and items.iloc[i].TPlon > -90 and items.iloc[i].TPlon < 40:
                        continue
                    
                    auroraintensity = 70
                    SH = items.iloc[n:i]
                    if len(SH) == 0 :
                        continue
                    for n, ccd in SH.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
                        top_mean = np.sum(ccd_strip[airglowlim+10:])/len(ccd_strip[airglowlim+10:])

                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        
                        if ccd_strip.item(top_max) > auroraintensity:
                            if top_mean > auroramean:
                                #create strip objects of the aurora columns
                                new_strip = CenterStrip(ccd)
                                new_strip.makeVerticalStrip()
                                ccd_strip = new_strip.strip

                                #sets the position coordinates of the max intensity point of strips with aurora
                                [lat,lon,altitude,intensity] = set_aurora_spec(new_strip,ccd,row,centercol)
                                #timestamp = ccd.EXPDate
                                #print(top_max,ccd_strip.item(top_max), top_mean)
                                #print(row, altitude, timestamp)

                                #list of aurora strip objects
                                aurorastripsSH.append(new_strip)
                                aurorastrips.append(new_strip) 

                n = i+1 #start number for next orbit
                nextorbit_startdate = items.iloc[n].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day", orbit_startdate)
                    
                    break
    save_strips(aurorastripsNH,'aurorastripsNH.mat','aurorastripsNH')
    save_strips(aurorastripsSH,'aurorastripsSH.mat','aurorastripsSH')
    save_strips(aurorastrips,'aurorastrips.mat','aurorastrips')

    return aurorastrips

# %% 
def get_all_altitudes(aurorastrips):
    "Get all altitudes from given list of aurora strip objects"
    allaltitudes = []
    for index, strip in enumerate(aurorastrips):
        altitude = strip.maxalt
        allaltitudes.append(altitude)
    return allaltitudes

def get_aurora_max(aurorastrips):
    "Get the peak points of each aurora cluster"
    peak_strips = [] # list of the peak point strips for all events.
    allaltitudes = get_all_altitudes(aurorastrips)
    n = 0                    
    for i in range(n,len(aurorastrips)-1):
        strip = aurorastrips[i]
        deltat = aurorastrips[i+1].time-strip.time
        
        if deltat < timedelta(minutes=5): #time between events
            continue
        else:
            event = allaltitudes[n:i] 
            if len(event) == 0:
                continue
            
            peak = max(aurorastrips[n:i],key=lambda x: x.maxalt)
            #ind = np.argmax(event)
            #pixelI = aurorastrips[n+ind].strip
            #rad = aurorastrips[n+ind].maxrow
            #print(rad, pixelI[rad],aurorastrips[n+ind].time)
            
            peak_strips.append(peak)
            n = i+1
    save_strips(peak_strips,'peakstrips.mat','peakstrips')
    
    return

def save_strips(aurorastrips,filename,structname):
    maxalt = []
    maxrow = []
    maxlat = []
    maxlon = []
    times = []
    intensities = []
    for strip in aurorastrips:
        maxrow.append(strip.maxrow)
        maxalt.append(strip.maxalt)
        maxlat.append(strip.maxlat)
        maxlon.append(strip.maxlon)
        timestamp = strip.time
        times.append(timestamp.strftime("%d/%m %H:%M:%S"))
        intensities.append(strip.maxI)

    peakdf = pd.DataFrame({'row': maxrow ,'alt': maxalt, 'lat': maxlat, 'lon': maxlon, 'maxI': intensities, 'time' : times})  
    scipy.io.savemat(filename, {structname: peakdf.to_dict('list')})

    return
# %%
def Main(items):
    start_time = DT.datetime(2023,2,15,0,0,0)
    stop_time = DT.datetime(2023,2,17,0,0,0)
    numdays = stop_time-start_time #number of days
    Tperiod = timedelta(minutes=100)
    auroralist = aurora_strips(items,numdays,Tperiod)
    get_aurora_max(auroralist)
    return 

# %%
items = pd.read_pickle('15to16febIR1')

# %%

