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

def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    if nheights == None:
        nheights = ccditem['NROW']
    d = ccditem['EXPDate']
    ts = load.timescale()  #loads earth rotation data
    t = ts.from_datetime(d)
    #gets the cameras position and attitude data
    ecipos = ccditem['afsGnssStateJ2000'][0:3] #uses the J2000 equinox epoch,for the global navigation satellite system
    q = ccditem['afsAttitudeState']
    quat = R.from_quat(np.roll(q, -1))
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

#%% Plot a chosen image
def saveSpecIm(ccditem):
    ccdimage = ccditem['ImageCalibrated']
    #saves the matrix to matlab-file
    scipy.io.savemat('aurorapic2',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
    return

# %% 
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

# %%
def IntensityPeak(aurorapeaks):
    collim = 10
    rowlim = 145
    peak_images = []
    sumofpeaks = 0
    
    for strip in aurorapeaks:
 
        im_part = strip.image[rowlim:,collim:]
        im_sum = np.sum(im_part) #sum the top part of image
        peak_images.append(im_sum)
        sumofpeaks = sumofpeaks + im_sum
    return peak_images, sumofpeaks/len(peak_images)
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

# %%  
def orbit_events(items, numdays, Tperiod):
    "gets the aurora peak maximums position for every orbit per day, save in list"
    n = 0
    centercol = 22
    airglowlim = 160
    auroraintensity = 50
    # list of the maximum values (one max per orbit)
    altmaxes = []
    maxlat = []
    maxlon = []
    maxtime = []
    maxtime_strings = []
    maxintensities = []
    aurorastrips = []
    aurorapeaks = []

    #list with all max altitudes (and times) with aurora, not only one per orbit
    allaltitudes = []
    alltimes_strings = []

    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        #this for loop goes through the images starting from the end of previous orbit
        for i in range(n, len(items)-1):
            #print(i)
            orbit_startdate = items.iloc[n].EXPDate

            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/2:
                continue
            if deltat > Tperiod/2:  #if this is True, next image will belong to next orbit.                         
                
                #creates orbit from index n to i
                orbit = items.iloc[n:i]
                
                #Temporary position lists for all the max intensity points in an orbit
                alt_orbit = []
                lat_orbit = []
                lon_orbit = []
                T_orbit = []
                intensities_orbit = []
                
                if len(orbit) == 0 :
                    continue
                else:
                    #time0 = time.time()
                    # This for loop goes through the images belonging to an orbit
                    for n, ccd in orbit.iterrows():
                        #ccdimage = ccd['ImageCalibrated']
                        #ccd_strip= ccdimage[:,centercol]
                        new_strip = CenterStrip(ccd)
                        new_strip.makeVerticalStrip() #creates strip object
                        ccd_strip = new_strip.strip

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        
                        timestamp = ccd.EXPDate

                        if ccd_strip.item(top_max) < auroraintensity:
                            #If no aurora we set the altitude for that strip max to 0
                            allaltitudes.append(0)  
                            alltimes_strings.append(timestamp.strftime("%d/%m %H:%M:%S"))

                        else:
                            #gets the position coordinates of the max intensity point of strips with aurora
                            TPgeo = col_pos(ccd,centercol)
                            [lat,lon,altitude] = TPgeo[row,:]
                            intensity = ccd_strip.item(row)
                            new_strip.maxlat = lat
                            new_strip.maxlon = lon
                            new_strip.maxalt = altitude
                            new_strip.maxI = intensity
                            #print(row, altitude, intensity, timestamp)
                            
                            alt_orbit.append(altitude)
                            lat_orbit.append(lat)
                            lon_orbit.append(lon)
                            T_orbit.append(timestamp)
                            intensities_orbit.append(intensity)
                            #lists of all maximum altitude values of each strip
                            allaltitudes.append(altitude)
                            alltimes_strings.append(timestamp.strftime("%d/%m %H:%M:%S"))

                            aurorastrips.append(new_strip)
                    #print(time.time()-time0)
                    
                    if len(alt_orbit) == 0 :
                        pass
                    else:
                        maxalt = (max(alt_orbit)) #the maximum altitude point of an orbit, will be two points for two peaks
                        #to add the second peak if two exists, use time or distance condition
                        altmaxes.append(maxalt)
                        ind = alt_orbit.index(maxalt)
                            
                        #lists with maximums for each orbit
                        maxlat.append(lat_orbit[ind]) #adds the latitude corresponding to max altitude
                        maxlon.append(lon_orbit[ind]) #adds the longitude corresponding to max altitude
                        maxtime.append(T_orbit[ind])
                        maxtime_strings.append(T_orbit[ind].strftime("%d/%m %H:%M:%S"))
                        maxintensities.append(intensities_orbit[ind])
                        aurorapeaks.append(aurorastrips[ind])
                        
  
                n = i+1 #start number for next orbit
                nextorbit_startdate = items.iloc[n].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)
                    
                    break
        #saving to matlab files
        scipy.io.savemat('altitudesfeb',{'altitudesfeb': allaltitudes, 'label':'altitudes'}) 
        scipy.io.savemat('alltimesfeb',{'alltimesfeb': alltimes_strings, 'label':'times'}) 
        scipy.io.savemat('maxaltfeb',{'maxaltfeb': altmaxes, 'label':'altitudes'}) 
        scipy.io.savemat('maxtimefeb',{'maxtimefeb': maxtime_strings, 'label':'times'}) 
        scipy.io.savemat('maxlatfeb',{'maxlatfeb': maxlat, 'label':'latitudes'}) 
        scipy.io.savemat('maxlonfeb',{'maxlonfeb': maxlon, 'label':'longitudes'})
        scipy.io.savemat('maxintensities',{'maxintensities': maxintensities, 'label':'intensities'})          

    return aurorapeaks

# %%
def Main():
    start_time = DT.datetime(2023,2,15,00,0,0)
    stop_time = DT.datetime(2023,2,16,0,0,0)
    channel = 'IR1'
    numdays = stop_time-start_time #number of days
    Tperiod = timedelta(minutes=100)

    #df = read_MATS_data(start_time,stop_time,filter={"TPlat":[50,90],'NROW': [0,400]},version=0.5,level='1b')
    #df.to_pickle('4weekfeb')
    #"change latitude filter depending on if you want to look at north or south pole."

    items = pd.read_pickle('15feb')
    items = items[items['channel'] == channel]
    #pic = items.iloc[100]
    auroralist = orbit_events(items,numdays,Tperiod)
    #IntensityEvent(auroralist)
    #KeogramAltitudes(items,channel)
    #saveSpecIm(pic)
    return auroralist
# %%
auroralist = Main()
IntensityPeak(auroralist)
# %%
