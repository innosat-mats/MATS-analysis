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
from Keogram import makeStripMatrix, getSatDates, 
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

#%% Plotting a chosen image
def saveSpecIm(ccditem):
    ccdimage = ccditem['ImageCalibrated']
    #saves the matrix to matlab-file
    scipy.io.savemat('aurorapic2',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
    return

# %% altitudes for all strips and keogram for it, even airglow
def KeogramAltitudes(items,channel):
    centercol = 22
    dates = []
    time_strings = []  #as strings
    altitudes= []
    for n, item in items.iterrows():
        ccdimage = item['ImageCalibrated']
        ccd_strip= ccdimage[:,centercol]  #creates strip object #full image [0:200,0:45]
        "finds the row of the max value of each center strip"
        row = np.argmax(ccd_strip)
        #if only some of the points above certain limit
        #row = np.argmax(ccd_strip[160:]) + 160 
        "calculates altitude of the max value"
        TPgeo = col_pos(item,centercol)
        [lat,lon,altitude] = TPgeo[row,:]
        #print(row,altitude)
        date = item.EXPDate
        dates.append(date)
        time_strings.append(date.strftime("%d/%m %H:%M:%S"))
        altitudes.append(altitude)
        print(row,altitude,date)

    matrix = makeStripMatrix(items,channel)
    scipy.io.savemat('keogram_mat',{'keogram_mat': matrix, 'label':'values'}) #saves to matlabfile
    scipy.io.savemat('alt15feb',{'alt15feb': altitudes, 'label':'altitudes'}) #saves to matlabfile
    scipy.io.savemat('time15feb',{'time15feb': time_strings, 'label':'times'}) #saves to matlabfile

    plt.plot(time_strings,altitudes, '.')
    plt.xticks(time_strings[::int(len(dates)/13)], rotation=30)
    plt.title('15 February')
    plt.ylabel('Altitude (km)')
    plt.grid()

# %%  gets altitude and latitudes of keograms for every orbit per day, save in list
def orbit_events(items, numdays, Tperiod):
    n = 0
    centercol = 22
    airglowlim = 160
    # list of the maximum values (one max per orbit)
    altmaxes = []
    maxlat = []
    maxlon = []
    maxtime = []
    maxtime_strings = []
    maxintensities = []

    #list with all altitudes above threshold
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
                altitudes_orbit = [] #contains altitudes for all the max intensity points in an orbit
                latitudes_orbit = []
                longitudes_orbit = []
                times_orbit = []
                intensities_orbit = []
                
                if len(orbit) == 0 :
                    continue
                else:
                    #time0 = time.time()
                    # for loop goes through the image strips in a orbit
                    for n, ccd in orbit.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip= ccdimage[:,centercol]
                        #finds the row of the max intensity value of each strip
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim                       
                        
                        #if the maximum intensity is lower than this threshold, skip strip
                        if row <= airglowlim + 3:
                            allaltitudes.append(0)
                            timestamp = ccd.EXPDate
                            alltimes_strings.append(timestamp.strftime("%d/%m %H:%M:%S"))
                            continue
                        else:
                        #threshold (a few points from airglow limit), points below will not be checked
                            TPgeo = col_pos(ccd,centercol)
                            [lat,lon,altitude] = TPgeo[row,:]
                            timestamp = ccd.EXPDate
                            #gets the position coordinates of the max intensity point of each strip
                            print(row, altitude,ccdimage.item(row,centercol), timestamp)
                            #Temporary position lists for each orbit
                            altitudes_orbit.append(altitude)
                            latitudes_orbit.append(lat)
                            longitudes_orbit.append(lon)
                            times_orbit.append(timestamp)
                            intensities_orbit.append(ccdimage.item(row,centercol))
                            #lists of all maximum altitude values of each strip
                            allaltitudes.append(altitude)
                            alltimes_strings.append(timestamp.strftime("%d/%m %H:%M:%S"))
                    #print(time.time()-time0)
                    if len(altitudes_orbit) == 0 :
                        pass
                    else:

                        maxalt = (max(altitudes_orbit)) #the maximum altitude point of an orbit, will be two points for two peaks
                        ind = altitudes_orbit.index(maxalt)

                        #lists with maximums for each orbit
                        altmaxes.append(maxalt)
                        maxlat.append(latitudes_orbit[ind]) #adds the latitude corresponding to max altitude
                        maxlon.append(longitudes_orbit[ind]) #adds the longitude corresponding to max altitude
                        maxtime.append(times_orbit[ind])
                        maxtime_strings.append(times_orbit[ind].strftime("%d/%m %H:%M:%S"))
                        maxintensities.append(intensities_orbit[ind])
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

        #plt.plot(maxtime_strings,altmaxes, '.') 
        #plt.xticks(maxtime_strings[::int(len(maxtime)/20)], rotation=30)
        #plt.title('15 february')
        #plt.ylabel('Altitude (km)')
        #plt.grid()             

    return

# %%
def Main():
    start_time = DT.datetime(2023,2,15,00,0,0)
    stop_time = DT.datetime(2023,2,22,00,0,0)
    channel = 'IR1'
    numdays = stop_time-start_time #number of days
    Tperiod = timedelta(minutes=100)

    #df = read_MATS_data(start_time,stop_time,filter={"TPlat":[50,90],'NROW': [0,400]},version=0.5,level='1b')
    #df.to_pickle('4weekfeb')
    #"change latitude filter depending on if you want to look at north or south pole."

    items = pd.read_pickle('3weekfeb')
    items = items[items['channel'] == channel]

    #orbit_events(items,numdays,Tperiod)
    #KeogramAltitudes(items,channel)
    #saveSpecIm(ccditem)
    return