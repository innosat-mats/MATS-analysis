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
from Keogram import makeStripMatrix
import time

def getTPtimes(objects):
    listoftimes = []
    for n, row in objects.iterrows():
        listoftimes.append(row.EXPDate.time())
    return listoftimes

def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    if nheights == None:
        nheights = ccditem['NROW']
    d = ccditem['EXPDate']
    ts = load.timescale()  #earth rotation data
    t = ts.from_datetime(d)
    ecipos = ccditem['afsGnssStateJ2000'][0:3] #uses the J2000 global navigation satellite system
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
        #lat,lon,rad = pos.frame_latlon(itrs)
        position = wgs84.geographic_position_of(posGC)
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
    
# %%
# Determine the main time span
start_time = DT.datetime(2023,2,15,00,0,0)
stop_time = DT.datetime(2023,2,22,00,0,0)
channel = 'IR1'
strip_dir = 'v'
centercol = 22

# %% read data in certain latitude span
df = read_MATS_data(start_time,stop_time,version=0.5,level='1b',filter={"TPlat":[50,90],'NROW': [0,400]})
df.to_pickle('3weekfeb')
"change latitude filter depending on if you want to look at north or south pole."

# %% Reads in data from file
items = pd.read_pickle('3weekfeb')
items = items[items['channel'] == channel]

#%% Plotting a chosen image
R_E = 6356.752  #earth polar radius
ccditem = items.iloc[40]
ccdimage = ccditem['ImageCalibrated']
ccd_strip= ccdimage[:,centercol]  #full image [0:200,0:45]
#finds the max value of the center strip, and the corresponding altitude
maxval = max(ccd_strip)
row = np.where(ccd_strip==maxval)[0][0]
TPgeo = col_pos(ccditem,centercol)
[lat,lon,height] = TPgeo[row,:]
altitude = height-R_E
print(lat,altitude)
#saves the matrix to matlab-file
scipy.io.savemat('aurorapic2',{'ccdimage': ccdimage, 'label':'intensity'}) #saves to matlabfile
#plotting
plt.pcolormesh(ccdimage)
plt.grid(True, which='minor', axis='both', linestyle='-', color='k')


# %% IGNORE latitudes for entire event
lat_list= []
TPlat_list = []
for n, item in items.iterrows():
    ccdimage = item['ImageCalibrated']
    ccd_strip= ccdimage[:,centercol]  #full image [0:200,0:45]
    #finds the max value of the center strip
    maxval = max(ccd_strip)
    row = np.where(ccd_strip==maxval)[0][0]
    #print(maxval, row)
    #calculated the latitude and altitude of the max value in center strip
    TPgeo = col_pos(item,centercol)
    [lat,lon,height] = TPgeo[row,:]
    #print(lat,height-R_E)
    lat_list.append(lat)
    TPlat_list.append(item.TPlat)
plt.plot(lat_list, '.')
plt.plot(TPlat_list,'.')
scipy.io.savemat('lat15feborb10',{'latlist': lat_list, 'label':'latitude'}) #saves to matlabfile

# %% altitudes for one event
times = getTPtimes(items)
times_strings = [dt.strftime("%H:%M:%S") for dt in times]  #as strings
altitudes= []
for n, item in items.iterrows():
    ccdimage = item['ImageCalibrated']
    ccd_strip= ccdimage[:,centercol]  #full image [0:200,0:45]
    #finds the row of the max value of each center strip
    row = np.argmax(ccd_strip[155:]) + 155 
    #calculates altitude of the max value
    TPgeo = col_pos(item,centercol)
    [lat,lon,altitude] = TPgeo[row,:]
    print(row,altitude)

    altitudes.append(altitude)
plt.plot(times_strings,altitudes, '.')
plt.xticks(times_strings[::int(len(times_strings)/10)], rotation=30)
plt.title('15 February Orbit 11')
plt.ylabel('Altitude (km)')
plt.grid()
matrix = makeStripMatrix(items,channel,strip_dir)
scipy.io.savemat('keogram_mat',{'keogram_mat': matrix, 'label':'values'}) #saves to matlabfile
scipy.io.savemat('alt15feb',{'alt15feb': altitudes, 'label':'altitudes'}) #saves to matlabfile
scipy.io.savemat('time15feb',{'time15feb': times_strings, 'label':'times'}) #saves to matlabfile

# %% load settings for orbit_altitudes function
start_time = DT.datetime(2023,2,15,00,0,0)
stop_time = DT.datetime(2023,2,16,00,0,0)
channel = 'IR1'
strip_dir = 'v'
centercol = 22
airglowlim = 155

items = pd.read_pickle('15febIR1')
items = items[items['channel'] == channel]
numdays = stop_time-start_time #number of days
Tperiod = timedelta(minutes=100)

# %%  gets altitude and latitudes of keograms for every orbit per day, save in list
def orbit_altitudes(items, numdays, Tperiod):
    n = 0
    # list of the maximum values (one max per orbit)
    altmaxes = []
    maxlat = []
    maxlon = []
    maxtime = []
    maxtime_strings = []

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
                
                if len(orbit) == 0 :
                    continue
                else:
                    time0 = time.time()
                    # for loops goes through the image strips in a orbit
                    for n, ccd in orbit.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip= ccdimage[:,centercol]
                        #finds the row of the max intensity value of each strip
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim                       
                        
                        if row <= airglowlim + 3:
                            continue
                        else:
                        #kolla bara altituden på minst två eller tre punkter upp,
                        #strunta i att kolla de punkterna under gränsen
                            TPgeo = col_pos(ccd,centercol)
                            [lat,lon,altitude] = TPgeo[row,:]
                            timestamp = ccd.EXPDate
                            #gets the position of the max intensity point of each strip
                            #print(row,altitude)
                            #Temporary position lists for each orbit
                            altitudes_orbit.append(altitude)
                            latitudes_orbit.append(lat)
                            longitudes_orbit.append(lon)
                            times_orbit.append(timestamp)

                            #lists of all maximum altitude values of each strip
                            allaltitudes.append(altitude)
                            alltimes_strings.append(timestamp.strftime("%d/%m %H:%M"))
                    print(time.time()-time0)
                    maxalt = (max(altitudes_orbit)) #the maximum altitude point of an orbit
                    ind = altitudes_orbit.index(maxalt)

                    #lists with maximums for each orbit
                    altmaxes.append(maxalt)
                    maxlat.append(latitudes_orbit[ind]) #adds the latitude corresponding to max altitude
                    maxlon.append(longitudes_orbit[ind]) #adds the longitude corresponding to max altitude
                    maxtime.append(times_orbit[ind])
                    maxtime_strings.append(times_orbit[ind].strftime("%d/%m %H:%M"))
                    
                n = i+1 #start number for next orbit
                nextorbit_startdate = items.iloc[n].EXPDate

                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)
                    break
        #saving to matlab files
        scipy.io.savemat('altitudesfeb',{'altitudes15feb': allaltitudes, 'label':'altitudes'}) 
        scipy.io.savemat('maxalt15feb',{'maxalt15feb': altmaxes, 'label':'altitudes'}) 
        scipy.io.savemat('maxtime15feb',{'maxtime15feb': maxtime_strings, 'label':'times'}) 
        scipy.io.savemat('maxlat15feb',{'maxlat15feb': maxlat, 'label':'latitudes'}) 
        scipy.io.savemat('maxlon15feb',{'maxlon15feb': maxlon, 'label':'longitudes'})

        #plt.plot(maxtime_strings,altmaxes, '.') 
        #plt.xticks(maxtime_strings[::int(len(maxtime)/20)], rotation=30)
        #plt.title('15 february')
        #plt.ylabel('Altitude (km)')
        #plt.grid()             

    return
# %%
