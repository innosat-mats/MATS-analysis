# %%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import findtangent
import scipy
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation as R
from skyfield.api import load
from mats_l1_processing.pointing import pix_deg
#from orbitsperday import getSatDates
from skyfield.units import Distance
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric, ICRF
import datetime as DT
from datetime import datetime, timedelta, timezone
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot, all_channels_plot
from Keogram import makeStripMatrix

def getTPtimes(objects):
    listoftimes = []
    for n, row in objects.iterrows():
        listoftimes.append(row.EXPDate.time())
    return listoftimes

def col_pos(ccditem, x, nheights=None, splineTPgeo=False):
    if nheights == None:
        nheights = ccditem['NROW']
    d = ccditem['EXPDate']
    ts = load.timescale()
    t = ts.from_datetime(d)
    ecipos = ccditem['afsGnssStateJ2000'][0:3]
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
        pos = Geocentric(position_au=Distance(
        m=TPpos[iy,:]).au, t=t)
        lat,lon,height = pos.frame_latlon(itrs)
        TPgeo[iy,0] = lat.degrees
        TPgeo[iy,1] = lon.degrees 
        TPgeo[iy,2] = height.km 
    if splineTPgeo:
        return CubicSpline(ypixels, TPgeo)
    else:
        return TPgeo 
    
# %%
# Determine the main time span
start_time = DT.datetime(2023,2,1,00,0,0)
stop_time = DT.datetime(2023,2,3,00,0,0)
channel = 'IR1'
strip_dir = 'v'
centercol = 22
filename = "15feborb11.pdf"

# %% read data in certain latitude span
df = read_MATS_data(start_time,stop_time,version=0.5,level='1b',filter={"TPlat":[50,90],'NROW': [0,400]})
df.to_pickle('1to2febIR1')
"change latitude filter depending on if you want to look at north or south pole."

# %% Reads in data from file
items = pd.read_pickle('1to2febIR1')
items = items[items['channel'] == channel]
times = getTPtimes(items)
times_strings = [dt.strftime("%H:%M:%S") for dt in times]  #as strings
R_E = 6356.752  #earth polar radius

#%% Plotting a chosen image
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


#%% latitudes for entire event
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

#%% altitudes for entire event
altitudes= []
for n, item in items.iterrows():
    ccdimage = item['ImageCalibrated']
    ccd_strip= ccdimage[:,centercol]  #full image [0:200,0:45]
    #finds the row of the max value of each center strip
    maxval = max(ccd_strip)
    row = np.where(ccd_strip==maxval)[0][0]
    #print(maxval, row)
    #calculates altitude of the max value
    TPgeo = col_pos(item,centercol)
    [lat,lon,height] = TPgeo[row,:]
    altitude = height-R_E
    altitudes.append(altitude)
    

# %% gets altitude and latitudes of keograms for every orbit per day, save in list
numdays = stop_time-start_time #number of days
Tperiod = timedelta(minutes=100)
def orbit_keograms(items, channel, strip_dir, numdays, Tperiod):
    n = 0
    altitudes_period = []
    maxlat = []
    timespan = []
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        maxorbits = 16 #assumed maximum possible orbits in one day

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
                altitudes_orbit = []
                latitudes_orbit = []
                times = getTPtimes(orbit)
                if len(orbit) == 0 :
                    continue
                else:
                    for n, ccd in orbit.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip= ccdimage[:,centercol]
                        #finds the row of the max intensity value of each strip
                        maxval = max(ccd_strip)
                        row = np.where(ccd_strip==maxval)[0][0]
                        print(maxval, row)
                        #gets altitude of the max value
                        TPgeo = col_pos(ccd,centercol)
                        [lat,lon,height] = TPgeo[row,:]
                        altitude = height-R_E
                        altitudes_orbit.append(altitude) #contains altitudes for each strip's max intensity point
                        latitudes_orbit.append(lat)
                maxalt = (max(altitudes_orbit)) #the maximum altitude point of an orbit
                ind = altitudes_orbit.index(maxalt)
                maxlat.append(latitudes_orbit[ind]) #gets the latitude corresponding to max altitude
                altitudes_period.append(altitudes_orbit)
                timespan.append(times)
                n = i+1 #start number for next orbit
                nextorbit_startdate = items.iloc[n].EXPDate
                #print(n)
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)
                    break
        plt.plot(altitudes_period,timespan, '.')            
    return


#%%plotting
plt.plot(times_strings,altitudes, '.')
plt.xticks(times_strings[::8], rotation=30)
plt.title('15 February Orbit 11')
plt.ylabel('Altitude (km)')
plt.grid()
scipy.io.savemat('alt15feborb11',{'alt15feborb11': altitudes, 'label':'altitudes'}) #saves to matlabfile
scipy.io.savemat('time15feborb11',{'time15feborb11': times_strings, 'label':'times'}) #saves to matlabfile


# %%
