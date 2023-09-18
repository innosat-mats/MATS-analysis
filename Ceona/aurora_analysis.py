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
    intensity = IntensityPeak(strip)
    #intensity integration of the peak strip image, if I want only a pixel value, ccd_strip.item(row)
    TPgeo = col_pos(ccd,centercol)
    [lat,lon,altitude] = TPgeo[row,:]
    strip.maxrow = row
    strip.maxlat = lat
    strip.maxlon = lon
    strip.maxalt = altitude
    strip.maxI = intensity
    return [lat,lon,altitude,intensity]

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
def IntensityPeak(aurorastrip):
    collow = 14
    coltop = 30
    
    rowlow = aurorastrip.maxrow - 4
    rowtop = aurorastrip.maxrow + 4
    if rowtop >= 180:
        im_part = aurorastrip.image[rowlow:180,collow:coltop]
    else:
        im_part = aurorastrip.image[rowlow:rowtop,collow:coltop]
    im_sum = np.sum(im_part) #sum the surrounding part of the peak
    return im_sum
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
def aurora_strips(items, numdays, Tperiod):
    "gets the aurora strips and the peak maximums, save in list"
    n = 0
    centercol = 22
    airglowlim = 160
    auroraintensity = 50
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
                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        

                        if ccd_strip.item(top_max) > auroraintensity:
                            #create strip objects of the aurora columns
                            new_strip = CenterStrip(ccd)
                            new_strip.makeVerticalStrip()
                            ccd_strip = new_strip.strip

                            #sets the position coordinates of the max intensity point of strips with aurora
                            [lat,lon,altitude,intensity] = set_aurora_spec(new_strip,ccd,row,centercol)
                            timestamp = ccd.EXPDate
                            #print(row, altitude, intensity, timestamp)

                            #list of aurora strip objects
                            aurorastripsNH.append(new_strip)
                            aurorastrips.append(new_strip)       

                if items.iloc[i].TPlat < 0: #south hemisphere
                    SH = items.iloc[n:i]
                    if len(SH) == 0 :
                        continue
                    for n, ccd in SH.iterrows():
                        ccdimage = ccd['ImageCalibrated']
                        ccd_strip = ccdimage[:,centercol]

                        #finds the row of the max intensity value of each strip, above airglow limit
                        row = np.argmax(ccd_strip[airglowlim:]) + airglowlim
                        #gives the row of the maximum 10 rows above the limit to check for aurora there.
                        top_max = np.argmax(ccd_strip[airglowlim+10:]) + airglowlim + 10        
                        
                        if ccd_strip.item(top_max) > auroraintensity:
                            #create strip objects of the aurora columns
                            new_strip = CenterStrip(ccd)
                            new_strip.makeVerticalStrip()
                            ccd_strip = new_strip.strip

                            #sets the position coordinates of the max intensity point of strips with aurora
                            [lat,lon,altitude,intensity] = set_aurora_spec(new_strip,ccd,row,centercol)
                            timestamp = ccd.EXPDate
                            #print(row, altitude, intensity, timestamp)

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
        
    return aurorastrips

# %% Get all altitudes from aurora strips
def get_all_altitudes(aurorastrips):
    #list with all max altitudes (and times) with aurora, not only one per orbit
    allaltitudes = []
    alltimes_strings = []
    for strip in aurorastrips:
        altitude = strip.maxalt
        timestamp = strip.time
        allaltitudes.append(altitude)
        alltimes_strings.append(timestamp.strftime("%d/%m %H:%M:%S"))
    scipy.io.savemat('altitudesfeb',{'altitudesfeb': allaltitudes, 'label':'altitudes'}) 
    scipy.io.savemat('alltimesfeb',{'alltimesfeb': alltimes_strings, 'label':'times'})                               
    return allaltitudes

# %% Not ready to use yet, goal is to get only the peak points
def get_aurora_max(aurorastrips):
    # list of the maximum values (one max per orbit)
    maxalt = []
    maxlat = []
    maxlon = []
    maxtime = []
    maxtime_strings = []
    maxintensities = []
    allaltitude = get_all_altitudes(aurorastrips)
    "Gets the peak points"                    
    for strip in aurorastrips:
        #sorted_indices = np.argsort(alt_orbit)[::-1]
        #ind1 = sorted_indices[0]
        #ind2 = sorted_indices[1]
        #strip.time
        #check so that the second maximum is not the same crossing
        if abs(Time[ind1]-Time[ind2]) > timedelta(minutes=6):  
            maxalt = [max1,max2]
            #lists with 1 or 2 maximums for each orbit
            maxlat.append([maxlat1,maxlat2]) #adds the latitude corresponding to max altitude
            maxlon.append([maxlon1,maxlon2]) #adds the longitude corresponding to max altitude
            maxtime.append([Time[ind1],Time[ind2]])
            maxtime_strings.append([Time[ind1].strftime("%d/%m %H:%M:%S"),Time[ind2].strftime("%d/%m %H:%M:%S")])
            maxintensities.append([intensities[ind1],intensities[ind2]])
            aurorapeaks.append([aurorastrips[ind1],aurorastrips[ind2]])

        else:
            maxalt = [maxalt1,0]
            maxlat.append([maxlat1,0]) #adds the latitude corresponding to max altitude
            maxlon.append([maxlon1,0]) #adds the longitude corresponding to max altitude
            maxtime.append([Time[ind1],0])
            maxtime_strings.append([Time[ind1].strftime("%d/%m %H:%M:%S"),0])
            maxintensities.append([intensities[ind1],0])
            aurorapeaks.append([aurorastrips[ind1],0])
        #to add the second peak if two exists, use time or distance condition
        altmaxes.append(maxalt)    
    #saving to matlab files
    scipy.io.savemat('maxaltfeb',{'maxaltfeb': maxalt, 'label':'altitudes'}) 
    scipy.io.savemat('maxtimefeb',{'maxtimefeb': maxtime_strings, 'label':'times'})         scipy.io.savemat('maxlatfeb',{'maxlatfeb': maxlat, 'label':'latitudes'}) 
    scipy.io.savemat('maxlonfeb',{'maxlonfeb': maxlon, 'label':'longitudes'})
    scipy.io.savemat('maxintensities',{'maxintensities': maxintensities, 'label':'intensities'})          

    return
# %%
def Main():
    start_time = DT.datetime(2023,2,8,0,0,0)
    stop_time = DT.datetime(2023,2,15,0,0,0)
    channel = 'IR1'
    numdays = stop_time-start_time #number of days
    Tperiod = timedelta(minutes=100)

    #df = read_MATS_data(start_time,stop_time,filter={"TPlat":[50,90],'NROW': [0,400]},version=0.5,level='1b')
    #df.to_pickle('4weekfeb')
    #"change latitude filter depending on if you want to look at north or south pole."

    items = pd.read_pickle('8to14febIR1')
    items = items[items['channel'] == channel]
    #pic = items.iloc[100]
    auroralist = aurora_strips(items,numdays,Tperiod)
    #IntensityEvent(auroralist)
    #KeogramAltitudes(items,channel)
    #saveSpecIm(pic)
    return print(len(auroralist))
# %%
auroralist = Main()

# %%
