import scipy
from datetime import timedelta
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from aurora_analysis import col_pos, makeStripMatrix, save_strips, set_aurora_spec
from Keogram import CenterStrip, getTPLatitudes, getSatDates

# %% Other extra functions and old aurora finder
def IntensityEvent(aurorastrips):
    "OLD FUNCTION : Adds all aurora image intensity integrations from an orbit. Saves the sums in list"
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
    """Function to get altitudes for all strips and keogram for it, above row 160"""
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

def plotKeogram(df, channels, strip_dir):
    #Given a list of ccds. Returns and saves (matlabfile) 2 plots, one from keogram-matrix and one with the TPlatitudes
    if len(channels)==1 :
        IR_list = df[df['channel'] == channels[0]]
        fig,axs = plt.subplots(nrows=2, ncols=1,sharex='all')
        latitudes = getTPLatitudes(IR_list)
        dates = getSatDates(IR_list)
        times_strings = [dt.strftime("%d/%m %H:%M") for dt in dates]
        matrix = makeStripMatrix(df,channels[0],strip_dir)
        #axs[0].imshow(matrix, origin = 'lower')
        axs[0].pcolormesh(dates,range(matrix.shape[0]),matrix)
        axs[0].set_title(f"28 Feb Channel {channels[0]}")
        axs[0].set_xlim(dates[0],dates[-1])
        axs[0].set_xticks(dates[::int(len(dates)/10)])
        axs[0].set_xticklabels(times_strings[::int(len(times_strings)/10)], rotation = 30)
               
        axs[1].set_xlabel('Time')
        axs[1].set_ylabel('Latitude')
        axs[1].plot(dates,latitudes,'.')
        #axs[1].set_xlim(dates[0],dates[-1])
        axs[1].grid(linestyle='-')
        axs[1].set_xticks(dates[::int(len(dates)/10)])
        axs[1].set_xticklabels(times_strings[::int(len(times_strings)/10)], rotation = 30) 
               
    else:
        fig,axs = plt.subplots(nrows=len(channels)+1, ncols=1)
        IR1_list = df[df['channel'] == channels[0]]
        # gets the latitudes from objects only for the first channel.
        latitudes = getTPLatitudes(IR1_list)
        dates = getSatDates(IR1_list)
        times_strings = [dt.strftime("%d/%m %H:%M") for dt in dates] 

        axs[len(channels)].set_xlabel('Time')
        axs[len(channels)].set_ylabel('Latitude')
        axs[len(channels)].plot(dates,latitudes,'.')
        axs[len(channels)].set_xlim(dates[0],dates[-1])
        axs[len(channels)].set_xticks(dates[::int(len(dates)/10)])
        axs[len(channels)].set_xticklabels(times_strings[::int(len(times_strings)/10)], rotation = 30) 

        for i in range(len(channels)):
            IR_objects =  df[df['channel'] == channels[i]]
            latitudes = getTPLatitudes(IR_objects)
            dates = getSatDates(IR_objects)
            matrix = makeStripMatrix(df,channels[i],strip_dir)  
            axs[i].pcolormesh(dates,range(matrix.shape[0]),matrix, vmax=1500)
            axs[i].set_title(f"Channel {channels[i]}")
       
    #plt.tight_layout()
    #plt.gcf().autofmt_xdate()
    scipy.io.savemat('keogram_mat',{'keogram_mat': matrix, 'label':'values'}) #saves to matlabfile
    scipy.io.savemat('lat28feb',{'lat28feb': latitudes, 'label':'altitudes'}) #saves to matlabfile
    scipy.io.savemat('time28feb',{'time28feb': times_strings, 'label':'times'}) #saves to matlabfile
    plt.show()
