# %%
# Functions to create overviews for one hemisphere only
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
from Keogram import makeStripMatrix, getTPLatitudes, getSatDates
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import time

# %% Saves keograms for every orbit per day on a pdf page.
def orbit_pdfHS(items, channel, strip_dir, filename, numdays):
    Tperiod = timedelta(minutes=100)
    #one orbit = ca 90 min
    pdf = PdfPages(filename)
    n = 0
    
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        maxorbits = 16 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=1, figsize=(8.27,40), gridspec_kw={'height_ratios': [2]*maxorbits})
        
        #mask_day = (items['EXPDate'] >= pd.to_datetime(starttime + timedelta(days=day-1),utc=True)) & (items['EXPDate'] <= pd.to_datetime(starttime + timedelta(days=day),utc=True))
        #CCD_day = items.loc[mask_day]

        #the first subplot number
        orbnum = 1
        
        #this for loop goes through the images starting from the end of previous orbit
        for i in range(n, len(items)-1):
            #print(i)
            orbit_startdate = items.iloc[n].EXPDate

            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/2: 
                continue
            if deltat > Tperiod/2: #if this is True, next image will belong to next orbit.                         
                #creates orbit from index n to i
                orbit = items.iloc[n:i]
                dates = getSatDates(orbit)
                times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                satlatitudes = getTPLatitudes(orbit)
                #gets the matrix corresponding to that orbit
                matrix = makeStripMatrix(orbit,channel,strip_dir)
                if len(orbit) == 0 :
                    continue
                if orbnum == 1:
                    #plotting latitude vs time for the first orbit
                    axs[0].plot(dates,satlatitudes,'.')
                    axs[0].set_xlabel('Time')
                    axs[0].set_ylabel('Latitude of tangent point')
                    axs[0].set_xlim(dates[0],dates[-1])
                    axs[0].grid(linestyle='-')
                    axs[0].set_title(dates[0].date(), fontsize=16)
                    axs[0].set_xticks(dates[::20])
                    axs[0].set_xticklabels(times_strings[::20], rotation = 30) 
                #plots the orbit found from n to current i
                if channel == 'IR2':
                    axs[orbnum].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                else:
                    axs[orbnum].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True) # rasterized makes a pixel image instead of vector graphic, less saving time
                axs[orbnum].set_title(f"Orbit {orbnum}")
                axs[orbnum].set_xticks(dates[::20])
                axs[orbnum].set_xticklabels(times_strings[::20], rotation = 30) 
                plt.tight_layout(h_pad=1)
                
                orbnum = orbnum + 1
                n = i+1 #start number for next orbit
                nextorbit_startdate = items.iloc[n].EXPDate
                #print(n)
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)

                    break
    
                
        time0 = time.time()   #seems to take around 225 s which is way to slow
        #fig.savefig('plot.png')
        pdf.savefig(fig)
        print(time.time()-time0)
        plt.close(fig)
        
    pdf.close()
    return
 # %% To run the code above
def Main():
    # Determine the main time span and settings for multiple plots
    start_time = DT.datetime(2023,2,15,00,00,0)
    stop_time = DT.datetime(2023,2,22,00,00,0)
    channel = 'IR1'
    strip_dir = 'v'
    filename = "3weekfebIR1_l1b.pdf"
    numdays = stop_time-start_time #number of days

    # puts data in file temp_data
    #df = read_MATS_data(start_time,stop_time,version=0.5,level='1b',filter={"TPlat":[50,90],'NROW': [0,400]})
    #df.to_pickle('3weekfeb')
    #"change latitude filter depending on if you want to look at north or south pole."
    # Reads in data from file
    items = pd.read_pickle('3weekfeb')
    items = items[items['channel'] == channel]

    orbit_pdfHS(items, channel, strip_dir, filename, numdays)

    return


