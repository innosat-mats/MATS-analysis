# %%
import datetime as DT
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
from Keogram import makeStripMatrix, getSatDates
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# %% Saves keograms for every orbit per day on a pdf page.
def orbit_pdf(items, channel, strip_dir, filename, numdays):
    Tperiod = timedelta(minutes=100)
    #one orbit = ca 90 min
    pdf = PdfPages(filename)
    n = 0
    
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        maxorbits = 15 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=2, figsize=(16.54,40), gridspec_kw={'height_ratios': [2]*(maxorbits)})
        #the first subplot number
        orbnum = 1
        subplotNum = 0
        orbcheck = False  #to check if previous orbit was the same as now

        #this for loop goes through the images starting from the end of previous orbit
        for i in range(n, len(items)-1):
            
            orbit_startdate = items.iloc[n].EXPDate
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6: 
                continue
            else:                          
                #creates orbit from index n to i
                if items.iloc[i].TPlat > 0: #north hemisphere
                    NH = items.iloc[n:i]
                    if len(NH) == 0 :
                        continue
                    dates = getSatDates(NH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix = makeStripMatrix(NH,channel,strip_dir)

                    if orbcheck == True:
                        subplotNum += 1
                        orbnum += 1
                    print('NH',subplotNum,orbnum, len(NH))

                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                    else:
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True) # rasterized makes a pixel image instead of vector graphic, less saving time
                
                    axs[subplotNum,0].set_title(f"Orbit {orbnum} Northern Hemisphere")
                    axs[subplotNum,0].set_xticks(dates[::20])
                    axs[subplotNum,0].set_xticklabels(times_strings[::20], rotation = 30) 
                    orbcheck = True

                if items.iloc[i].TPlat < 0: #south hemisphere
                    SH = items.iloc[n:i]
                    if len(SH) == 0 :
                        continue
                    dates = getSatDates(SH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix = makeStripMatrix(SH,channel,strip_dir)

                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                    else:
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True) # rasterized makes a pixel image instead of vector graphic, less saving time
                    print('SH',subplotNum,orbnum, len(SH))

                    axs[subplotNum,1].set_title(f"Orbit {orbnum} Southern Hemisphere")
                    axs[subplotNum,1].set_xticks(dates[::20])
                    axs[subplotNum,1].set_xticklabels(times_strings[::20], rotation = 30) 
                    orbnum += 1
                    subplotNum += 1
                    orbcheck = False
                
                suptitle = fig.suptitle(f"{orbit_startdate.date()} Ch IR1", fontsize=16)   
                suptitle.set_y(0.999)
                
                n = i+1 #start number for next orbit, to avoid unnecessary search
                nextorbit_startdate = items.iloc[n].EXPDate
                #comparing the day at start of the new orbit with the active orbits start.
                if orbit_startdate.day != nextorbit_startdate.day:
                    #then we want to quit this for loop and start a new day
                    print("new day")
                    print(orbit_startdate)
                    
                    break

        plt.tight_layout(h_pad=1)
        pdf.savefig(fig)
        plt.close(fig)
        
    pdf.close()
    return pdf
 # %% To run the code above
def Main():
    # Determine the main time span and settings for multiple plots
    start_time = DT.datetime(2023,2,8,00,00,0)
    stop_time = DT.datetime(2023,2,15,00,00,0)
    channel = 'IR1'
    strip_dir = 'v'
    filename = "2weekfebIR1_l1b.pdf"
    numdays = stop_time-start_time #number of days

    items = pd.read_pickle('8to14febIR1')
   
    orbit_pdf(items, channel, strip_dir, filename, numdays)

    return

# %%
