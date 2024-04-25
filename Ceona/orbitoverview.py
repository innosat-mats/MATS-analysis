# %%
# functions to create overview including plots for both hemispheres.
# and overview with aurora points plotted on the keograms
import datetime as DT
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
from Keogram import makeStripMatrix, getSatDates, get_stripRow
from matplotlib.backends.backend_pdf import PdfPages

# %% Saves keograms for every orbit per day on a pdf page.
def orbit_pdf(items, channel, filename, numdays):
    "Saves keograms for every orbit per day on a pdf page"
    Tperiod = timedelta(minutes=100)
    #one orbit = ca 90 min
    pdf = PdfPages(filename)
    n = 0
    
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        maxorbits = 16 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=2, figsize=(16.54,40), gridspec_kw={'height_ratios': [2]*(maxorbits)})
        #the first subplot number
        orbnum = 1
        subplotNum = 0
        orbcheck = False  #to check if previous orbit was on the same hemisphere as now

        #this for loop goes through the images starting from the end of previous day
        for i in range(n, len(items)-1):
            
            orbit_startdate = items.iloc[n].EXPDate
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6 and i < len(items)-2: 
                continue
            else:                          
                #creates orbit from index n to i
                if items.iloc[i].TPlat > 0: #north hemisphere
                    NH = items.iloc[n:i+1]
                    if len(NH) == 0 :
                        continue
                    dates = getSatDates(NH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix, stripslist = makeStripMatrix(NH)

                    if orbcheck == True:
                        subplotNum += 1
                        orbnum += 1

                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                    else:
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=0, vmax=480) # rasterized makes a pixel image instead of vector graphic, less saving time
                    
                    axs[subplotNum,0].set_title(f"Orbit {orbnum} Northern Hemisphere")
                    axs[subplotNum,0].set_xticks(dates[::20])
                    axs[subplotNum,0].set_xticklabels(times_strings[::20], rotation = 30) 
                    orbcheck = True

                elif items.iloc[i].TPlat < 0: #south hemisphere
                    SH = items.iloc[n:i+1]
                    if len(SH) == 0 :
                        continue
                    dates = getSatDates(SH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix, stripslist = makeStripMatrix(SH)

                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                    else:
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=0, vmax=480) # rasterized makes a pixel image instead of vector graphic, less saving time

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

# %%
def overview_points(items, channel, allrows, filename, numdays):
    "Create the overviews with added red dots for the aurora strips max points"
    Tperiod = timedelta(minutes=100)
    #one orbit = ca 90 min
    pdf = PdfPages(filename)
    n = 0
    
    # loop that goes through number of days
    for day in range(1,numdays.days+1):
        maxorbits = 16 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=2, figsize=(16.54,40), gridspec_kw={'height_ratios': [2]*(maxorbits)})
        #the first subplot number
        orbnum = 1
        subplotNum = 0
        orbcheck = False  #to check if previous orbit was the same as now

        #this for loop goes through the images starting from the end of previous day
        for i in range(n, len(items)-1):
            
            orbit_startdate = items.iloc[n].EXPDate
            
            #checks the time change for each image
            deltat = items.iloc[i+1].EXPDate-items.iloc[i].EXPDate
            if deltat < Tperiod/6 and i < len(items)-2: 
                continue
            else:                          
                #creates orbit from index n to i
                if items.iloc[i].TPlat > 0: #north hemisphere
                    NH = items.iloc[n:i+1]

                    if len(NH) == 0 :
                        continue
                    dates = getSatDates(NH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix, stripslist = makeStripMatrix(NH)
                    if orbcheck == True:
                        subplotNum += 1
                        orbnum += 1

                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                        axs[subplotNum,0].scatter(dates,allrows[n:i+1], marker='.',s=5, color="red")

                    else:
                        axs[subplotNum,0].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=0, vmax=480) # rasterized makes a pixel image instead of vector graphic, less saving time
                        axs[subplotNum,0].scatter(dates,allrows[n:i+1], marker='.',s=5, color="red")

                    axs[subplotNum,0].set_title(f"Orbit {orbnum} Northern Hemisphere")
                    axs[subplotNum,0].set_xticks(dates[::20])
                    axs[subplotNum,0].set_xticklabels(times_strings[::20], rotation = 30) 
                    axs[subplotNum,0].set_ylabel('Row')
                    axs[subplotNum,0].set_xlabel('Time')
                    orbcheck = True

                if items.iloc[i].TPlat < 0: #south hemisphere
                    SH = items.iloc[n:i+1]
                    if len(SH) == 0 :
                        continue
                    dates = getSatDates(SH)
                    times_strings = [dt.strftime("%H:%M") for dt in dates]  #as strings
                    #gets the matrix corresponding to that hemisphere
                    matrix, stripslist = makeStripMatrix(SH)
                    
                    #plots the orbit found from n to current i
                    if channel == 'IR2':
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=-20, vmax=260) # rasterized makes a pixel image instead of vector graphic
                        axs[subplotNum,1].scatter(dates, allrows[n:i+1], marker='.',s = 5, color="red")

                    else:
                        axs[subplotNum,1].pcolormesh(dates,range(matrix.shape[0]),matrix, rasterized = True, vmin=0, vmax=480) # rasterized makes a pixel image instead of vector graphic, less saving time
                        axs[subplotNum,1].scatter(dates,allrows[n:i+1], marker='.', s = 5,color="red")
                        

                    axs[subplotNum,1].set_title(f"Orbit {orbnum} Southern Hemisphere")
                    axs[subplotNum,1].set_xticks(dates[::20])
                    axs[subplotNum,1].set_ylabel('Row')
                    axs[subplotNum,1].set_xlabel('Time')

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
    start_time = DT.datetime(2023,4,1,00,00,0)
    stop_time = DT.datetime(2023,4,8,00,00,0)
    channel = 'IR1'
    filename = "apr1W_IR1.pdf"
    numdays = stop_time-start_time #number of days
    items = pd.read_pickle(r'MatsData\1to7aprIR1')
    #orbit_pdf(items, channel, strip_dir, filename, numdays)
    
    #Run this to read in all the strips, and to get the row parameter for each strip.
    allstrips = pd.read_pickle(r'MatsData\apr1Wallstrips')
    allrows = get_stripRow(allstrips)
    overview_points(items,channel,allrows,filename,numdays)
    return

# %%
