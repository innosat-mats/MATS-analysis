# %%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd 
from datetime import timedelta
import matplotlib.pyplot as plt 
from mats_utils.geolocation.coordinates import TPpos, satpos
from Keogram import makeStripMatrix
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Determine the main time span and settings for multiple plots
start_time = DT.datetime(2023,2,27,00,00,0)
stop_time = DT.datetime(2023,2,28,23,59,0)
channel = 'IR1'
strip_dir = 'v'
filename = "february.pdf"
numdays = stop_time-start_time #number of days
Tperiod = timedelta(minutes=100)

def getTPLatitudes(objects):
    TPlat_list= []
    for n, CCD in objects.iterrows():
        TPlat_list.append(CCD.TPlat)
    return TPlat_list


def getSatDates(objects):
    listofdates = []
    for n, row in objects.iterrows():
        listofdates.append(row.EXPDate)
    return listofdates

# %% puts data in file temp_data
df = read_MATS_data(start_time,stop_time,version=0.5,level='1a',filter={"TPlat":[50,90]})
df.to_pickle('28febdata')
"change latitude filter depending on if you want to look at north or south pole."

#%% Reads in data from file
items = pd.read_pickle('28febdata')
items = items[items['channel'] == channel]

# %% Saves keograms for every orbit per day on a pdf page.
def orbit_pdf(items, channel, strip_dir, filename, numdays, Tperiod):
    daystart = start_time
    starttime = items.iloc[0].EXPDate
    startlat = items.iloc[0].TPlat
    pdf = PdfPages(filename)
   
    # loop that goes through each day
    for day in range(1,numdays.days+1):
        maxorbits = 15 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=1, figsize=(8.27,40), gridspec_kw={'height_ratios': [2]*maxorbits})
        
        #makes a list of the CCD items from the day specified in mask range.
        mask_day = (items['EXPDate'] >= pd.to_datetime(starttime + timedelta(days=day-1),utc=True)) & (items['EXPDate'] <= pd.to_datetime(starttime + timedelta(days=day),utc=True))
        CCD_day = items.loc[mask_day]
        print("len of days",len(CCD_day))
        #print(CCD_day.iloc[0].EXPDate)
        #print(CCD_day.iloc[-1].EXPDate)

        n = 0
        #the first plot number
        orbnum = 1

        #the while loop goes through the list of CCDs of one day, as long it is not empty
        while n < len(CCD_day):
            exit = False
            #this for loop goes through the images starting from the end of previous orbit
            for i in range(n, len(CCD_day)-1):
                #print(i)
                #print(CCD_day.iloc[i].satlat)

                #checks the time change for each image
                deltat = CCD_day.iloc[i+1].EXPDate-CCD_day.iloc[i].EXPDate
                if deltat < Tperiod/2:
                    continue
                if deltat > Tperiod/2:  #if this is True, next image will belong to next orbit. 
                    print(n,i)
                    #print(CCD_day.iloc[i].EXPDate)
                    orbit = CCD_day.iloc[n:i] 
                    n = i+1 #start number for next orbit
                    exit = True
                    break
            if not exit:
                print('slut')
                break

            #plots here
            print(len(orbit))
            dates = getSatDates(orbit)
            times_strings = [dt.strftime("%H:%M:%S") for dt in dates]  #as strings
            satlatitudes = getTPLatitudes(orbit)
            matrix = makeStripMatrix(orbit,channel,strip_dir)

            axs[0].set_xlabel('Time')
            axs[0].set_ylabel('Latitude')
            axs[0].plot(dates,satlatitudes,'.')
            axs[0].set_xlim(dates[0],dates[-1])
            axs[0].grid(linestyle='-')
            axs[0].set_title((starttime + timedelta(days=day-1)).date(), fontsize=16)
            print(orbnum)
           
            axs[orbnum].pcolormesh(dates,range(matrix.shape[0]),matrix)  #vmax = 4500
            axs[orbnum].set_title(f"Orbit {orbnum}")
            #axs[orbnum].set_xticks(dates)
            #axs[orbnum].set_xticklabels(times_strings) 
            # Test adding a comment for git               
            axs[orbnum].set_xlabel('Time')
            plt.tight_layout(h_pad=1)
            orbnum = orbnum + 1

        pdf.savefig(fig)
        plt.close()
        plt.show()
    pdf.close()
  
    return
# %%
