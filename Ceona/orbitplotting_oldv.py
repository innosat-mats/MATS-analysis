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
start_time = DT.datetime(2023,2,19,00,00,0)
stop_time = DT.datetime(2023,2,20,00,00,0)
channel = 'IR1'
strip_dir = 'v'
filename = "day.pdf"
numdays = stop_time-start_time #number of days

def getSatLatitudes(objects):
    satlat_list= []
    for n, CCD in objects.iterrows():
        satlat_list.append(CCD.satlat)
    return satlat_list


def getSatDates(objects):
    listofdates = []
    for n, row in objects.iterrows():
        listofdates.append(row.EXPDate)
    return listofdates

# %% puts data in file temp_data
df = read_MATS_data(start_time,stop_time,version=0.5,level='1a',filter={"satlat":[40,90]})
df.to_pickle('polar_data')
"change latitude filter depending on if you want to look at north or south pole."

#%% Reads in data from file
items = pd.read_pickle('polar_data')
items = items[items['channel'] == channel]

# %% Saves keograms for every orbit per day on a pdf page.
def orbit_Keogram(items, channel, strip_dir, filename, numdays):
    start = items.iloc[0].EXPDate  #instead of start_time ??
    pdf = PdfPages(filename)
    
    # loop that goes through each day
    for day in range(1,numdays.days+1):
        maxorbits = 15 #assumed maximum possible orbits in one day
        fig, axs = plt.subplots(nrows=maxorbits, ncols=1, figsize=(8.27,20))
        
        #makes a list of the CCD items from the day specified in mask range.
        mask_day = (items['EXPDate'] >= pd.to_datetime(start + timedelta(days=day-1),utc=True)) & (items['EXPDate'] <= pd.to_datetime(start + timedelta(days=day),utc=True))
        CCD_day = items.loc[mask_day]
        print("days",len(CCD_day))
        print(CCD_day.iloc[-1].EXPDate)
        print(CCD_day.iloc[0].EXPDate)
        #tells us if the first pictures on that day increases or decreases in latitude
        direction = np.sign(CCD_day.iloc[1].satlat - CCD_day.iloc[0].satlat)
        n = 0
        #the first plot number
        orbnum = 1

        #the while loop goes through the list of CCDs of one day, as long it is not empty
        while n < len(CCD_day):
            changed = False
            exit = False
            #this for loop goes through the images starting from the end of previous orbit
            for i in range(n, len(CCD_day)-1):
                #print(i)
                #print(CCD_day.iloc[i].satlat)

                #checks the latitude change for each image
                latitude_change = CCD_day.iloc[i+1].satlat-CCD_day.iloc[i].satlat
                if latitude_change == 0:
                    continue
                if not changed and np.sign(latitude_change) != direction:
                    #print('TRUE')
                    changed = True
                elif changed:
                    if np.sign(latitude_change) == direction:
                        print(n,i)
                        #print(CCD_day.iloc[i].EXPDate)
                        orbit = CCD_day.iloc[n:i] 
                        n = i+1
                        exit = True
                        break
            if not exit:
                break
            #plot here
            #print(len(orbit))
            dates = getSatDates(orbit)
            satlatitudes = getSatLatitudes(orbit)
            matrix = makeStripMatrix(orbit,channel,strip_dir)
            
            #times = [dt.strftime("%H:%M:%S") for dt in dates]

            axs[0].set_xlabel('Time')
            axs[0].set_ylabel('Latitude')
            axs[0].plot(dates,satlatitudes,'.')
            axs[0].set_xlim(dates[0],dates[-1])
            axs[0].grid(linestyle='-')
            axs[0].set_title((start + timedelta(days=day-1)).date(), fontsize=16)
            print(orbnum)
            if orbnum >= maxorbits: #handles the case when we have more pictures left, but not enough subplots
                break
            else:
                axs[orbnum].pcolormesh(dates,range(matrix.shape[0]),matrix, vmax=1200)
                axs[orbnum].set_title(f"Orbit {orbnum}")
                axs[orbnum].set_xlabel('Time')
                plt.tight_layout(h_pad=1)
                orbnum = orbnum + 1
            
        plt.gcf().autofmt_xdate()
        pdf.savefig(fig)
        plt.close()
        plt.show()
    pdf.close()
  
    return
# %%
