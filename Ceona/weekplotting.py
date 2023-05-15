#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from datetime import timedelta
import pandas as pd 
import matplotlib.pyplot as plt 
from Keogram import getLatitudes, makeStripMatrix
from matplotlib.backends.backend_pdf import PdfPages

# Determine the main time span and settings for multiple plots
start_time = DT.datetime(2023,2,17,00,00,0)
stop_time = DT.datetime(2023,2,24,00,00,0)
tdelta = stop_time-start_time #number of days
weeks = tdelta.days//7
channel = 'IR1'
strip_dir = 'v'
filename = "weekplot.pdf"

# %% puts data in file temp_data
df = read_MATS_data(start_time,stop_time,version=0.5,level='1a',filter={"TPlat":[65,75]})
df.to_pickle('temp_data')
#%% Reads in data from file
items = pd.read_pickle('temp_data')

# %% Creates a figure with keograms for each day
def week_Keogram(items, channel, strip_dir, filename):

    pdf = PdfPages(filename)
    for week in range(weeks):
        fig, axs = plt.subplots(nrows=8, ncols=1, figsize=(8.27,11.69))
        start = start_time
        for day in range(week*7+1,(week+1)*7+1): 
            print(day)
            print(start)
            mask_day = (items['EXPDate'] >= pd.to_datetime(start + timedelta(days=day-1),utc=True)) & (items['EXPDate'] <= pd.to_datetime(start + timedelta(days=day),utc=True))
            #makes a list of the CCD items from the day specified in mask range.
            CCD_day = items.loc[mask_day]
            print(len(CCD_day))
            latitudes, dates = getLatitudes(CCD_day,channel)
            matrix = makeStripMatrix(CCD_day,channel,strip_dir)  

            axs[0].set_xlabel('Time')
            axs[0].set_ylabel('Latitude')
            axs[0].plot(dates,latitudes,'.')
            #axs[0].set_xlim(dates[0],dates[-1])
            axs[0].grid(linestyle='-')
            axs[day].pcolormesh(latitudes,range(matrix.shape[0]),matrix)
            axs[day].set_title(f"Day {day}")
            axs[day].set_xlabel('Latitude') 
            
            #Ändra till latitud begränsning istället för tid på något sätt
            #CCD_nighttime = CCD_day.between_time('19:00','20:00')
        pdf.savefig(fig)
        plt.close()
        start = start + timedelta(days=day)
        pdf.close()
  
    return