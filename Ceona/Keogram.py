#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from datetime import timedelta
import pandas as pd 
from pandas import Timestamp
import sys
import numpy as np
import math
import matplotlib.pyplot as plt 
from mats_utils.geolocation.coordinates import TPpos

#%%
class CenterStrip:
    def __init__(self, CCDobject):
        self.image = CCDobject['IMAGE']
        self.strip = []
        self.latitude = TPpos(CCDobject)[0]  #the first position of TPpos gives the latitude
        self.time =  pd.to_datetime(CCDobject['EXPDate'])

    "Makes a strip object from the image"
    def makeVerticalStrip(self):
        "finds the center pixel"
        center = math.ceil(len(self.image[0])/2)
        self.strip = self.image[:,int(center)]
        return  self.strip
    
    def makeHorizontalStrip(self):
        center = math.ceil(len(self.image)/2)
        self.strip = self.image[int(center),:]
        return  np.transpose(self.strip)

#%%  Make a keogram of a specific channel and CCD-objects
def makeStripMatrix(df, channel_type, strip_dir):

    IR_list = df[df['channel'] == channel_type]
    strips_matrix = []  #matrix of concatenated strips

    #creates a matrix from vertical strips
    if strip_dir == 'v':
        #iterates through the CCDobjects (each row) and creates a strip
        for index, row in IR_list.iterrows():
            new_strip = CenterStrip(row) #creates strip object
            new_strip.makeVerticalStrip()
            strips_matrix.append(new_strip.strip)
        strips_matrix = np.array(strips_matrix)

    #creates a matrix from horizontal strips
    if strip_dir== 'h':
        #iterates through the CCDobjects (each row) and creates a strip
        for index, row in IR_list.iterrows():
            new_strip = CenterStrip(row)  #creates strip object
            new_strip.makeHorizontalStrip()
            strips_matrix.append(new_strip.strip)
        strips_matrix = np.array(strips_matrix)
    return np.transpose(strips_matrix)
    
def getLatitudes(df, channel_type):
    IR_list = df[df['channel'] == channel_type]
    listoflatitudes= []
    listofexpdates = []
    for index, row in IR_list.iterrows():
        listoflatitudes.append(TPpos(row)[0])
        listofexpdates.append(row['EXPDate'])
    return  listoflatitudes, listofexpdates


# %% Read data
start_time = DT.datetime(2023,2,18,19,30,0)
stop_time = DT.datetime(2023,2,18,23,30,0)
df = read_MATS_data(start_time,stop_time,version=0.5,level='1a')
#Can add filter in read Mats data filter = {"TPlat":[-20,20]}

# %%  Settings to run
channels = ['IR1']
strip_dir = 'v'

# %% Makes a plot of the matrix made with makeStripMatrix()
def plotKeogram(df, channels, strip_dir):

    #gives two plots, one with the matrix and one with the latitudes
    if len(channels)==1 :
        fig,axs = plt.subplots(nrows=2, ncols=1)
        latitudes, dates = getLatitudes(df,channels[0])
        matrix = makeStripMatrix(df,channels[0],strip_dir)
        axs[0].pcolormesh(dates,range(matrix.shape[0]),matrix)
        axs[0].set_title(f"Channel {channels[0]}")
        axs[1].set_xlabel('Time')
        axs[1].set_ylabel('Latitude')
        axs[1].plot(dates,latitudes,'.')
        axs[1].set_xlim(dates[0],dates[-1])
        axs[1].grid(linestyle='-')
    else:
        fig,axs = plt.subplots(nrows=len(channels)+1, ncols=1)
        latitudes, dates = getLatitudes(df,channels[0])
        axs[len(channels)].set_xlabel('Time')
        axs[len(channels)].set_ylabel('Latitude')
        axs[len(channels)].plot(dates,latitudes,'.')
        axs[len(channels)].set_xlim(dates[0],dates[-1])
        for i in range(len(channels)):
            latitudes, dates = getLatitudes(df,channels[i])
            matrix = makeStripMatrix(df,channels[i],strip_dir)  
            axs[i].pcolormesh(dates,range(matrix.shape[0]),matrix)
            axs[i].set_title(f"Channel {channels[i]}")
            plt.gcf().autofmt_xdate()
            plt.tight_layout()
       
    plt.tight_layout()
    plt.gcf().autofmt_xdate()
    plt.show()

# %%
start_time = DT.datetime(2023,2,22,19,30,0)
stop_time = DT.datetime(2023,2,24,19,30,0)
tdelta = stop_time-start_time #number of days
channel = 'IR2'
strip_dir = 'v'

df = read_MATS_data(start_time,stop_time,version=0.5,level='1a')

# %% Creates a figure with keograms for each day
def days_Keogram(df, channel, strip_dir):
    fig, axs = plt.subplots(nrows=tdelta.days, ncols=1)
    for day in range(1,tdelta.days+1): 
        print(day)
        mask_day = (df['EXPDate'] >= pd.to_datetime(start_time + timedelta(days=day-1),utc=True)) & (df['EXPDate'] <= pd.to_datetime(start_time + timedelta(days=day),utc=True))
        CCD_day = df.loc[mask_day]
        mask_night = (df['EXPDate']>= pd.to_datetime(timedelta(hours=19))) & (df['EXPDate']<= pd.to_datetime(timedelta(hours=23)))
        CCD_nighttime = CCD_day.loc[mask_night]
        print(len(CCD_day))
        latitudes, dates = getLatitudes(CCD_nighttime,channel)
        matrix = makeStripMatrix(CCD_nighttime,channel,strip_dir)  
        axs[day].pcolormesh(latitudes,range(matrix.shape[0]),matrix)
        axs[day].set_title(f"Channel {channel}")
    plt.gcf().autofmt_xdate()
    plt.tight_layout()
    plt.show()
    return fig 


# %%
