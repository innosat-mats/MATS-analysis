#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd 
import numpy as np
import math
import matplotlib.pyplot as plt 
from mats_utils.geolocation.coordinates import TPpos, satpos
from matplotlib.backends.backend_pdf import PdfPages

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

# Make a keogram of a specific channel and CCD-objects
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

#Get a list of latitudes and dates for TP
def getTimePos(df):
    listoflatitudes= []
    listofexpdates = []
    for index, row in df.iterrows():
        listoflatitudes.append(TPpos(row)[0])
        listofexpdates.append(row['EXPDate'])
    return listoflatitudes, listofexpdates

# Makes a plot of the matrix made with makeStripMatrix()
def plotKeogram(df, channels, strip_dir):
    
    #gives at least 2 plots, one from image-matrix and one with the latitudes
    if len(channels)==1 :
        IR_list = df[df['channel'] == channels[0]]
        fig,axs = plt.subplots(nrows=2, ncols=1)
        latitudes, dates = getTimePos(IR_list)
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
        IR1_list = df[df['channel'] == channels[0]]
        # gets the latitudes from objects only for the first channel.
        latitudes, dates = getTimePos(IR1_list)
        axs[len(channels)].set_xlabel('Time')
        axs[len(channels)].set_ylabel('Latitude')
        axs[len(channels)].plot(dates,latitudes,'.')
        axs[len(channels)].set_xlim(dates[0],dates[-1])
        for i in range(len(channels)):
            IR_objects =  df[df['channel'] == channels[i]]
            latitudes, dates = getTimePos(IR_objects)
            matrix = makeStripMatrix(df,channels[i],strip_dir)  
            axs[i].pcolormesh(dates,range(matrix.shape[0]),matrix, vmax=1500)
            axs[i].set_title(f"Channel {channels[i]}")
            plt.gcf().autofmt_xdate()
            plt.tight_layout()
       
    plt.tight_layout()
    plt.gcf().autofmt_xdate()
    plt.show()

# %%  Settings to run for normal plot with various channels
start_time = DT.datetime(2023,2,19,18,30,0)
stop_time = DT.datetime(2023,2,20,18,30,0)
channels = ['IR1']  #list of what channels we want to plot
strip_dir = 'v'  #vertical 'v' or horiozontal 'h' direction of the strip
df = read_MATS_data(start_time,stop_time,version=0.5,level='1a')

# %%
