#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd 
import sys
import numpy as np
import math
import matplotlib.pyplot as plt 
from matplotlib import ticker
from mats_utils.geolocation.coordinates import TPpos, satpos
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot, all_channels_plot

#%%
class CenterStrip:
    def __init__(self, CCDobject, width):
        self.image = CCDobject['IMAGE']
        self.width = width
        self.strip = self.makestrip()
        self.latitude = TPpos(CCDobject)[0] 
        self.time =  pd.to_datetime(CCDobject['EXPDate'])

    "Makes a strip object from the image"
    def makestrip(self):
        "Number of columns in the image"
        #center = math.ceil(len(self.image[0])/2)
        #self.strip = self.image[:,int(center-self.width/2):int(center+self.width/2)]

        #changed to horizontal strip
        center = math.ceil(len(self.image)/2)
        self.strip = self.image[int(center-self.width/2):int(center+self.width/2),:]
       
        return  np.transpose(self.strip)

#%%  Make a keogram of a specific channel and CCD-objects
def makeStripMatrix(df, channel_type, width):

    IR_list = df[df['channel'] == channel_type]
    strips_matrix = []  #matrix of concatenated strips
    listofstrips = []   #list of strip objects

    #iterates through the CCDobjects (each row) and creates a strip
    for index, row in IR_list.iterrows():
        new_strip = CenterStrip(row, width)  #creates strip object
        listofstrips.append(new_strip) 
        for i in range(width):
            strips_matrix.append(new_strip.strip[:,i])
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


# %%
start_time = DT.datetime(2023,2,8,19,15,0)
stop_time = DT.datetime(2023,2,8,20,30,0)
df = read_MATS_data(start_time,stop_time)

#Can add filter in read Mats data filter = {"TPlat":[-20,20]}

# %%
width = 1
channels = ['IR1']

# %%
def plotKeogram(df, channels, width):
    fig,axs = plt.subplots(nrows=2, ncols=1)
    #axs = axs.flatten()

    #for i,ax in enumerate(axs):
    i = 0
    matrix = makeStripMatrix(df,channels[i],width)
    latitudes, dates = getLatitudes(df,channels[i])
   
    axs[0].pcolormesh(dates,range(matrix.shape[0]),matrix)
    axs[0].set_title(f"Channel {channels[i]}")
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Latitude')
    axs[1].plot(dates,latitudes,'.')
    axs[1].set_xlim(dates[0],dates[-1])

    plt.tight_layout()
    plt.gcf().autofmt_xdate()
    plt.show()

# %%
