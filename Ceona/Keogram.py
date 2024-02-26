#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

"""Class Centerstrip function library with keogram matrix funtions"""
class CenterStrip:
    def __init__(self, CCDobject):
        self.image = CCDobject['ImageCalibrated']   #for L1b otherwise 'IMAGE'
        self.strip = []
        self.latitude = CCDobject.TPlat  
        self.time =  pd.to_datetime(CCDobject['EXPDate'])
        self.maxrow = 0 #row of max intensity point
        self.maxalt = 0 #altitude of max intensity point
        self.maxlat = 0 #geodetic latitude of max intensity point
        self.maxlon = 0 #geodetic longitude of max intensity point
        self.maxI = 0 #full image integrated intensities
        #magnetic aacgm coordinates of the max intensity point
        self.MagLT = 0
        self.Maglat = 0
        self.Maglon = 0

    def makeVerticalStrip(self):
        "Makes a vertical strip object from the image"
        #finds the center pixel
        center = math.ceil(len(self.image[0])/2)
        self.strip = self.image[:,int(center)]
        return  self.strip
    
    def makeHorizontalStrip(self):
        center = math.ceil(len(self.image)/2)
        self.strip = self.image[int(center),:]
        return  np.transpose(self.strip)

def makeStripMatrix(IR_list, strip_dir ='v'):
    "Creates the keogram matrix and corresponding list of strips, from specific channel and list of CCD-items"
    strips_matrix = []  #matrix of concatenated strips
    strips_list = []
    #creates a matrix from vertical strips
    if strip_dir == 'v':
        #iterates through the CCDobjects (each panda row) and creates a strip
        for index, row in IR_list.iterrows():
            new_strip = CenterStrip(row) #creates strip object
            new_strip.makeVerticalStrip()
            strips_matrix.append(new_strip.strip)
            strips_list.append(new_strip)
        strips_matrix = np.array(strips_matrix)

    #creates a matrix from horizontal strips
    if strip_dir== 'h':
        #iterates through the CCDobjects (each panda row) and creates a strip
        for index, row in IR_list.iterrows():
            new_strip = CenterStrip(row)  #creates strip object
            new_strip.makeHorizontalStrip()
            strips_matrix.append(new_strip.strip)
            strips_list.append(new_strip)
        strips_matrix = np.array(strips_matrix)
    return np.transpose(strips_matrix) , strips_list

def getTPLatitudes(objects):
    "Return list of latitudes from list of CCD-items"
    TPlat_list= []
    for n, CCD in objects.iterrows():
        TPlat_list.append(CCD.TPlat)
    return TPlat_list

def getSatDates(objects):
    "Return list of dates from list of CCD-items"
    listofdates = []
    for n, row in objects.iterrows():
        listofdates.append(row.EXPDate)
    return listofdates

def get_stripRow(strips):
    "Get all row values from given list of strip objects (Centerstrip class)"
    allrows = strips['row']
    return allrows
# %%
