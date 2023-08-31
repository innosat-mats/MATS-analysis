#hot pixels removal test

#%%
from mats_utils.rawdata.read_data import read_MATS_data
from database_generation.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
from mats_utils.rawdata.calibration import  calibrate_dataframe
import datetime as DT
import copy
import matplotlib.pyplot as plt
import numpy as np
from mats_utils.statistiscs.images_functions import create_imagecube
import pandas as pd
from mats_utils.geolocation.altitude_correction import rows_to_altitudes

def find_hot_pixels(image):
    # Create an array of zeros with the same shape as the image
    hot_pixels = np.zeros_like(image, dtype=bool)
    hotdiffmat=np.zeros_like(image, dtype=float)

    # Get the dimensions of the image
    height, width = image.shape

    # Iterate over each pixel in the image
    for i in range(1, height - 1):
        for j in range(1, width - 1):
            # Get the neighboring pixels for the current pixel
            neighbors = image[i-1:i+2, j-1:j+2]
            
            # Exclude the current pixel from the neighbors
            neighbors = np.delete(neighbors, 4)  # Assuming a 3x3 neighborhood
            
            # Check if the current pixel is significantly higher than its neighbors
            if image[i, j] > np.max(neighbors) and image[i, j] > np.mean(neighbors)+200:
                hot_pixels[i, j] = True
            
            hotdiffmat[i,j]=image[i,j]-neighbors.mean()

    return hot_pixels, hotdiffmat


def check_hot_pix(significantly_high_pixels, image):

    indicies= np.where(significantly_high_pixels)
    pixelvalues=[]
    neighborsvalues=[]
    
    for [i,j] in zip(indicies[0],indicies[1]):


        # Get the neighboring pixels for the current pixel
        neighbors = image[i-1:i+2, j-1:j+2]
            
        # Exclude the current pixel from the neighbors
        neighbors = np.delete(neighbors, 4)  # Assuming a 3x3 neighborhood
            

        pixelvalues.append(image[i,j])
        neighborsvalues.append(np.mean(neighbors))

    pixelvaluearray=np.array(pixelvalues)
    neighborsvaluearray=np.array(neighborsvalues)

    return  indicies, pixelvaluearray, neighborsvaluearray

def compensate_hot_pix(CCDitem, hotdiffmat):
    
    image_compensated=CCDitem['IMAGE']-hotdiffmat
    return image_compensated


#%%

starttime=DT.datetime(2023,4,10,0,0,0)
endtime=DT.datetime(2023,4,10,0,15,0)


filter_IR1={'TPsza':[97,150],'CCDSEL': [1,1], 'NRBIN': [2, 2],'NCBINCCDColumns': [40, 40],'NCOL':[43,43], 'NROW':[187,187], 'NRSKIP':[109,109]} 
ir1l1a=read_MATS_data(starttime, endtime,filter_IR1, level='1a',version='0.5')

print(len(ir1l1a))

#Calibrate
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)


ir1l1b=calibrate_dataframe(ir1l1a, instrument)


#%%
ccddf=copy.deepcopy(ir1l1a)
channel=ccddf['channel'].unique()[0]
print(channel)
#%%
# Find hot pixels using the mean image
imagecube= create_imagecube(ccddf, image_specification='IMAGE')
image=imagecube.mean(0)
hotpix, hotdiffmat=find_hot_pixels(image)

fig,ax=plt.subplots(2,1, figsize=(10,10))
sp=plot_CCDimage(image, fig, ax[0], title=channel)
sp=plot_CCDimage(hotpix, fig, ax[1], title=channel)


# %%
# Now investibgate relation between hot pixels and their neighbors

#ccddf.loc[:,'ImageCalibrated']= ccddf.apply(check_hot_pix, args=(fixaltvec, 'ImageCalibrated'), axis=1)
# check_hot_pix

fig, ax=plt.subplots(1,1, figsize=(10,10))
fig2, ax2=plt.subplots(10,1, figsize=(3,30))

# Initialize an empty list to store rows 
rows = []

fullpixelvalues=[]
fullneighborsvalues=[]

nrhotpix=sum(sum(hotpix))
ivalues= np.empty((len(ccddf), nrhotpix))
jvalues= np.empty((len(ccddf), nrhotpix))
pixvalues= np.empty((len(ccddf), nrhotpix))
neighborvalues= np.empty((len(ccddf), nrhotpix))



for irow, row in ccddf[:10].iterrows():
    # Access row data using dot notation or index
    indicies, pixelvaluesarray, neighborsvaluearray=check_hot_pix(hotpix, row.IMAGE)
    
    fullpixelvalues.append(pixelvaluesarray)
    fullneighborsvalues.append(neighborsvaluearray)   

    for index in range(len(indicies[0])):
        ivalues[irow,index]=indicies[0][index]
        jvalues[irow,index]=indicies[1][index]
        pixvalues[irow,index]=pixelvaluesarray[index]   
        neighborvalues[irow,index]=neighborsvaluearray[index]   



for index in range(10):

    ax2[index].scatter(neighborvalues[:,index], pixvalues[:,index])
    ax2[index].set_xlabel('neighborvalue')
    ax2[index].set_ylabel('pixvalue')
    ax2[index].set_title('Pixvalue vs. Neighborvalue')


fullpixelvaluesarray=np.array(fullpixelvalues).flatten()
fullneighborsvaluesarray=np.array(fullneighborsvalues).flatten()
sp=ax.scatter(fullneighborsvaluesarray, fullpixelvaluesarray)
ax.set_ylabel('Pixel value')
ax.set_xlabel('Mean of neighbors')

poly=np.polyfit(fullneighborsvaluesarray,fullpixelvaluesarray, 1)
ax.plot(fullneighborsvaluesarray, np.polyval(poly, fullneighborsvaluesarray), 'r-')


# %%
# Now compensate for hot pixels
#ccddf.loc[:,'image_hotpix']= ccddf.apply(compensate_hot_pix, args=(hotdiffmat))


ccddf['ImageHotPixComp'] = ccddf.apply(lambda row: compensate_hot_pix(row, hotdiffmat), axis=1)
#%%
for row in ccddf.itertuples(index=False):
    fig, ax=plt.subplots(2,1)
    sp=plot_CCDimage(row.ImageHotPixComp, fig, ax[0], title=channel)
    sp=plot_CCDimage(row.IMAGE, fig, ax[1], title=channel)

# fixaltvec=np.arange(61000, 107000, 1000)
# ccddf.loc[:,'image_fixalt']= ccddf.apply(rows_to_altitudes, args=(fixaltvec, 'ImageCalibrated'), axis=1)

# %%
