# Author: Linda Megner
# Created: October 26, 2023
# Script to investigate the flickering in UV2

#%%
from mats_utils.rawdata.read_data import read_MATS_data
from database_generation.experimental_utils import plot_CCDimage
from database_generation.flatfield import read_flatfield
from mats_l1_processing.instrument import Instrument
import datetime as DT
from mats_utils.rawdata.calibration import  calibrate_dataframe
from mats_l1_processing.L1_calibration_functions import calculate_flatfield, bin_image_with_BC, meanbin_image_with_BC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from lindas_own_functions import rename_CCDitem_entries
from database_generation import flatfield as flatfield
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_l1_processing.L1_calibrate import L1_calibrate

from mats_l1_processing.L1_calibration_functions import (
    get_true_image,
    desmear_true_image,
    subtract_dark,
    flatfield_calibration,
    get_linearized_image,
    combine_flags,
    make_binary,
    flip_image,
    handle_bad_columns,
    artifact_correction,
)

#%%

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output'
starttime=DT.datetime(2023,5,7,1,0,0)
endtime=DT.datetime(2023,5,7,2,0,0)


filter_UV2={'CCDSEL': [6,6]} 
df=read_MATS_data(starttime, endtime,filter_UV2, level='1a',version='0.6')
print(len(df))

# %%
for index, CCD in df[:2].iterrows():
    plot_image(CCD, outpath=data_folder)
    #CCD['ImageCalibrated']


CCDitems=df[:min(500,len(df))]

#%%
# Calibrate the data
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)
def calibrate(CCDitem, instrument):

    (

    image_lsb, 
    image_bias_sub, 
    image_desmeared, 
    image_dark_sub, 
    image_calib_nonflipped, 
    image_calib_flipped, 
    image_calibrated, 
    errors

    ) = L1_calibrate(CCDitem, instrument)

    CCDitem["image_lsb"] = image_lsb
    CCDitem["image_bias_sub"] = image_bias_sub
    CCDitem["image_desmeared"] = image_desmeared
    CCDitem["image_dark_sub"] = image_dark_sub
    CCDitem["image_calib_nonflipped"] = image_calib_nonflipped
    CCDitem["image_calib_flipped"] = image_calib_flipped
    CCDitem["image_calibrated"] =image_calibrated
    CCDitem["errors"] = errors

    return image_calibrated

def wrap_get_true_image(CCDitem):
    image_bias_sub, error=get_true_image(CCDitem)
    return image_bias_sub
def wrap_get_linearized_image(CCDitem):    
    image_linear, error = get_linearized_image(CCDitem, CCDitem["image_bias_sub"], force_table=True)
    return image_linear
def wrap_desmear_true_image(CCDitem):
    image_desmeared, error=desmear_true_image(CCDitem, CCDitem["image_linear"])
    return image_desmeared
def wrap_subtract_dark(CCDitem):
    image_dark_sub, error=subtract_dark(CCDitem, CCDitem["image_desmeared"])
    return image_dark_sub
def wrap_flatfield_calibration(CCDitem):
    image_calib_nonflipped, error=flatfield_calibration(CCDitem, CCDitem["image_dark_sub"])
    return image_calib_nonflipped
def wrap_flip_image(CCDitem):
    image_calib_flipped=flip_image(CCDitem, CCDitem["image_calib_nonflipped"])
    return image_calib_flipped
def wrap_artifact_correction(CCDitem):
    image_calibrated, error=artifact_correction(CCDitem, CCDitem["image_calib_flipped"])
    return image_calibrated





  
#%%
CCDitems=rename_CCDitem_entries(CCDitems)
CCDitems['ImageCalibrated']=CCDitems.apply(lambda CCDitem: calibrate(CCDitem,instrument ),axis=1)

CCDitems["CCDunit"] =CCDitems.apply(lambda CCDitem: instrument.get_CCD(CCDitem["channel"]),axis=1)  

#Calibration step by step
CCDitems["image_bias_sub"]=CCDitems.apply(lambda CCDitem: wrap_get_true_image(CCDitem) ,axis=1)

#%%
CCDitems["image_linear"]=CCDitems.apply(lambda CCDitem : wrap_get_linearized_image(CCDitem) ,axis=1)
CCDitems["image_desmeared"]=CCDitems.apply(lambda CCDitem: wrap_desmear_true_image(CCDitem) ,axis=1)   
CCDitems["image_dark_sub"]=CCDitems.apply(lambda CCDitem: wrap_subtract_dark(CCDitem) ,axis=1)
CCDitems["image_calib_nonflipped"]=CCDitems.apply(lambda CCDitem: wrap_flatfield_calibration(CCDitem) ,axis=1)
CCDitems["image_calib_flipped"]=CCDitems.apply(lambda CCDitem: wrap_flip_image(CCDitem) ,axis=1)      
CCDitems["image_calibrated"]=CCDitems.apply(lambda CCDitem: wrap_artifact_correction(CCDitem) ,axis=1)   



#%%
ystart=130
CCDitems['IMAGEMean']=CCDitems['IMAGE'].apply(lambda x: np.mean(x[ystart,:]))
#CCDitems['ImageCalibratedMean']=CCDitems['ImageCalibrated'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_bias_subMean']=CCDitems['image_bias_sub'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_linearMean']=CCDitems['image_linear'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_desmearedMean']=CCDitems['image_desmeared'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_dark_subMean']=CCDitems['image_dark_sub'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_calib_nonflippedMean']=CCDitems['image_calib_nonflipped'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_calib_flippedMean']=CCDitems['image_calib_flipped'].apply(lambda x: np.mean(x[ystart,:]))
CCDitems['image_calibratedMean']=CCDitems['image_calibrated'].apply(lambda x: np.mean(x[ystart,:]))










fig2, ax2=plt.subplots(8,1, figsize=(8, 20))
ax2[0].plot(CCDitems['EXPDate'],CCDitems['IMAGEMean'])
ax2[0].set_title(CCDitems.iloc[0].channel+' Mean value of image')
ax2[1].plot(CCDitems['EXPDate'],CCDitems['image_bias_subMean'])
ax2[1].set_title('Mean value of image_bias_sub')
ax2[2].plot(CCDitems['EXPDate'],CCDitems['image_linearMean'])
ax2[2].set_title('Mean value of image_linear')
ax2[3].plot(CCDitems['EXPDate'],CCDitems['image_desmearedMean'])
ax2[3].set_title('Mean value of image_desmeared')
ax2[4].plot(CCDitems['EXPDate'],CCDitems['image_dark_subMean'])
ax2[4].set_title('Mean value of image_dark_sub')
ax2[5].plot(CCDitems['EXPDate'],CCDitems['image_calib_nonflippedMean'])
ax2[5].set_title('Mean value of flat field compensated image')
ax2[6].plot(CCDitems['EXPDate'],CCDitems['image_calib_flippedMean'])
ax2[6].set_title('Mean value of image_calib_flipped')
ax2[7].plot(CCDitems['EXPDate'],CCDitems['image_calibratedMean'])
ax2[7].set_title('Mean value of image_calibrated')
fig2.autofmt_xdate()



i=7
fig, ax=plt.subplots(8,1, figsize=(8, 20))
plot_CCDimage(CCDitems.iloc[i]['IMAGE'], fig, ax[0], title='Original image')
plot_CCDimage(CCDitems.iloc[i]['image_bias_sub'], fig, ax[1], title='Bias subtracted image')    
plot_CCDimage(CCDitems.iloc[i]['image_linear'], fig, ax[2], title='Linearized image')
plot_CCDimage(CCDitems.iloc[i]['image_desmeared'], fig, ax[3], title='Desmeared image')
plot_CCDimage(CCDitems.iloc[i]['image_dark_sub'], fig, ax[4], title='Dark current subtracted image')    
plot_CCDimage(CCDitems.iloc[i]['image_calib_nonflipped'], fig, ax[5], title='Flat field compensated image') 
plot_CCDimage(CCDitems.iloc[i]['image_calib_flipped'], fig, ax[6], title='Flipped image')   
plot_CCDimage(CCDitems.iloc[i]['image_calibrated'], fig, ax[7], title='Artifact corrected image')   




#simple_plot(CCDitems,data_folder)

# %%
#FOURIER ANALYSIS

from scipy.fft import fft, fftfreq
import numpy as np
# Calculate the Fourier Transform
signal = CCDitems['ImageMean'].values

yf = fft(signal)
N=len(signal)
T=6

xf = fftfreq(N, T)[:N//2]

import matplotlib.pyplot as plt

plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))

plt.grid()
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Fourier Spectrum of ImageMean')

plt.figure(figsize=(12, 6))
counts, bins=np.histogram(CCDitems['ImageMean'],bins=100)
plt.stairs(counts,bins)
plt.xlabel('Counts')
plt.ylabel('Nr of occurences')


#%%
#Fourier analysis version 1
time = CCDitems.index.values
sample_rate = 1/6 # Calculate the sample rate
fourier_transform = np.fft.fft(signal)
frequencies = np.fft.fftfreq(len(signal), 1 / sample_rate)
plt.figure(figsize=(12, 6))
plt.plot(frequencies, np.abs(fourier_transform))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Fourier Spectrum of ImageMean')
plt.grid()
plt.show()# %%


# %%
#Fourier analysis version 2

from scipy.fft import fft, fftfreq
import numpy as np

# Number of sample points

N = 600

# sample spacing

T = 1.0 / 800.0

x = np.linspace(0.0, N*T, N, endpoint=False)

#y = np.random.normal(size=N)
y=np.sin(50000000.0 * 2.0*np.pi*x)

yf = fft(y)

xf = fftfreq(N, T)[:N//2]

import matplotlib.pyplot as plt

plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))

plt.grid()

plt.show()



plt.plot(y)

plt.grid()

plt.show()
# %%
