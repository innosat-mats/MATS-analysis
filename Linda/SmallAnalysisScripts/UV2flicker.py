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
#%%

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output'
starttime=DT.datetime(2023,5,7,1,0,0)
endtime=DT.datetime(2023,5,7,2,0,0)


filter_UV2={'CCDSEL': [6,6]} 
df=read_MATS_data(starttime, endtime,filter_UV2, level='1b',version='0.5')
print(len(df))

# %%
for index, CCD in df[:2].iterrows():
    plot_image(CCD, outpath=data_folder)
    #CCD['ImageCalibrated']

    
#%%
df['ImageMean']=df['ImageCalibrated'].apply(lambda x: np.mean(x[130,:]))
#df['ImageMean']=df['ImageCalibrated'].mean()

plt.plot( df['EXPDate'],df['ImageMean'])
#simple_plot(df,data_folder)
# %%
from scipy.fft import fft, fftfreq
import numpy as np
# Calculate the Fourier Transform
signal = df['ImageMean'].values

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
counts, bins=np.histogram(df['ImageMean'],bins=100)
plt.stairs(counts,bins)
plt.xlabel('Counts')
plt.ylabel('Nr of occurences')


#%%
#Fourier analysis version 1
time = df.index.values
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
