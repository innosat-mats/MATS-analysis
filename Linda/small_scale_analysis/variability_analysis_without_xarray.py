
# %%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import argparse
from datetime import date, timedelta
from mats_utils.plotting.plotCCD import all_channels_plot, make_ths
#from mats_utils.daily_preview.temp_nadirs import NADIR_geolocation, average_stacking
import numpy as np
import pandas as pd
import multiprocessing
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import ceil
import cartopy.crs as ccrs
import sys
import os
from multiprocessing import Manager
import matplotlib
import xarray as xr
import scipy
import proplot as pplt
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
from database_generation.experimental_utils import plot_CCDimage
import pywt
import cv2
import numpy as np

from variability_analysis_functions import calculate_wavelike_score, calculate_energy_score, analyze_wavelet_coefficients, wavelet_denoising_and_thresholding
from variability_analysis_functions import generate_movies
# %%
# ASCENDING MAP PLOT /// UPDATED TO READ FROM L1b_v05


#axs.format(lonlabels='b')
#axs[0].format(latlabels='l')
#axs[3].format(latlabels='l')


#monthstrings=["022023_2","032023_1","032023_2","042023_1"]
monthstrings=["022023_1","032023_1","042023_1","022023_2","032023_2",False]
monthstrings=["022023","032023","042023"]
monthstrings=["022023"]

#monthstrings=["022023_1","032023_1",False,False,False,False]
#monthstrings=[False,"032023_1",False,False,False,False]
# dayglow
dmin,dmax = 0, 95


directory = '/Users/lindamegner/MATS/MATS-retrieval/data/variability/'
os.makedirs(directory, exist_ok=True)

#subtitles=['Feb (1st-14th)']#, 'Feb (15st-28th)','March (1st-14th)', 'March (15st-31st)','April (1st-14th)', 'April (15st-31st)']
#for monthstring in monthstrings:
monthstring="022023"
week=2
#    for offset in [1]:
#        for j in [1]:
#            week = j + offset
#            print(f'Processing {monthstring}_{week}')
datasource='pickled'
if datasource=='bjorn':
    df=pd.read_pickle(f'/Users/lindamegner/MATS/MATS-retrieval/data/BjornsData/{monthstring}_{week}.pk1')
    df.ImageCalibrated = df.ImageCalibrated.div(df.TEXPMS/1000)

elif datasource=='aws':
    start_time = DT.datetime(2023, 2, 4, 9, 0, 0)
    stop_time = DT.datetime(2023, 2, 5, 9, 0, 0)
    df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
    df['TPlocaltime'] = pd.to_datetime(df['TPlocaltime'],utc=True)
    #pickle the dataframe
    name='df_'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(directory+name+'.pkl')
    print('pickling:'+name)
elif datasource=='pickled':
    df=pd.read_pickle(f'{directory}df_20230204_090000_20230204_103000.pkl') #one orbit
    #df=pd.read_pickle(f'{directory}df_20230204_090000_20230205_090000.pkl') # one day


df_saveme=df.copy()



#%%
#df=df[:100]

# convert string df['TPlocaltime'] to datetime

#df=pd.read_pickle(f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL1/alts_80to90/finished/{monthstring}_{week}.pk1')

midday=DT.time(12, 0, 0)

df = df_saveme[(df_saveme['channel'] == 'IR2')]

# ascending
ascending = True
if ascending:
    df = df[(df['TPlocaltime'].dt.time > midday)]
else:
    df = df[(df['TPlocaltime'].dt.time < midday)]

#reset index
df.reset_index(drop=True, inplace=True)
#df.ImageCalibrated= df.ImageCalibrated.apply(np.mean, axis=1)

add_field_with_subtracted_rolling_mean2(df, 'ImageCalibrated', 'im_detrended', window_before=3, window_after=3, skipbefore=0, skipafter=0)

df['im_diff'] = df.ImageCalibrated.diff()

# Initialize 'im_detrended' with NaN arrays of the same shape as 'ImageCalibrated'
#df['im_detrended'] = df['ImageCalibrated'].apply(lambda x: np.full_like(x, np.nan))
#for i in range(1, len(df)):
#    df.at[i, 'im_detrended'] = scaled_diff_2d(df.at[i, 'ImageCalibrated'], df.at[i-1, 'ImageCalibrated'])
    
#cut away when the time difference is more than 20 seconds
df = df[df['TPlocaltime'].diff().dt.total_seconds() < 10]
# drop first row
df = df.drop(df.index[0])
#df=df[1:-1]
#df['im_detrended'].iloc[0] = df.im_detrended.iloc[1] #hack to avaid nans



#plot_CCDimage(df['ImageCalibrated_minus_rolling_mean'].iloc[0])
#df['im_detrended_filtered'] = df.im_detrended.apply(lambda x: scipy.ndimage.median_filter(x, size=(3, 3)))
df['im_detrended_filtered'] = df.im_detrended.apply(lambda x: scipy.ndimage.median_filter(x, size=(3, 3)))
df['im_diff_filtered'] = df.im_diff.apply(lambda x: scipy.ndimage.median_filter(x, size=(3, 3)))
#cut away when the mean is an outlier, ie when the mean is more than 3 sigma from the mean
#df['im_detrended_filtered_mean'] = df.im_detrended_filtered.apply(np.mean)
#df = df[np.abs(df['im_detrended_filtered_mean'] - df['im_detrended_filtered_mean'].mean()) < 3 * df['im_detrended_filtered_mean'].std()]

irow=60
# icol=20

# # divide with the absolute mean of the image

# df['im_std_vertical']= df.im_detrended_filtered.apply(lambda x: np.std(x[irow,:]))
# df['im_std_horizontal']= df.im_detrended_filtered.apply(lambda x: np.std(x[:,icol]))
df['ImageCalibrated_hormean']= df.ImageCalibrated.apply(lambda x: np.mean(x[irow,:]))
# df['ImageCalibrated_vermean']= df.ImageCalibrated.apply(lambda x: np.mean(x[:,icol]))
df['fourier_score'] = df['im_detrended_filtered'].apply(lambda x: calculate_wavelike_score(x[:90,:], plot=False))
df['variance_score'] = df['im_detrended_filtered'].apply(lambda x: np.var(x[:90,:], axis=1).mean())
df['energy_score'] = df['im_detrended'].apply(lambda x: calculate_energy_score(x[:90,:], level=4))
df['fourier_score_diff'] = df['im_diff_filtered'].apply(lambda x: calculate_wavelike_score(x[:90,:], plot=False))
df['variance_score_diff'] = df['im_diff_filtered'].apply(lambda x: np.var(x[:90,:], axis=1).mean())
df['energy_score_diff'] = df['im_diff'].apply(lambda x: calculate_energy_score(x[:90,:], level=4))
df['aurora_indicator']=df.ImageCalibrated.apply(lambda x: np.mean(x[180,:]))





#%%


fig = pplt.figure(figwidth='14cm')
axs = fig.subplots(ncols=1, nrows=5, proj='cyl')
fig.format(suptitle='Dayglow observations from IR1 (~17:30 LT)',abc='a.',abcborder=False)
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

m=axs[0].scatter(df.TPlon,df.TPlat,c=df.ImageCalibrated_hormean,s=0.5,cmap='Thermal')
axs[0].format(title='Airglow at 86 km')
#add colorbar to the subplot
axs[0].colorbar(m,label='Airglow intensity')

m=axs[1].scatter(df.TPlon,df.TPlat,c=np.log(df.variance_score),s=0.5,cmap='Thermal')
#m=axs[1].scatter(df.TPlon,df.TPlat,c=df.variance_score,s=0.3,cmap='Thermal', vmin=0, vmax=50)
axs[1].format(title='Variability score')
axs[1].colorbar(m,label='Variance score')

m=axs[2].scatter(df.TPlon,df.TPlat,c=df.fourier_score,s=0.5,cmap='Thermal')
axs[2].format(title='Fourier score')
axs[2].colorbar(m,label='Fourier score')

m=axs[3].scatter(df.TPlon,df.TPlat,c=np.log(df.energy_score),s=0.5,cmap='Thermal', vmin=-3, vmax=-0.1)
axs[3].format(title='Energy score')
axs[3].colorbar(m,label='Energy score')

m=axs[4].scatter(df.TPlon,df.TPlat,c=np.log(df.energy_score_diff),s=0.5,cmap='Thermal', vmin=-3, vmax=-0.1)
axs[4].format(title='Energy score diff')
axs[4].colorbar(m,label='Energy score diff')



#BELOW HERE ONLY FOR TESTING


 #%%

#field_of_choise='im_detrended'
#generate_movies(df[100:110], '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/temp', '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output', movienameprefix='movie_'+field_of_choise,field_of_choise=field_of_choise, clim=[-10,10])
field_of_choise='im_detrended_filtered'

# generate a random name for the folder
import random
import string
imagefolder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/temp/'+''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
generate_movies(df[10:-10], imagefolder, '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output', movienameprefix='movie_'+field_of_choise,field_of_choise=field_of_choise, clim=[-15,15], cmap='YlGnBu_r')
#%%
# om programmet crashar och du vill generera film af de bilder du redan har:
from mats_utils.plotting.animate import generate_gif
fullimagefolder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/temp/336RL/movies/CCDSEL4'
generate_gif(fullimagefolder, '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/movie1.gif')

#%%

#sort the dataframe by wavelike score
df = df.sort_values(by='fourier_score_diff', ascending=False)
df.reset_index(drop=True, inplace=True)
    
clim = [-20, 20]
for index, dfrow in df[:10].iterrows():
    plot_CCDimage(dfrow['im_detrended_filtered'], title=f'variance: {dfrow["variance_score"]:.2g}, fourier: {dfrow["fourier_score"]:.2g}, energy: {dfrow["energy_score"]:.2g}', clim=clim, cmap='viridis')  
    wenergy=calculate_energy_score(dfrow['im_detrended'], level=5, printout=True)
    #plt.plot(wenergy, label=str(index))
    #plt.legend()

#%%
dfshort=df[:10]
dfshort['im_detrended_filtered'].apply(lambda x: calculate_wavelike_score(x[:100,:], plot=True))


#%%

#df_with_ssa=pd.read_pickle('../output/temp/df_with_ssa.pkl')
#df_no_ssa=pd.read_pickle('../output/temp/df_no_ssa.pkl')  
#df_no_ssa=df

for index, dfrow in df_no_ssa[:5].iterrows():
    #fig, ax = plt.subplots(2,1)
    image=dfrow['im_detrended'][:100,:]
    #plot_CCDimage(image, fig=fig, axis=ax[0],title='Without SSA')
    wenergy=analyze_wavelet_coefficients(image, level=5)
    print('wenergy without ssa:',wenergy)
    plt.plot(wenergy, label='without ssa', color='blue')


for index, dfrow in df_with_ssa[:5].iterrows():
    #fig, ax = plt.subplots(2,1)
    image=dfrow['im_detrended'][:100,:]
    #plot_CCDimage(image, fig=fig, axis=ax[0], title='With SSA')
    wenergy=analyze_wavelet_coefficients(image, level=5)
    print('wenergy with ssa:',wenergy)
    plt.plot(wenergy, label='with ssa', color='red')

plt.legend()
plt.show()
# Example usage
image=df_no_ssa['im_detrended'].iloc[0][:100,:]

denoisedimage, wavescore=wavelet_denoising_and_thresholding(image, level=5, threshold=0.6, plot=False)
#%%
#
# 
for index, dfrow in df[:5].iterrows():
    #fig, ax = plt.subplots(2,1)
    image=dfrow['im_detrended_filtered'][:100,:]
    #plot_CCDimage(image, fig=fig, axis=ax[0], title='With SSA')
    wenergy=analyze_wavelet_coefficients(image, level=4)
    print('wenergy high rated:',wenergy)
    plt.plot(wenergy, label='high rated', color='red')
for index, dfrow in df[-5:].iterrows():
    #fig, ax = plt.subplots(2,1)
    image=dfrow['im_detrended_filtered'][:100,:]
    #plot_CCDimage(image, fig=fig, axis=ax[0], title='With SSA')
    wenergy=analyze_wavelet_coefficients(image, level=4)
    print('wenergy low rated:',wenergy)
    plt.plot(wenergy, label='low rated', color='blue')



plt.legend()
plt.show()
# 
# 
#  
# 
# 
#%%

# image=df_no_ssa['im_detrended'].iloc[0][:100,:]
# plot_CCDimage(image, title='original image')
# coeffs = pywt.wavedec2(image, 'db1', level=5)
    
# # Flatten the list of coefficients
# flattened_coeffs = [item for sublist in coeffs for item in (sublist if isinstance(sublist, tuple) else [sublist])]

# coeffs[0]=coeffs[0]*0

# denoised_image = pywt.waverec2(coeffs, 'db1')
# plot_CCDimage(denoised_image, title='denoised image')

    
# # Calculate the energy of all coefficients
# energies = [np.sum(np.square(coeff)) for coeff in flattened_coeffs]
    
# # Calculate total energy
# total_energy = sum(energies)
    
# # Calculate the proportion of energy in each coefficient
# energy_ratios = [energy / total_energy for energy in energies]

# print('energy ratios:',energy_ratios)
# %%

#plot df['aurora_indicator'] against df['TPlat']
fig, ax = plt.subplots()
ax.scatter(df['TPlat'],df['aurora_indicator'])
ax.set_xlabel('Latitude')
ax.set_ylabel('Aurora indicator')

image=df['im_detrended_filtered'].iloc[0]
plot_CCDimage(image, title='original image')

# %%
