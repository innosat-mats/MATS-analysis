
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
from mats_utils.rawdata.read_data import load_multi_parquet
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
#dmin,dmax = 0, 95



directory = '/home/linda/data/variability/'
os.makedirs(directory, exist_ok=True)



start = DT.datetime(2023, 2, 15, 0, 0, 0)
stop= DT.datetime(2023, 4, 17, 0, 0, 0)

tstarttime = start
while tstarttime < stop:
    tstoptime = tstarttime + DT.timedelta(days=14)
    try:
        df= read_MATS_data(tstarttime, tstoptime, level='1a', version='0.7', pritfilesys=False)
        #df=load_multi_parquet('/home/linda/MATS_data/CROPD_v0.6',start, stop)
    except:
        print('found no data for intervalstarttime', tstarttime)

    midday=DT.time(12, 0, 0)

    #select channel
    df = df[(df['channel'] == 'IR1')]
    # ascending
    ascending = True
    if ascending:
        df = df[(df['TPlocaltime'].dt.time > midday)]
    else:
        df = df[(df['TPlocaltime'].dt.time < midday)]

    #reset index

    df.reset_index(drop=True, inplace=True)


    add_field_with_subtracted_rolling_mean2(df, 'ImageCalibrated', 'im_detrended', window_before=3, window_after=3, skipbefore=0, skipafter=0)

    df['im_diff'] = df.ImageCalibrated.diff()

    #cut away when the time difference is more than 20 seconds
    df = df[df['TPlocaltime'].diff().dt.total_seconds() < 10]
    # drop first row
    df = df.drop(df.index[0])



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





    fig = pplt.figure(figwidth='14cm')
    axs = fig.subplots(ncols=1, nrows=5, proj='cyl')
    fig.format(suptitle='Dayglow observations from IR1 (~17:30 LT)',abc='a.',abcborder=False)
    axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)

    m=axs[0].scatter(df.TPlon,df.TPlat,c=df.ImageCalibrated_hormean,s=0.3,cmap='Thermal')
    axs[0].format(title='Airglow at 86 km')
    #add colorbar to the subplot
    axs[0].colorbar(m,label='Airglow intensity')

    m=axs[1].scatter(df.TPlon,df.TPlat,c=np.log(df.variance_score),s=0.3,cmap='Thermal')
    #m=axs[1].scatter(df.TPlon,df.TPlat,c=df.variance_score,s=0.3,cmap='Thermal', vmin=0, vmax=50)
    axs[1].format(title='Variability score')
    axs[1].colorbar(m,label='Variance score')

    m=axs[2].scatter(df.TPlon,df.TPlat,c=df.fourier_score,s=0.3,cmap='Thermal')
    axs[2].format(title='Fourier score')
    axs[2].colorbar(m,label='Fourier score')

    m=axs[3].scatter(df.TPlon,df.TPlat,c=np.log(df.energy_score),s=0.3,cmap='Thermal', vmin=-3, vmax=-0.1)
    axs[3].format(title='Energy score')
    axs[3].colorbar(m,label='Energy score')

    m=axs[4].scatter(df.TPlon,df.TPlat,c=np.log(df.energy_score_diff),s=0.3,cmap='Thermal', vmin=-3, vmax=-0.1)
    axs[4].format(title='Energy score diff')
    axs[4].colorbar(m,label='Energy score diff')

    #save figure to file named after date
    if ascending:
        fig.savefig(directory + 'airglow_ascending'+tstarttime.strftime('%Y%m%d')+'.png')
    else:
        fig.savefig(directory + 'airglow_descending'+tstarttime.strftime('%Y%m%d')+'.png')




# #field_of_choise='im_detrended'
# #generate_movies(df[100:110], '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/temp', '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output', movienameprefix='movie_'+field_of_choise,field_of_choise=field_of_choise, clim=[-10,10])
# field_of_choise='im_detrended_filtered'

# # generate a random name for the folder
# import random
# import string
# imagefolder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/temp/'+''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
# generate_movies(df[100:120], imagefolder, '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output', movienameprefix='movie_'+field_of_choise,field_of_choise=field_of_choise, clim=[-15,15])



# %%
