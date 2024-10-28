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
def stack(image):
    stacked_im=np.stack(image.values)
    return stacked_im


def generate_movies(df, data_folder, outputfolder,movienameprefix='movie', field_of_choise='ImageCalibrated', clim=None):
    """
    Generate a movie from a dataframe with a certain channel

    Parameters
    ----------
    df : pandas dataframe
        Dataframe with the data
    channel : str
        Channel to be used
    data_folder : str
        Folder where the images are stored
    outputfolder : str
        Folder where the movie is to be stored
    moviename : str, optional
        Name of the movie. The default is 'movie.gif'.
    field_of_choise : str, optional
        Field to be used. The default is 'ImageCalibrated'.
    clim : list, optional
        custom color limits for the plot

    """
    from mats_utils.plotting.plotCCD import orbit_plot
    from mats_utils.plotting.animate import generate_gif
    import os
    

    imagedir=data_folder+'/movies/'
    # Create directory if it does not exist
    if not os.path.exists(imagedir):
        os.makedirs(imagedir)

    

    orbit_plot(df,imagedir,nbins=7,cmap='magma', plothistogram=False, field_of_choise=field_of_choise, useplotCCDimage=True, ranges=clim)
    
    for CCDno in range(0, 8):
        CCDs = df[df['CCDSEL'] == CCDno]
        if len(CCDs) > 0:
            imagedir_CCDno=f"{imagedir}CCDSEL{str(CCDno)}"   
    
            generate_gif(imagedir_CCDno, outputfolder+'/'+movienameprefix+'_CCDSEL'+str(CCDno)+'.gif')

    print('Gifs are to be found in ',outputfolder+'/'+movienameprefix+'_CCDSEL*.gif')
    

    return



def calculate_wavelike_score(image, plot=False):

    # Apply Fourier Transform
    f = np.fft.fft2(image)
    fshift = np.fft.fftshift(f)
    magnitude_spectrum = 20 * np.log(np.abs(fshift))
    #magnitude_spectrum = f
    # Calculate the wavelike score as the sum of high-frequency components
    rows, cols = magnitude_spectrum.shape
    crow, ccol = rows // 2 , cols // 2
    #mask center of the image
    maskwidth=1
    magnitude_spectrum[crow-maskwidth:crow+maskwidth, ccol-maskwidth:ccol+maskwidth] = np.nan
    # mask the edges 
    edge=1
    magnitude_spectrum[-edge:, :] = np.nan
    magnitude_spectrum[:edge, :] = np.nan
    magnitude_spectrum[:, -edge:] = np.nan
    magnitude_spectrum[:, :edge] = np.nan

    mid_freq_magnitude = np.nanmean(magnitude_spectrum)
    #magnitude_spectrum[crow+2:crow+30, :ccol-1]
    #high_freq_magnitude = magnitude_spectrum[crow-30:crow+30, ccol-30:ccol+30]
    wavelike_score = mid_freq_magnitude

    if plot:
        
        # Display the result
        fig, ax = plt.subplots(1,2)
        sp=ax[0].imshow(magnitude_spectrum, origin="lower")#, cmap='gray')
        ax[0].set_title('Magnitude Spectrum')
        #add colorbar
        fig.colorbar(sp, ax=ax[0])
        sp=ax[1].imshow(image, origin="lower")
        ax[1].set_title('Image')
        fig.colorbar(sp, ax=ax[1])
        plt.show()


    return wavelike_score



def rolling_mean_images(df, column_name, window):
    """ Calculate the rolling mean for a column in a DataFrame.

    CALL example: df = rolling_mean_images(df, 'ImageCalibrated', window=3)

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the column to calculate the rolling mean for.
    column_name : str
        The name of the column to calculate the rolling mean for.
    window : int
        The size of the window to use for the rolling mean calculation.

    Returns
    -------
    pandas.DataFrame
        The DataFrame with the rolling mean column added.

    Raises      
    ------
    ValueError
        If the window size is not an odd number.
    """

    # Ensure the window size is an odd number
    if window % 2 == 0:
        raise ValueError("Window size must be an odd number.")
    
    # Initialize the rolling mean column with NaN values
    df['rollingmean'] = [np.nan] * len(df)
    
    # Calculate the rolling mean for each position
    half_window = window // 2
    for i in range(half_window, len(df) - half_window):
        window_values = df[column_name].iloc[i - half_window:i + half_window + 1]
        df.at[i, 'rollingmean'] = np.mean(window_values, axis=0)
    
    return df


# Function to calculate the scaled difference for 2D arrays
def scaled_diff_2d(current, previous):
    if previous is None:
        return np.nan
    mean_val = (current + previous) / 2
    return (current - previous) / mean_val

def wavelet_denoising(image, plot=False):
    # Perform wavelet transform
    coeffs = pywt.wavedec2(image, 'db1', level=2)
    print('coeffs',coeffs)
    
    # Apply thresholding
    threshold = 0.5 * np.max(coeffs[-1])
    coeffs[1:] = [list(pywt.threshold(detail, value=threshold, mode='soft') for detail in level) for level in coeffs[1:]]
    print('coeffs after thresholding',coeffs)
    # Reconstruct the image
    denoised_image = pywt.waverec2(coeffs, 'db1')

    if plot:
            # Display the result
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.imshow(image, cmap='gray')
        plt.title('Original Image')
        
        plt.subplot(1, 2, 2)
        plt.imshow(denoised_image, cmap='gray')
        plt.title('Denoised Image')
        plt.show()


    return denoised_image

def wavelet_denoising_and_thresholding(image, wavelet='db1', level=2, threshold=0.04, plot=False):
    # Perform wavelet transform
    coeffs = pywt.wavedec2(image, wavelet, level=level)
    print('coeffs',coeffs)
    
    # Apply thresholding to the detail coefficients
    threshold_value = threshold * np.max(coeffs[-1])
    #coeffs[1:] = [list(pywt.threshold(detail, value=threshold, mode='soft') for detail in level) for level in coeffs[1:]]
    coeffs[1:] = [(pywt.threshold(detail, value=threshold_value, mode='soft') for detail in level) for level in coeffs[1:]]
    

    print('coeffs after thresholding',coeffs)
    # Reconstruct the image
    denoised_image = pywt.waverec2(coeffs, wavelet)

    # Calculate wave score based on the remaining coefficients
    wave_score = np.sum([np.sum(np.abs(detail)) for detail in coeffs[1:]])

    if plot:
            # Display the result
        plt.figure(figsize=(10,10))
        plt.subplot(2, 1, 1)
        plt.imshow(image, aspect="auto",origin="lower")
        plt.title('Original Image')
        
        plt.subplot(2, 1, 2)
        plt.imshow(denoised_image, aspect="auto",origin="lower")
        plt.title('Denoised Image')
        plt.show()

    
    return denoised_image, wave_score



def analyze_wavelet_coefficients(image, wavelet='db1', level=5):
    # Perform wavelet transform
    coeffs = pywt.wavedec2(image, wavelet, level=level)
    
    # Flatten the list of coefficients
    flattened_coeffs = [item for sublist in coeffs for item in (sublist if isinstance(sublist, tuple) else [sublist])]
    
    # Calculate the energy of all coefficients
    energies = [np.sum(np.square(coeff)) for coeff in flattened_coeffs]
    
    # Calculate total energy
    total_energy = sum(energies)
    
    # Calculate the proportion of energy in each coefficient
    energy_ratios = [energy / total_energy for energy in energies]
    
    return energy_ratios

def calculate_energy_score(image, level=4, printout=False):

    energy_ratios = analyze_wavelet_coefficients(image, level=level)
    if printout:
        print('Energy ratios:',energy_ratios)
    return sum(energy_ratios[1:7])