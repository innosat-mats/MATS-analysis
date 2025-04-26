import os
import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import argparse
from matplotlib import pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, satpos
from mats_l1_processing.pointing import pix_deg
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
import pdb
import netCDF4 as nc
import glob

'''
- script to plot the MATS profiles
'''

# NLC profiles for MATS
def prepare_profile(ch):
    # This function averages some columns and calculates tangent heights
    
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))
    # multiply with factor to get right values (depends on version)
    profile = np.array(image[:, col-2:col+2].mean(axis=1) * 1e12)

    common_heights = np.arange(60000, 110250, 250)
    profile = np.interp(common_heights, heights, profile)
    
    return common_heights, profile

# Define RGB colors for plotting
rgb_color_green = (0/255, 150/255, 130/255)
rgb_color_blue = (70/255, 100/255, 170/255)
rgb_color_black = (64/255, 64/255, 64/255)

def process_and_plot(file_path, output_folder):
    with open(file_path, 'rb') as file:
        df = pickle.load(file)

    for i in range(3): # len(df)
        heights, profile = prepare_profile(df.iloc[i])
        plt.plot(profile,heights, color=rgb_color_black, marker='', linestyle='-', linewidth=2, label=f'MATS')
    
        plt.tick_params(axis='both', which='major', length=8, direction='in')
        plt.tick_params(axis='both', which='minor', length=4, direction='in')
        plt.minorticks_on()
        plt.xlabel('Limb radiance (UV2)\n[ph·m$^{-2}$·s$^{-1}$·sr$^{-1}$·nm$^{-1}$]')
        plt.ylabel('Tangent altitude (km)')
        #plt.xlim(-0.25e15, 2e15)
        #plt.ylim(60, 92)
        plt.legend()
        plt.tight_layout()

        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_path = os.path.join(output_folder, f'{base_name}_profile_{i+1}.png')
        plt.savefig(output_path)
        print(f"Saved plot for profile {i+1} from {file_path} to {output_path}")
        plt.clf()  
        
def process_multiple_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    
    pkl_files = glob.glob(os.path.join(input_folder, '*.pkl'))

    for pkl_file in pkl_files:
        process_and_plot(pkl_file, output_folder)

input_folder = '../output/images/all/'
output_folder = '../output/profils/mats/'

process_multiple_files(input_folder, output_folder)
