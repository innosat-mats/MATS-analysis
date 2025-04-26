import os
import pandas as pd
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pickle
from mats_utils.geolocation.coordinates import col_heights
import json
import pdb
from matplotlib.ticker import AutoMinorLocator

'''
- profiles of MATS and Odin used in thesis

'''

# Function to prepare MATS profile
def prepare_profile(ch):
    # This function averages some columns and calculates tangent heights
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL']/2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    # Multiply with factor to get right values (depends on version?) --> for version 0.6 = 1e12
    profile = np.array(image[:, col-2:col+2].mean(axis=1) * 1e12)

    common_heights = np.arange(60000, 110250, 250)
    profile = np.interp(common_heights, heights, profile)
    return common_heights, profile

# Load MATS data
csv_data = pd.read_csv('/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/osiris_project/git/Analysis/output/radiance_profiles/close_encounters/mats_IDs_February_23.csv')
# csv_data = pd.read_csv('/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/osiris_project/git/Analysis/output/radiance_profiles/close_encounters/mats_IDs_December_22.csv')

# csv_data = pd.read_csv('/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/osiris_project/git/Analysis/scripts/filtered_result_df.csv')

file_path_mats = '/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/mats_project/MATS-analysis/Aglaja/output/images/Feb_23/Feb_23_new.pkl'
# file_path_mats = '/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/mats_project/MATS-analysis/Aglaja/output/images/Dec_22/Dec_22.pkl'

with open(file_path_mats, 'rb') as file:
    df = pickle.load(file)

with open('config.json') as config_file:
    config = json.load(config_file)

file_path = config['osiris_file_path']
csv_file_path = config['csv_file_path']

osiris_data_files = os.listdir(file_path)
csv_data_odin = pd.read_csv(csv_file_path)

odin_ids_to_plot = csv_data_odin['Odin ID'].unique()
csv_data_odin['Time'] = pd.to_datetime(csv_data_odin['Time'])
scan_time_dict = {row['Odin ID']: row['Time'].strftime('%Y-%m-%d %H:%M') for index, row in csv_data_odin.iterrows()}

def mjd_to_datetime(mjd):
    mjd_epoch = datetime(1858, 11, 17)  # MJD epoch
    return mjd_epoch + timedelta(days=float(mjd))

rgb_color_green = (0/255, 150/255, 130/255)
rgb_color_blue = (70/255, 100/255, 170/255)
rgb_color_black = (64/255, 64/255, 64/255)

mats_color = rgb_color_blue
osiris_color = rgb_color_green

group_ids = csv_data['group_id'].unique()
for group_id in group_ids:
    group_data = csv_data[csv_data['group_id'] == group_id]
    filtered_data = []

    for index, row in group_data.iterrows():
        matching_data = df[df['EXPDate'] == row['Time']]
        
        if not matching_data.empty:
            filtered_data.extend(matching_data.values)

    filtered_df = pd.DataFrame(filtered_data, columns=df.columns)
    mats_profiles = []
    label_added = False
    
    plt.figure(figsize=(6, 5.5))  # Width = 6 inches, Height = 8 inches
    plt.rcParams.update({'font.size': 12.5})
    
    for i in range(len(filtered_df)):
        heights, profile = prepare_profile(filtered_df.iloc[i])
        mats_profiles.append(profile)
        if not label_added:
            plt.plot(profile, heights / 1000, alpha=0.4, color=rgb_color_black, label='MATS measurement')
            label_added = True  
        else:
            plt.plot(profile, heights / 1000, alpha=0.5, color=rgb_color_black) 

    mean_mats_profile = np.mean(mats_profiles, axis=0)
    # mean_mats_profile = mean_mats_profile / 10000

    plt.plot(mean_mats_profile, heights / 1000, color=mats_color, linestyle='solid', linewidth=2, label=f'Mean MATS measurements')

    group_odin_ids = csv_data_odin[csv_data_odin['group_id'] == group_id]['Odin ID'].unique()
    for scan_num in group_odin_ids:
        for i in osiris_data_files:
            if i.endswith(".nc"):
                try:
                    osiris_data = nc.Dataset(os.path.join(file_path, i))
                except OSError as e:
                    print(f"Error opening file {i}: {e}")
                    continue
                
                tangent_latitude = osiris_data.variables['tangent_latitude'][:]
                tangent_altitude = osiris_data.variables['tangent_altitude'][:] / 1000
                scan_num_time_dim = osiris_data.variables['scan_number_time_dim'][:]
                radiance = osiris_data.variables['radiance'][:]
                wavelength = osiris_data.variables['wavelength'][:]

                if scan_num not in scan_num_time_dim:
                    continue
                
                indices = np.where(scan_num_time_dim == scan_num)[0]
                radiance_NLC = radiance[:, 72:80]
                radiance_NLC_mean = np.mean(radiance_NLC, axis=1)
                rad = radiance_NLC_mean[indices]
                rad = rad*10**4
                alt = tangent_altitude[indices]
                time = scan_time_dict.get(scan_num, "Unknown Time")

                plt.plot(rad, alt, color=osiris_color, marker='', linestyle='-', linewidth=2, label=f'OSIRIS measurements')
                
    plt.tick_params(axis='both', which='major', length=8, direction='in')
    plt.tick_params(axis='both', which='minor', length=4, direction='in')
    plt.minorticks_on()
    plt.xlabel('Limb radiance\n[ph·m$^{-2}$·s$^{-1}$·sr$^{-1}$·nm$^{-1}$]')
    plt.ylabel('Tangent altitude (km)')
    # data_min = min(mean_mats_profile)
    # data_max = max(mean_mats_profile) 
    # threshold = 1e15
    # plt.xlim(data_min - threshold, data_max + threshold)
    plt.xlim(-0.5e15, 2e15) 
    plt.ylim(60, 92)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/mats_project/MATS-analysis/Aglaja/output/profils/mats_and_odin/profiles_{group_id}_{time}_feb_v5.pdf')  
    plt.clf()

