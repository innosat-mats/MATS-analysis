import os
import pandas as pd
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pickle
import json
from mats_utils.geolocation.coordinates import col_heights
import logging

logging.basicConfig(level=logging.INFO)

# Function to prepare MATS profile
def prepare_profile(ch):
    image = np.stack(ch.ImageCalibrated)
    col = int(ch['NCOL'] / 2)
    cs = col_heights(ch, col, 10, spline=True)
    heights = np.array(cs(range(ch['NROW'])))

    # Apply scaling (1e12) based on the data version
    profile = np.array(image[:, col-2:col+2].mean(axis=1) * 1e12)

    common_heights = np.arange(60000, 110250, 250)
    profile = np.interp(common_heights, heights, profile)

    return common_heights, profile

def load_data():
    with open('config.json') as config_file:
        config = json.load(config_file)

    csv_data = pd.read_csv(config['csv_file_path'])
    osiris_data_files = [f for f in os.listdir(config['osiris_file_path']) if f.endswith(".nc")]
    
    return csv_data, osiris_data_files, config

def load_osiris_data(file_path, scan_num):
    try:
        osiris_data = nc.Dataset(file_path)
    except OSError as e:
        logging.error(f"Error opening file {file_path}: {e}")
        return None

    scan_num_time_dim = osiris_data.variables['scan_number_time_dim'][:]
    
    if scan_num not in scan_num_time_dim:
        return None

    indices = np.where(scan_num_time_dim == scan_num)[0]
    radiance_NLC = osiris_data.variables['radiance'][:, 72:80]
    radiance_NLC_mean = np.mean(radiance_NLC, axis=1)
    rad = radiance_NLC_mean[indices] * 1e4
    alt = osiris_data.variables['tangent_altitude'][indices] / 1000  # Convert to km

    return rad, alt

def plot_profiles(group_id, filtered_df, mean_mats_profile, heights, osiris_profiles, save_path):
    plt.figure(figsize=(5, 6.5))
    plt.rcParams.update({'font.size': 14})
    
    # Plot MATS profiles
    for profile in filtered_df:
        plt.plot(profile, heights / 1000, alpha=0.5, color='black', label='MATS')

    # Plot mean MATS profile
    plt.plot(mean_mats_profile, heights / 1000, color='blue', linestyle='solid', linewidth=2, label=f'Mean MATS')

    # Plot OSIRIS profiles
    for rad, alt in osiris_profiles:
        plt.plot(rad, alt, color='red', linewidth=2, label='OSIRIS')

    plt.tick_params(axis='both', which='major', length=8, direction='in')
    plt.tick_params(axis='both', which='minor', length=4, direction='in')
    plt.minorticks_on()
    plt.xlabel('Limb radiance\n[ph·m$^{-2}$·s$^{-1}$·sr$^{-1}$·nm$^{-1}$]')
    plt.ylabel('Tangent altitude [km]')
    plt.xlim(-0.25e15, 2.25e15)
    plt.ylim(60, 90)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(f'{save_path}_group_{group_id}.pdf')
    plt.clf()

def process_data():
    csv_data, osiris_files, config = load_data()
    df_path = config['mats_pkl_file']
    
    with open(df_path, 'rb') as file:
        df = pickle.load(file)
    
    odin_ids_to_plot = csv_data['Odin ID'].unique()
    scan_time_dict = {row['Odin ID']: row['Time'] for index, row in csv_data.iterrows()}

    group_ids = csv_data['group_id'].unique()

    for group_id in group_ids:
        group_data = csv_data[csv_data['group_id'] == group_id]
        filtered_df = df[df['EXPDate'].isin(group_data['Time'])]

        if filtered_df.empty:
            logging.info(f"No matching data for group_id {group_id}. Skipping...")
            continue
        
        # Prepare MATS profiles
        mats_profiles = [prepare_profile(filtered_df.iloc[i])[1] for i in range(len(filtered_df))]
        mean_mats_profile = np.mean(mats_profiles, axis=0)
        common_heights = prepare_profile(filtered_df.iloc[0])[0]

        # Process OSIRIS data for the group
        osiris_profiles = []
        group_odin_ids = group_data['Odin ID'].unique()

        for scan_num in group_odin_ids:
            for file in osiris_files:
                rad, alt = load_osiris_data(os.path.join(config['osiris_file_path'], file), scan_num)
                if rad is not None and alt is not None:
                    osiris_profiles.append((rad, alt))

        plot_profiles(group_id, mats_profiles, mean_mats_profile, common_heights, osiris_profiles, config['output_folder'])

if __name__ == "__main__":
    process_data()
