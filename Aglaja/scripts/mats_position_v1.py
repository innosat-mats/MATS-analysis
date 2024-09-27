import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import os
from mats_utils.geolocation.coordinates import col_heights, TPpos
import glob
from joblib import Parallel, delayed

'''
- find all the TPs of the MATS satellite to compare later to Odin
- use image id, time, latitude, longitutde and altitude 
- save the data in a csv file

TPpos Returns:
        TPlat: latitude of satellite (degrees) (Latitude of tangent point at time of measurement [-90, 90])
        TPlon: longitude of satellite (degrees) (Longitude of tangent point at time of measurement, east of Greenwich meridian [0,360])
        TPheight: Altitude in metres (Altitude of tangent point at time of measurement)
        Time: EXPDate --> Time of exposure (yyyy-mm-dd hh:mm:ss.ss). (Timestamp)

'''

def process_pkl_file(file_path, output_csv_path):
    df = pd.read_pickle(file_path)
    df = df[df['channel'] == 'UV2'].dropna().reset_index(drop=True)

    TPlat, TPlon, TPheight = zip(*df.apply(TPpos, axis=1))
    lst_time = df['EXPDate'].tolist()
    lst_id = list(range(len(df)))

    df_out = pd.DataFrame({
        'Image ID': lst_id,
        'Time': lst_time,
        'Latitude': TPlat,
        'Longitude': TPlon,
        'Altitude': TPheight
    })

    df_out.to_csv(output_csv_path, index=False)
    print(f"Processed {file_path} and saved to {output_csv_path}")

def process_multiple_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    pkl_files = glob.glob(os.path.join(input_folder, "*.pkl"))

    Parallel(n_jobs=-1)(delayed(process_pkl_file)(pkl_file, os.path.join(output_folder, f"{os.path.splitext(os.path.basename(pkl_file))[0]}_mats.csv")) for pkl_file in pkl_files)

input_folder = '../output/images/all/'  
output_folder = '../output/mats_location/all/'  

process_multiple_files(input_folder, output_folder)
