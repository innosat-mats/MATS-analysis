import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import argparse
from matplotlib import pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import col_heights, TPpos
from mats_l1_processing.pointing import pix_deg
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
import pdb

dftop = pd.read_pickle('NLC_test_file.pkl')
df = dftop[dftop['channel'] == 'UV2'].dropna().reset_index(drop=True)

# TPlat,TPlon,TPheight = TPpos(df.iloc[1])

lst_lat = []
lst_lon = []
lst_alt = []
lst_time = []
lst_id = []

for i in range(len(df)):
    TPlat,TPlon,TPheight = TPpos(df.iloc[i])
    lst_lat.append(TPlat)
    lst_lon.append(TPlon)
    lst_alt.append(TPheight)
    time = df['EXPDate'].iloc[i]
    lst_time.append(time)
    lst_id.append(i)
'''
Returns:
        TPlat: latitude of satellite (degrees) (Latitude of tangent point at time of measurement [-90, 90])
        TPlon: longitude of satellite (degrees) (Longitude of tangent point at time of measurement, east of Greenwich meridian [0,360])
        TPheight: Altitude in metres (Altitude of tangent point at time of measurement)
        Time: EXPDate --> Time of exposure (yyyy-mm-dd hh:mm:ss.ss). (Timestamp)
'''
# Time
# df['EXPDate']
midday=DT.time(12, 0, 0)
df['EXPDate'].dt.time > midday

# Write CSV Table out of the Lists (include image number)
df_out = pd.DataFrame(list(zip(lst_id, lst_time, lst_lat, lst_lon, lst_alt)), columns =['Image ID', 'Time', 'Latitude', 'Longitude', 'Altitude'])
df_out.to_csv('NLC_test_file_out.csv', index=False)

pdb.set_trace()