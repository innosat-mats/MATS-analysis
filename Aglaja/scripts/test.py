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

# Criteria
'''
- SH 22/23 (Nov-March)
start: 2022 11 1 0 0
stop: 2023 3 31 0 0
python get_data.py --start_time 2022 11 1 0 0 --stop_time 2023 3 31 0 0 --version 0.6 --channel UV2 --ncdf_out mats_data.nc

- UV channel
- latitude < -40°
ztan= 70-10 km

'''

# for testing: 
#starttime=datetime(2023,1,8,19,0)
#stoptime=datetime(2023,1,8,19,0)


# load images
channel='UV'
starttime=datetime(2023,2,11,0,0)
stoptime=datetime(2023,2,12,0,0)


df=read_MATS_data(starttime,stoptime,level="1b", version="0.6")

df.columns

# plot iamge
plot_image(df.iloc[3], save=False)
plt.show(block=True)

'''
simple_plot(df, './', custom_cbar=True, ranges=[500, 20000])
orbit_plot(df, './')
'''



pdb.set_trace()
