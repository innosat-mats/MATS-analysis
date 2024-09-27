import numpy as np
import pandas as pd
from datetime import datetime, timezone
import datetime as DT
import pickle
import argparse
from matplotlib import pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
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

# load images
channel='UV'
starttime=datetime(2023,2,11,19,0)
stoptime=datetime(2023,2,11,20,0)
df=read_MATS_data(starttime,stoptime, level="1b", version="0.6", filter={"CCDSEL":[6,6], "TPlat": [-90,-40]}) # df.shape[0] --> number of rows; df.shape[1] --> number of columns

