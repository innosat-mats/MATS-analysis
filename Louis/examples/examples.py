##%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import pandas as pd

# times for start and stop
start_time = DT.datetime(2023, 2, 12, 6, 0, 0)
stop_time = DT.datetime(2023, 2, 12, 6, 2, 0)

# filter
# filter={'CCDSEL': [7,7]}

# read in measurments
CCDitems = read_MATS_data(start_time, stop_time, level='1a', version=0.5, filter=None)

# how many items
print(len(CCDitems))

# print all available variables
print(CCDitems.keys())

#%% showing all the keys and their type
pd.set_option('display.max_rows', 100)
CCDitems.dtypes

#%%
for i in range(len(CCDitems)):
    print("\n\n")
    print(CCDitems.iloc[i]['CCDSEL'])
    print(CCDitems.iloc[i]['flipped'])