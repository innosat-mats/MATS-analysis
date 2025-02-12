

"""
Script to Download and Plot Data for Comparison with Peter Dalin's NLC Observations 2024

This script downloads data from specified sources and generates plots to compare with 
Peter Dalin's Noctilucent Clouds (NLC) observations for the year 2024.

Author: Linda Megner
Date: 2024-02-11

Från Peter Dalin:
Jag antar att NLC syns ganska bra på bilder den 25 juni vid 10:47:14:261947 (UV1) och 10:47:15:417938 (IR1).

Och vi ser starka NLC från ballongen mellan kl. 10:00 och 12:00 UT den 25 juni som ligger över Greenland.

"""

import matplotlib.pyplot as plt
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT

# times for start and stop
start_time = DT.datetime(2024, 6, 25, 10, 47, 0)
stop_time = DT.datetime(2023, 6, 25, 10, 48, 0)

# filter
filter={'channel': ['IR1','UV1'] }

#%%
# read in measurements

df = read_MATS_data(start_time, stop_time,filter, level='1b',version='1.0')
