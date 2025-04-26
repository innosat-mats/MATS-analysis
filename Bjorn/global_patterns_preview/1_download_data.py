# THIS IS AN UPDATED VERSION OF SCRIPT MAIN 
# Only used to download and save weeks of data
# now uses version 5 of data

from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import numpy as np
import pandas as pd
import gc
import os

level = '1b'

# THIS IS DONE IN WEEKS NOW - TO PREVENT CRASHING

day_couples=([1,8],[8,15],[15,22],[22,1])


for month in [2, 3, 4]:
    for i in range(len(day_couples)):
        for CCDno in [2, 3]:
            if i == 3:
                print(f'Start: Month {month}, {day_couples[i][0]}; End: Month {month+1}, {day_couples[i][1]}')
                start_time = DT.datetime(2023, month, day_couples[i][0], 0, 0, 0)
                stop_time = DT.datetime(2023, month+1, day_couples[i][1], 0, 0, 0)

            else:
                print(f'Start: Month {month}, {day_couples[i][0]}; End: Month {month}, {day_couples[i][1]}')
                start_time = DT.datetime(2023, month, day_couples[i][0], 0, 0, 0)
                stop_time = DT.datetime(2023, month, day_couples[i][1], 0, 0, 0)
            
            # load data
            filter = {'CCDSEL': [CCDno, CCDno], 'TPlat': [-70, 70]}
            data = read_MATS_data(start_time, stop_time, level=level, version='0.5', filter=filter)
            
            # save data
            monthstring = f'0{month}2023_{i+1}'
            pathname = f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL{CCDno}'

            # Check whether the specified path exists or not
            isExist = os.path.exists(pathname)
            if not isExist:

                # Create a new directory because it does not exist
                os.makedirs(pathname)
                print("DIRECTORY CREATED")

            data.to_pickle(f'{pathname}/{monthstring}.pk1')
            print(f'Saved file: {pathname}/{monthstring}.pk1')
            gc.collect()