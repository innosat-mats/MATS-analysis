# THIS IS AN UPDATED VERSION OF ADD ALTITUDES
# NOW ON V0.5 data and WEEK BY WEEK

import datetime as DT
import numpy as np
import pandas as pd
import gc
from mats_utils.geolocation.altitude_correction import rows_to_altitudes
import sys
import multiprocessing
import os

# THIS IS DONE IN WEEKS NOW - TO PREVENT CRASHING
day_couples=([1,8],[8,15],[15,22],[22,1])

# ALTITUDES
fixaltitudes=np.array([80,81,82,83,84,85,86,87,88,89,90])*1000

# ASCENDING/DESCENDING 
ascending = True


# FUNCTIONS -------------------
def generate_heights_save(CCDs,part):

    for element in range(0, len(CCDs)):
        CCDs.ImageCalibrated.iloc[element]=rows_to_altitudes(CCDs.iloc[element], fixaltvec=fixaltitudes, imagefield='ImageCalibrated')

    if ascending:
        CCDs.to_pickle(f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL{CCDno}/alts_80to90/{monthstring}/{monthstring}_part{part}.pk1')
    else:
        CCDs.to_pickle(f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL{CCDno}/alts_80to90/{monthstring}_desc/{monthstring}_part{part}.pk1')

def parallel_plotting(part):

    files_per_part = 500

    if int(len(CCDs)) > files_per_part:
            
        if part == 0:
            start_point = 0
        else:
            start_point = part*files_per_part-1
        try:
            if (part+1)*files_per_part < int(len(CCDs)):
                generate_heights_save(CCDs[start_point:(part+1)*files_per_part-1],part)
            else:
                generate_heights_save(CCDs[start_point:int(len(CCDs))-1],part)
        except KeyboardInterrupt:
            sys.exit()
    else:
        try:
            generate_heights_save(CCDs, part)
        except KeyboardInterrupt:
            sys.exit()


for month in [3, 4, 5]:
    for day_couple in range(len(day_couples)):
        for CCDno in [2, 3]:
 
            # load data
            monthstring = f'0{month}2023_{day_couple+1}'
            pathname = f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL{CCDno}'

            CCDs=pd.read_pickle(f'{pathname}/{monthstring}.pk1')

            # GET RID OF FULL-FRAMES
            CCDs = CCDs[(CCDs['NCOL'] < 1000)]
            CCDs = CCDs[(CCDs['NROW'] < 300)]

            # local time stuff
            midday=DT.time(12, 0, 0)
            CCDs['TPlocaltime'] = pd.to_datetime(CCDs['TPlocaltime'],utc=True)
            CCDs['EXPDate'] = pd.to_datetime(CCDs['EXPDate'],utc=True)
            for elem in range(0,len(CCDs.EXPDate)):
                CCDs['EXPDate'].iloc[elem] = CCDs['EXPDate'].iloc[elem].to_pydatetime()

            if ascending:
                CCDs = CCDs[(CCDs['TPlocaltime'].dt.time > midday)]
            else:
                CCDs = CCDs[(CCDs['TPlocaltime'].dt.time < midday)]

            # to check path
            pathname = f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL{CCDno}/alts_80to90/{monthstring}'

            # Check whether the specified path exists or not
            isExist = os.path.exists(pathname)
            if not isExist:
                # Create a new directory because it does not exist
                os.makedirs(pathname)
                print("! DIRECTORY CREATED !")

            print(f'-------- Adding altitudes to {monthstring} (CCD{CCDno}) --------')
            
            # parallel processing stuff
            files_per_part = 500
            sets = int(np.floor(len(CCDs)/files_per_part))
            parts = list(np.arange(0, sets))
            pool = multiprocessing.Pool(13)
            pool.map(parallel_plotting, parts)
            CCDs=[]

            gc.collect()