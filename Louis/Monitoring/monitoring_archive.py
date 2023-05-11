#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from mats_utils.statistiscs.images_functions import create_imagecube
from tqdm import tqdm
from mats_utils.geolocation.coordinates import NADIR_geolocation 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import warnings


#%%

start_time = DT.datetime(2023, 4, 29, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 29, 6, 0, 0)

filter={'CCDSEL': [1,7]}

df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')

print(len(df1a_tot))
df1a = df1a_tot





#%%
colors = ['b','g','r','m','orange','purple','cyan']
channels = {1:'IR1',2:'IR4',3:'IR3',4:'IR2',5:'UV1',6:'UV2',7:'NADIR'}
sampling_rates={'IR1':timedelta(seconds=6),
                'IR2':timedelta(seconds=6),
                'IR3':timedelta(seconds=6),
                'IR4':timedelta(seconds=6),
                'UV1':timedelta(seconds=6),
                'UV2':timedelta(seconds=6),
                'NADIR':timedelta(seconds=2)}





plt.figure()
for i in range(1,8):
    channel = channels[i]
    dates = df1a[df1a['CCDSEL']==i]['EXPDate']
    
    plt.scatter(dates[::10],i*np.ones_like(dates[::10]),color=colors[i-1],marker='o',label=channels[i])
    for j in range(len(dates)-1):
        if dates.iloc[j+1]-dates.iloc[j] > sampling_rates[channel]:
            plt.scatter(dates.iloc[j]+sampling_rates[channel],i,color='black',marker='o')
plt.xlabel('EXPDate')
plt.ylabel('CCDSEL')
plt.legend()
plt.show()

plt.figure()
for i in range(1,8):
    SZA = df1a[df1a['CCDSEL']==i]['nadir_sza']
    dates = df1a[df1a['CCDSEL']==i]['EXPDate']
    channel = channels[i]
    plt.scatter(SZA[::10],i*np.ones_like(SZA[::10]),color=colors[i-1],marker='o',label=channels[i])
    for j in range(len(dates)-1):
        if dates.iloc[j+1]-dates.iloc[j] > sampling_rates[channel]:
            plt.scatter(SZA.iloc[j],i,color='black',marker='o')
plt.xlabel('nadir_sza')
plt.ylabel('CCDSEL')
plt.legend()
plt.show()









# %%

nb_expected_images = {}
dropped_im_ratio = {}


# slice orbite by orbit

df = df1a[df1a['CCDSEL']==7]
orb_times = []
orb_start = df.iloc[0]['EXPDate']
orb_end = df.iloc[0]['EXPDate']
for i in range(len(df)-1):
    if df.iloc[i+1]['EXPDate']-df.iloc[i]['EXPDate'] > timedelta(seconds = 100):
        orb_end = df.iloc[i]['EXPDate']
        orb_times.append([orb_start,orb_end])
        orb_start = df.iloc[i+1]['EXPDate']
orb_times.append([orb_start,df.iloc[-1]['EXPDate']])


for channel in ['IR1','IR2','IR3','IR4','UV1','UV2','NADIR']:
    if channel in ['IR1','IR2','IR3','IR4']:
        nb_expected_images[channel] = (stop_time - start_time)/sampling_rates[channel]
        dropped_im_ratio[channel] = 1-len(df1a[df1a['channel']==channel])/nb_expected_images[channel]
    if channel == 'NADIR':
        nadir_measurements = []
        start_sza = 97
        stop_sza = 95
        df = df1a.sort_values('EXPDate')
        orb_start = df.iloc[0]['EXPDate']
        orb_end = df.iloc[0]['EXPDate']
        for i in range(len(df)-1):
            if df.iloc[i]['nadir_sza']<start_sza and df.iloc[i+1]['nadir_sza']>start_sza:
                orb_start = df.iloc[i]['EXPDate']             
            elif df.iloc[i]['nadir_sza']>stop_sza and df.iloc[i+1]['nadir_sza']<stop_sza:
                orb_end = df.iloc[i]['EXPDate']
                nadir_measurements.append([orb_start,orb_end])
        if orb_end < orb_start:
            orb_end = df.iloc[-1]['EXPDate']
            nadir_measurements.append([orb_start,orb_end])

        nadir_duration = timedelta(seconds=0)
        for i in range(len(nadir_measurements)):
            nadir_duration = nadir_duration + (nadir_measurements[i][1]-nadir_measurements[i][0])
        nb_expected_images[channel] = nadir_duration/sampling_rates[channel]
        dropped_im_ratio[channel] = 1-len(df1a[df1a['channel']==channel])/nb_expected_images[channel]

           
        
            
    



# %%
import boto3

s3 = boto3.resource('s3')
bucket = s3.Bucket('ops-payload-level1a-v0.5')
total_size = sum([obj.size for obj in bucket.objects.all()])
print(total_size)


# %%
from typing import Any, Optional, Union
from datetime import datetime, date

import boto3
import re


BotoClient = Any


def get_object_keys(
    s3_client: BotoClient,
    bucket: str,
    prefix: str = "",
) -> list[str]:
    """List objects and return only keys"""

    results = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix)
    keys = [f["Key"] for f in results["Contents"]]

    # AWS sometimes returns partial data, loop until we get all:
    while results["IsTruncated"]:
        results = s3_client.list_objects_v2(
            Bucket=bucket,
            ContinuationToken=results["NextContinuationToken"],
            Prefix=prefix,
        )
        keys.extend([f["Key"] for f in results["Contents"]])

    return keys


def get_object_size(
    s3_client: BotoClient,
    bucket: str,
    prefix: str = "",
) -> list[str]:
    """List objects and return only keys"""

    results = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix)
    sizes = [f["Size"] for f in results["Contents"]]

    # AWS sometimes returns partial data, loop until we get all:
    while results["IsTruncated"]:
        results = s3_client.list_objects_v2(
            Bucket=bucket,
            ContinuationToken=results["NextContinuationToken"],
            Prefix=prefix,
        )
        sizes.extend([f["Size"] for f in results["Contents"]])

    return sizes


def get_object_date(
    s3_client: BotoClient,
    bucket: str,
    prefix: str = "",
) -> list[str]:
    """List objects and return only keys"""

    results = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix)
    keys = [f["LastModified"] for f in results["Contents"]]

    # AWS sometimes returns partial data, loop until we get all:
    while results["IsTruncated"]:
        results = s3_client.list_objects_v2(
            Bucket=bucket,
            ContinuationToken=results["NextContinuationToken"],
            Prefix=prefix,
        )
        keys.extend([f["LastModified"] for f in results["Contents"]])

    return keys




def parse_object_keys(
    objects: list[str],
    prefix: str = "",
) -> Union[list[datetime], list[date]]:
    "Parse object keys to dates or datetimes"

    # Strip prefix and '/' if it exist:
    if (length := len(prefix)) > 0:
        objects = [o[length+1:] for o in objects]



    try:
        # First try hourly partitioning ...
        return [
            datetime.strptime(o.split("/MATS")[0], "%Y/%M/%d/%H")
            for o in objects
        ]
    except ValueError:
        # ... then try daily
        return [
            datetime.strptime(o.split("/MATS")[0], "%Y/%M/%d").date()
            for o in objects
        ]
    


def parse_object_date(
    objects: list[str],
    prefix: str = "",
) -> Union[list[datetime], list[date]]:
    "Parse object keys to dates or datetimes"

    # Strip prefix and '/' if it exist:
    if (length := len(prefix)) > 0:
        objects = [o[length+1:] for o in objects]

    dates = []

    for o in objects:
        if len(o)>9 and o[-8:] == '.parquet':
            # Define a regex pattern for matching dates in the format YYYY/MM/DD
            date_pattern = re.compile(r"\d{4}/\d{1,2}/\d{1,2}(?:/\d{1,2})?")

            # Search for the pattern in the string
            date_match = re.search(date_pattern, o)

            # If a match is found, extract the date
            if date_match:
                date_str = date_match.group()  
                #print(date_str)
                try:
                    date = datetime.strptime(date_str, "%Y/%m/%d/%H")
                except:
                    date = datetime.strptime(date_str, "%Y/%m/%d")
                if date > DT.datetime(2000,1,1):
                    dates.append(date)           
            
    
    return dates


# %%

# AWS_files = ['ops-payload-level0-source',
#              'ops-payload-level0-v0.1',
#              'ops-payload-level0-v0.2',
#              'ops-payload-level0-v0.3',
#              'ops-payload-level1a-v0.1',
#              'ops-payload-level1a-v0.2',
#              'ops-payload-level1a-v0.3',
#              'ops-payload-level1a-v0.4',
#              'ops-payload-level1a-v0.5',
#              'ops-payload-level1b-v0.3',
#              'ops-payload-level1b-v0.4',
#              'ops-payload-level1a-pm-v0.1',
#              'ops-payload-level1a-pm-v0.2']

# s3_client = boto3.client('s3')

# plt.figure()
# for i in range(len(AWS_files)):
#     print(AWS_files[i])
#     dates = get_object_date(s3_client,bucket=AWS_files[i])   
#     plt.plot(dates,i*np.ones(len(dates)),'o',label=AWS_files[i])
# plt.xlabel("Last modification date")
# plt.legend()
# plt.show()


# AWS_files = ['ops-platform-level1a-source',
#              'ops-platform-level1a-source-v0.1',
#              'ops-platform-level1a-v0.1',
#              'ops-platform-level1a-v0.2',
#              'ops-platform-level1a-v0.3']

# s3_client = boto3.client('s3')

# plt.figure()
# for i in range(len(AWS_files)):
#     print(AWS_files[i])
#     dates = get_object_date(s3_client,bucket=AWS_files[i])   
#     plt.plot(dates,i*np.ones(len(dates)),'o',label=AWS_files[i])
# plt.xlabel("Last modification date")
# plt.legend()
# plt.show()




#%%

AWS_files = ['ops-payload-level0-source',
             'ops-payload-level0-v0.1',
             'ops-payload-level0-v0.2',
             'ops-payload-level0-v0.3',
             'ops-payload-level1a-v0.1',
             'ops-payload-level1a-v0.2',
             'ops-payload-level1a-v0.3',
             'ops-payload-level1a-v0.4',
             'ops-payload-level1a-v0.5',
             'ops-payload-level1b-v0.3',
             'ops-payload-level1b-v0.4',
             'ops-payload-level1a-pm-v0.1',
             'ops-payload-level1a-pm-v0.2']

s3_client = boto3.client('s3')

plt.figure()
for i in range(len(AWS_files)):
    print(AWS_files[i])
    keys = get_object_keys(s3_client,bucket=AWS_files[i]) 
    dates = parse_object_date(keys)   
    print(f"{len(keys)} parquet files ({len(dates)} with dates)")
    if type(dates) != type(None): 
        plt.plot(dates,i*np.ones(len(dates)),'o',label=AWS_files[i])
plt.xlabel("Parquet date")
plt.title("ops payload")
plt.legend()
plt.show()



AWS_files = ['ops-platform-level1a-source',
             'ops-platform-level1a-source-v0.1',
             'ops-platform-level1a-v0.1',
             'ops-platform-level1a-v0.2',
             'ops-platform-level1a-v0.3']

s3_client = boto3.client('s3')

plt.figure()
for i in range(len(AWS_files)):
    print(AWS_files[i])
    keys = get_object_keys(s3_client,bucket=AWS_files[i]) 
    dates = parse_object_date(keys)     
    print(f"{len(keys)} parquet files ({len(dates)} with dates)") 
    if type(dates) != type(None): 
        plt.plot(dates,i*np.ones(len(dates)),'o',label=AWS_files[i])
plt.xlabel("Parquet date")
plt.title("ops platform")
plt.legend()
plt.show()


#%%

AWS_files = ['ops-payload-level0-source',
             'ops-payload-level0-v0.1',
             'ops-payload-level0-v0.2',
             'ops-payload-level0-v0.3',
             'ops-payload-level1a-v0.1',
             'ops-payload-level1a-v0.2',
             'ops-payload-level1a-v0.3',
             'ops-payload-level1a-v0.4',
             'ops-payload-level1a-v0.5',
             'ops-payload-level1b-v0.3',
             'ops-payload-level1b-v0.4',
             'ops-payload-level1a-pm-v0.1',
             'ops-payload-level1a-pm-v0.2']

s3_client = boto3.client('s3')

plt.figure()
for i in range(len(AWS_files)):
    print(AWS_files[i])
    keys = get_object_keys(s3_client,bucket=AWS_files[i]) 
    dates = parse_object_date(keys)   
    sizes = np.array(get_object_size(s3_client,bucket=AWS_files[i]))

    print(keys[0])
    print(f"{len(keys)} parquet files ({len(dates)} with dates)")
    if type(dates) != type(None) and len(sizes)==len(dates): 
        plt.plot(dates,i*np.ones(len(dates))+sizes/np.max(sizes),'o',label=AWS_files[i])
plt.xlabel("Parquet date")
plt.title("ops payload")
plt.legend()
plt.show()

# %%
AWS_files = ['ops-payload-level0-source',
             'ops-payload-level0-v0.1',
             'ops-payload-level0-v0.2',
             'ops-payload-level0-v0.3',
             'ops-payload-level1a-v0.1',
             'ops-payload-level1a-v0.2',
             'ops-payload-level1a-v0.3',
             'ops-payload-level1a-v0.4',
             'ops-payload-level1a-v0.5',
             'ops-payload-level1b-v0.3',
             'ops-payload-level1b-v0.4',
             'ops-payload-level1a-pm-v0.1',
             'ops-payload-level1a-pm-v0.2',
             'ops-platform-level1a-source',
             'ops-platform-level1a-source-v0.1',
             'ops-platform-level1a-v0.1',
             'ops-platform-level1a-v0.2',
             'ops-platform-level1a-v0.3']