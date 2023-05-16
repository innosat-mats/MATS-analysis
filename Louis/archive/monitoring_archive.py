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

#%%
start_time = datetime(2023, 5, 7, 21, 0, 0)
stop_time = datetime(2023, 5, 8, 0, 0, 0)

sampling = 'custom'

custom_period = timedelta(minutes=2)

if sampling == 'day':
    start=start_time.replace(hour=0, minute=0, second=0, microsecond=0)
    end=stop_time.replace(hour=0, minute=0, second=0, microsecond=0)
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).days + 1)

elif sampling == 'orbit':
    orbit_filter = {'CCDSEL': [1,7],'satlat':[-1.0,+1.0]}
    df = read_MATS_data(start_time, stop_time,filter=orbit_filter,level='1a',version='0.5')
    df = df[~np.isnan(df['satlat'])].sort_values('EXPDate')
    time_sampling = [start_time.replace(tzinfo=timezone.utc)]
    for i in range(len(df)-1):
        if df.iloc[i]['satlat']<0.0 and df.iloc[i+1]['satlat']>0.0:
            time_sampling.append(df.iloc[i]['EXPDate'])  
    time_sampling.append(stop_time.replace(tzinfo=timezone.utc))    

if sampling == 'custom':
    start=start_time
    end=stop_time
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/custom_period.total_seconds() + 1,tz=timezone.utc)


#%%

filter={'CCDSEL': [1,7]}


df0 = []
df1a = []
df1b = []


df0 = read_MATS_data(start_time, stop_time,level='0',version='0.3')
df1a = read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
df1b = read_MATS_data(start_time, stop_time,filter,level='1b',version='0.4')

#%%

image_processing = [[df0,df1a],[df1a,df1b],[df0,df1b]]
processing_labels = ["L0-v0.3 --> L1a-v0.5","L1a-v0.5 --> L1b-v0.4","L0-v0.3 --> L1b-v0.4"]


image_processing = [[df1a,df1b]]
processing_labels = ["L1a-v0.5 --> L1b-v0.4"]




lim_red = 0
lim_orange = 0.5
lim_blue = 0.8
lim_green = 0.99999
width = 0.9

fig, ax = plt.subplots()
n = len(image_processing)
for i in range(n):
        ax.hlines(y=processing_labels[i],xmin=time_sampling[0],xmax=time_sampling[0],color='white') # some invisible line to have a working plot
        for ind_time in range(len(time_sampling)-1):
            start = time_sampling[ind_time]
            end = time_sampling[ind_time+1]
            df_origin = image_processing[i][0]
            df_processed = image_processing[i][1]
            nb_im_origin = sum((start<df_origin['EXPDate']) & (df_origin['EXPDate']<end))  
            nb_im_processed = sum((start<df_processed['EXPDate']) & (df_processed['EXPDate']<end)) 
            # print(start,end)
            # print(nb_im_origin)  
            # print(nb_im_processed)         
            processing_rate = 0.0
            color = 'white'
            if nb_im_origin > 0:
                processing_rate = nb_im_processed/nb_im_origin            
            
                if (lim_red<=processing_rate) and (processing_rate<=lim_orange):
                    color = 'red'
                elif (lim_orange<processing_rate) and (processing_rate<=lim_blue):
                    color = 'orange'
                elif (lim_blue<processing_rate) and (processing_rate<=lim_green):
                    color = 'tab:blue'
                elif (lim_green<processing_rate):
                    color = 'green'
            #ax.hlines(y=AWS_files[i],xmin=day,xmax=day+timedelta(time_sampling=1),color=color)       
            ax.add_patch(Rectangle((start,i-width*0.5),end-start,width,color=color))

# legend
legend_elements = [Patch(facecolor='red',label=f"{lim_red*100:.0f}% <= success rate <= {lim_orange*100:.0f}%"),
                   Patch(facecolor='orange',label=f"{lim_orange*100:.0f}% < success rate <= {lim_blue*100:.0f}%"),
                   Patch(facecolor='tab:blue',label=f"{lim_blue*100:.0f}% < success rate < 100%"),
                   Patch(facecolor='green',label=f"success rate = 100%")]
ax.set_xlabel("Parquet date")
ax.set_title(f"Daily processing success rate | filter:{filter}")
ax.legend(handles=legend_elements,loc='upper left')
plt.show()





#%%
# given a version, how many images are produced vs the number that was expected knowing the sampling rate

channels = {1:'IR1',2:'IR4',3:'IR3',4:'IR2',5:'UV1',6:'UV2',7:'NADIR'}
sampling_rates={'IR1':timedelta(seconds=6),
                'IR2':timedelta(seconds=6),
                'IR3':timedelta(seconds=6),
                'IR4':timedelta(seconds=6),
                'UV1':timedelta(seconds=6),
                'UV2':timedelta(seconds=6),
                'NADIR':timedelta(seconds=2)}

lim_red = 0
lim_orange = 0.5
lim_blue = 0.8
lim_green = 0.99
width = 0.9

start_sza = 97 # sza for which nadir measurement starts
stop_sza = 95 # sza for which nadir measurement ends

df = df1a.sort_values('EXPDate') # dataframe being studied
df_loc = df1a # dataframe from the same period with geolocation 

fig, ax = plt.subplots()
# iteration over channels
for channel_ind in range(1,8):
    
    channel = channels[channel_ind]

    # compute the number of expected images in each time intervall
    nb_expected_images = []
    if channel in ['IR1','IR2','IR3','IR4','UV1','UV2']:
        for i in range(len(time_sampling)-1):
            #print((time_sampling[i+1]-time_sampling[i]).total_seconds()/sampling_rates[channel].total_seconds())
            nb_expected_images.append((time_sampling[i+1]-time_sampling[i])/sampling_rates[channel])

    if channel == 'NADIR':
        # list having the start and end times of all the NADIR measurement windows
        nadir_measurements = [] 
        try :
            df_loc = df1a.sort_values('EXPDate')
            orb_start = df_loc.iloc[0]['EXPDate']
            orb_end = df_loc.iloc[0]['EXPDate']
            for i in range(len(df_loc)-1):
                if df_loc.iloc[i]['nadir_sza']<start_sza and df_loc.iloc[i+1]['nadir_sza']>start_sza:
                    orb_start = df_loc.iloc[i+1]['EXPDate']             
                elif df_loc.iloc[i]['nadir_sza']>stop_sza and df_loc.iloc[i+1]['nadir_sza']<stop_sza:
                    orb_end = df_loc.iloc[i]['EXPDate']
                    nadir_measurements.append([orb_start,orb_end])
            if orb_end < orb_start:
                orb_end = df_loc.iloc[-1]['EXPDate']
                nadir_measurements.append([orb_start,orb_end])
        except :
            nadir_measurements = [[time_sampling[0],time_sampling[-1]]]
        
        # compute number of expected nadir images
        for i in range(len(time_sampling)-1):            
            nadir_duration = timedelta(seconds=0)
            start = time_sampling[i]
            end = time_sampling[i+1]
            for j in range(len(nadir_measurements)):
                start_win = nadir_measurements[j][0]
                end_win = nadir_measurements[j][1]
                if start_win<end and end_win>=start:
                    nadir_duration += min(end-start_win,end_win-start,end_win-start_win,end-start)
           
            nb_expected_images.append(nadir_duration/sampling_rates[channel])
            
    # if channel in ['UV1','UV2']:
    #     nb_expected_images = np.zeros(len(time_sampling)-1)

    
    # compute the ratio nb of images/expected number of images

    ax.hlines(y=channel,xmin=time_sampling[0],xmax=time_sampling[0],color='white') # some invisible line to have a working plot
    for i in range(len(time_sampling)-1):
        start = time_sampling[i]
        end = time_sampling[i+1]
        color = 'white'
        if nb_expected_images[i] > 0:
            ratio = len(df[(df['CCDSEL']==channel_ind) & (start<=df['EXPDate']) & (df['EXPDate']<end)])/nb_expected_images[i]
            
            if ratio == 0.0:
                color ='black'
            elif (lim_red<ratio) and (ratio<=lim_orange):
                color = 'red'
            elif (lim_orange<ratio) and (ratio<=lim_blue):
                color = 'orange'
            elif (lim_blue<ratio) and (ratio<=lim_green):
                color = 'tab:blue'
            elif (lim_green<ratio):
                color = 'green'
            #print(ratio)    
        ax.add_patch(Rectangle((start,channel_ind-1-width*0.5),end-start,width,color=color))

# legend
legend_elements = [Patch(facecolor='black',label=f"0 == ratio"),    
                   Patch(facecolor='red',label=f"{lim_red*100:.0f}% <= ratio <= {lim_orange*100:.0f}%"),
                   Patch(facecolor='orange',label=f"{lim_orange*100:.0f}% < ratio <= {lim_blue*100:.0f}%"),
                   Patch(facecolor='tab:blue',label=f"{lim_blue*100:.0f}% < ratio <= {lim_green*100:.1f}%"),
                   Patch(facecolor='green',label=f"{lim_green*100:.1f}% < ratio")]
ax.set_xlabel("Date")
ax.set_title("nb of images/expected nb of images (L0 v0.3)")
ax.legend(handles=legend_elements,loc='upper left')
plt.show()