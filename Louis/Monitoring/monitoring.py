#%% Import modules
#%matplotlib qt5
import datetime as DT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from typing import Any, Optional, Union
from datetime import datetime, date

import boto3
import re
import pyarrow.parquet as pq  # type: ignore
from pyarrow import fs
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch


#%%

# defining study window 

start_time = DT.datetime(2022, 12, 1, 0, 0, 0)
stop_time = DT.datetime(2023, 5, 10, 0, 0, 0)

start=start_time.replace(hour=0, minute=0, second=0, microsecond=0)
end=stop_time.replace(hour=0, minute=0, second=0, microsecond=0)

days = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).days + 1)


#%%
# defining usefull functions
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





#%%
# plotting all the parquet files dates
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
    #print(keys[6])
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

# number of parquet files per day
AWS_files = ['ops-payload-level0-v0.3',
             'ops-payload-level1a-v0.5',
             'ops-payload-level1b-v0.4']


lim_red = 0
lim_orange = 5
lim_blue = 20
lim_green = 23
width = 0.9

s3_client = boto3.client('s3')

fig, ax = plt.subplots()
for i in range(len(AWS_files)):
    print(AWS_files[i])
    day_parquet = np.zeros(len(days))
    keys = get_object_keys(s3_client,bucket=AWS_files[i]) 
    dates = parse_object_date(keys)   
    print(f"{len(keys)} parquet files ({len(dates)} with dates)")
    ax.hlines(y=AWS_files[i],xmin=days[0],xmax=days[0],color='white')
    if type(dates) != type(None): 
        for ind_day in range(len(days)):
            day = days[ind_day]
            for j in range(len(dates)):
                if dates[j]>day and dates[j]<day + timedelta(days=1):
                    day_parquet[ind_day] += 1
            nb_files = day_parquet[ind_day]
            color = 'white'
            if (lim_red<nb_files) and (nb_files<=lim_orange):
                color = 'red'
            elif (lim_orange<nb_files) and (nb_files<=lim_blue):
                color = 'orange'
            elif (lim_blue<nb_files) and (nb_files<=lim_green):
                color = 'tab:blue'
            elif (lim_green<nb_files):
                color = 'green'      
            ax.add_patch(Rectangle((day,i-width*0.5),timedelta(days=1),width,color=color))

legend_elements = [Patch(facecolor='red',label=f"{lim_red:.0f} < files/day <= {lim_orange:.0f}"),
                   Patch(facecolor='orange',label=f"{lim_orange:.0f} < files/day <= {lim_blue:.0f}"),
                   Patch(facecolor='tab:blue',label=f"{lim_blue:.0f} < files/day <= {lim_green:.0f}"),
                   Patch(facecolor='green',label=f"24 files/day")]        
        
ax.set_xlabel("Parquet date")
ax.set_title("Number of .parquet files generated per day")
ax.legend(handles=legend_elements,loc='upper left')
plt.show()


#%%

# checking the daily processing rate 

AWS_processing = np.array([['ops-payload-level0-v0.3','ops-payload-level1a-v0.5'],['ops-payload-level1a-v0.5','ops-payload-level1b-v0.4'],['ops-payload-level0-v0.3','ops-payload-level1b-v0.4']])
AWS_labels = ["L0-v0.3 --> L1a-v0.5","L1a-v0.5 --> L1b-v0.4","L0-v0.3 --> L1b-v0.4"]

lim_red = 0
lim_orange = 0.5
lim_blue = 0.8
lim_green = 0.999
width = 0.9

fig, ax = plt.subplots()
n = np.shape(AWS_processing)[0]
for i in range(n):
    origin = AWS_processing[i,0] # bucket with the files before processing
    processed = AWS_processing[i,1] # bucket with the processed files
    day_parquet_origin = np.zeros(len(days)) # number of parquet files per day before processing
    print(origin)    
    keys = get_object_keys(s3_client,bucket=origin) 
    dates = parse_object_date(keys)   
    print(f"{len(keys)} parquet files ({len(dates)} with dates)")
    
    if type(dates) != type(None): 
        for ind_day in range(len(days)):
            day = days[ind_day]
            nb_files = 0
            for j in range(len(dates)):
                if dates[j]>day and dates[j]<day + timedelta(days=1):
                    nb_files += 1
            day_parquet_origin[ind_day] = nb_files

    day_parquet_processed = np.zeros(len(days)) # number of parquet files per day after processing
    print(processed)    
    keys = get_object_keys(s3_client,bucket=processed) 
    dates = parse_object_date(keys)   
    print(f"{len(keys)} parquet files ({len(dates)} with dates)")
    ax.hlines(y=AWS_labels[i],xmin=days[0],xmax=days[0],color='white') # some invisible line to have a working plot
    if type(dates) != type(None): 
        for ind_day in range(len(days)):
            day = days[ind_day]
            nb_files = 0
            processing_rate = 0.0
            for j in range(len(dates)):
                if dates[j]>day and dates[j]<day + timedelta(days=1):
                    nb_files += 1
            day_parquet_processed[ind_day] = nb_files
            if day_parquet_origin[ind_day] != 0:
                processing_rate = day_parquet_processed[ind_day]/day_parquet_origin[ind_day]
            color = 'white'
            if (lim_red<processing_rate) and (processing_rate<=lim_orange):
                color = 'red'
            elif (lim_orange<processing_rate) and (processing_rate<=lim_blue):
                color = 'orange'
            elif (lim_blue<processing_rate) and (processing_rate<=lim_green):
                color = 'tab:blue'
            elif (lim_green<processing_rate):
                color = 'green'
            #ax.hlines(y=AWS_files[i],xmin=day,xmax=day+timedelta(days=1),color=color)       
            ax.add_patch(Rectangle((day,i-width*0.5),timedelta(days=1),width,color=color))

# legend
legend_elements = [Patch(facecolor='red',label=f"{lim_red*100:.0f}% < success rate <= {lim_orange*100:.0f}%"),
                   Patch(facecolor='orange',label=f"{lim_orange*100:.0f}% < success rate <= {lim_blue*100:.0f}%"),
                   Patch(facecolor='tab:blue',label=f"{lim_blue*100:.0f}% < success rate < 100%"),
                   Patch(facecolor='green',label=f"success rate = 100%")]
ax.set_xlabel("Parquet date")
ax.set_title("Daily processing success rate")
ax.legend(handles=legend_elements,loc='upper left')
plt.show()




#%%

# session = boto3.session.Session(profile_name="mats")
# credentials = session.get_credentials()

# s3 = fs.S3FileSystem(
#         secret_key=credentials.secret_key,
#         access_key=credentials.access_key,
##        region=session.region_name,
#         connect_timeout=10,
#         session_token=credentials.token)


# meta = pq.read_metadata

# #%%
# for i in range(len(AWS_files)):
#     plt.figure()
#     print(AWS_files[i])
#     day_data = np.zeros(len(days))
#     day_parquet = np.zeros(len(days))
#     keys = get_object_keys(s3_client,bucket=AWS_files[i]) 
#     dates = parse_object_date(keys)   
#     sizes = np.array(get_object_size(s3_client,bucket=AWS_files[i]))
#     print(f"{len(keys)} parquet files ({len(dates)} with dates)")

#     if type(dates) != type(None) and len(sizes)==len(dates): 
#         for ind_day in range(len(days)):
#             day = days[ind_day]
#             for j in range(len(dates)):
#                 if dates[j]>day and dates[j]<day + timedelta(days=1):
#                     day_parquet[ind_day] += 1
    
#     if type(dates) != type(None) and len(sizes)==len(dates): 
#         for ind_day in range(len(days)):
#             day = days[ind_day]
#             for j in range(len(dates)):
#                 if dates[j]>day and dates[j]<day + timedelta(days=1):
#                     day_data[ind_day] += sizes[j]
#         plt.plot(days[(lim_white<day_parquet)&(day_parquet<=lim_red)],day_data[(lim_white<day_parquet)&(day_parquet<=lim_red)],label=AWS_files[i],color='red',marker='o',linestyle='')
#         plt.plot(days[(lim_red<day_parquet)&(day_parquet<=lim_orange)],day_data[(lim_red<day_parquet)&(day_parquet<=lim_orange)],color='orange',marker='o',linestyle='')
#         plt.plot(days[(lim_orange<day_parquet)&(day_parquet<=lim_blue)],day_data[(lim_orange<day_parquet)&(day_parquet<=lim_blue)],color='tab:blue',marker='o',linestyle='')
#         plt.plot(days[lim_green<day_parquet],day_data[lim_green<day_parquet],color='green',marker='o',linestyle='')
        
#     plt.xlabel("Parquet date")
#     plt.ylabel("Data size")
#     plt.title(AWS_files[i])
#     plt.show()

#     plt.figure()
#     plt.title(AWS_files[i])
#     plt.hist(day_parquet,30)
#     plt.show()

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

# %%
