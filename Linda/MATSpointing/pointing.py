#%%
from pyarrow import fs, schema, string
import pyarrow.dataset as ds
import boto3
import numpy as np
import pandas as pd
from datetime import datetime, timezone,timedelta
#import Geoidlib
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
import io
#from tangentlib import *
from PIL import Image
import mats_utils.geolocation.coordinates as coordinates
#matplotlib widget

#%%
session = boto3.session.Session(profile_name="mats")
credentials = session.get_credentials()

s3 = fs.S3FileSystem(
    secret_key=credentials.secret_key,
    access_key=credentials.access_key,
    region=session.region_name,
    session_token=credentials.token)

dataset = ds.dataset(
    "ops-platform-level1a-v0.3/ReconstructedData",
    filesystem=s3,
    )

table = dataset.to_table(
    filter=(
        ds.field('time') > pd.to_datetime('2023-02-01T0:0:00z').to_datetime64()
    ) & (
        ds.field('time') < pd.to_datetime('2023-05-17T12:0z').to_datetime64()
    )
)


df = table.to_pandas()
df=df[0::20].reset_index(drop=True)
df.time=pd.to_datetime(df.time)
df=df.sort_values('time').reset_index(drop=True)
#%%
plt.figure()
plt.plot(df.time)
plt.figure()
plt.plot(df.time,df.afsTangentH_wgs84)
plt.gcf().autofmt_xdate()
plt.ylabel('Nominal Tangent height')
plt.xlabel('Time')
plt.ylim(80, 100)

# %%
df.columns.tolist()

# %%
