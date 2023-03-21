#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from matplotlib import pyplot as plt
from mats_utils.geolocation import coordinates
import numpy as np

#%% Select on time where we sample entire CCD
#start_time = DT.datetime(2023,1,9,18,0,0)
#stop_time = DT.datetime(2023,1,9,18,30,0)

#df = read_MATS_data(start_time,stop_time,version='0.5',level='1a')
#df.to_pickle('data')
df = pd.read_pickle('/home/olemar/Projects/Universitetet/MATS/MATS-analysis/OleM/data/desmear_profiles.pkl')
CCDitems = dataframe_to_ccd_items(df)
# %%
def plot_profile(df,n):
    col=int(df.iloc[n]['NCOL']/2)
    cs=coordinates.col_heights(df.iloc[n],col,10,spline=True)
    #%%
    plt.figure()
    plt.ylim([50,110])
    plt.title ('IR2 compared to Sasktran')
    plt.ylabel('Altitude km')
    plt.xlabel('Radiance Witts')
    plt.semilogx(np.vstack(df['ImageCalibrated'][n])[:,col-2:col+2].mean(axis=1),cs(np.arange(df['NROW'][n]))/1000)
# %%
plot_profile(df,1)

# %%
