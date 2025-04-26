#%%
from mats_utils.rawdata.read_data import load_multi_parquet
import datetime as DT
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.plotting.plotCCD import plot_image
from matplotlib import pyplot as plt
import numpy as np
from mats_utils.geolocation.altitude_correction import rows_to_altitudes
from mats_utils.geolocation.coordinates import col_heights
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
#%%
def create_altitudeprof_for_mid_of_image(CCDitem, common_heights = np.arange(60,110,0.25)*1000):
    # Function from Björns script, edited to not multiply be 10**12
    # This function averages some columns and
    # calculates tangent heights
    
    image = np.stack(CCDitem.ImageCalibrated)
    col = int(CCDitem['NCOL']/2) 
    #cs = col_heights(ch, col, 10, spline=True) # TEST BETTER 
    cs = col_heights(CCDitem, col, spline=True)
    heights = np.array(cs(range(CCDitem['NROW'])))

    profile = np.array(image[:, col-2:col+2].mean(axis=1))

    # set heights
    
    profile=np.interp(common_heights,heights,profile)
    return  profile

#%%


start = DT.datetime(2023, 3, 5, 0, 0, 0)
stop= DT.datetime(2023, 3, 6, 15, 0, 0)

#df=load_multi_parquet('/home/linda/MATS_data/CROPD_v0.6',start, stop)
filter_UV2={'CCDSEL': [6,6]} 
filter_IR1={'CCDSEL': [1,1], 'TPlat': [0, 20]}
filter_IR3={'channel': 'IR3'}

#df=load_multi_parquet('/home/linda/MATS_data/CROPD_v0.6',start, stop, filt=filter_IR1)

df=read_MATS_data(start, stop, level='1b', version='0.6', filter=filter_IR3)


#%%
# ALTITUDES
#fixaltitudes=np.array([80,81,82,83,84,85,86,87,88,89,90])*1000
#df['image_fixalt']= df.apply(rows_to_altitudes, axis=1, fixaltvec=fixaltitudes)

#2d height conversion
#df['IMAGE']=df['ImageCalibrated'] # hack to not have to send parameters to the function rows_to_altitudes
#df['image_fixalt']= df.apply(rows_to_altitudes, axis=1)

#altvec, profile=prepare_profile(df.iloc[0])
common_heights = np.arange(60,110,0.25)*1000
df['image_fixalt']= df.apply(create_altitudeprof_for_mid_of_image, axis=1)

#%%
#Use the function col_heights to add the lowest altidude point of the image for a certain column
icolmid=int(df['ImageCalibrated'].loc[0].shape[1]/2)
icolleft=0
icolright=int(df['ImageCalibrated'].loc[0].shape[1]-1)
df['lowpoint']=df.apply(lambda x: col_heights(x, icolmid)[0], axis=1)
df['roll_in_m']=df.apply(lambda x: col_heights(x, icolright)[0]-col_heights(x, icolleft)[0], axis=1)

#ccddf.loc[:,'image_fixalt']= ccddf.apply(rows_to_altitudes, args=(fixaltvec, 'ImageCalibrated'), axis=1)
    




#%%
#Plot the mean profile as function of afsTangentH_wgs84

#pointing_heights = np.array(df['afsTangentH_wgs84'])
#pointing_heights=np.array(df['TPheight'])

#dfm=df[df['lowpoint']>63000]
df['TPlocaltime_str'] = pd.to_datetime(df['TPlocaltime'],utc=True)
dfm=df[df['TPlocaltime_str'].dt.time>DT.time(12, 0, 0)]

dfm=df
 #               CCDs = CCDs[(CCDs['TPlocaltime'].dt.time > midday)]

#%%
dataadir='/Users/lindamegner/MATS/MATS-retrieval/data'
localfolder='pointing_heights'
load_data=True
if load_data:
    df=pd.read_pickle(dataadir+'/'+localfolder+'/df_IR1_0.6_5-6march.pkl')


#%%
pointing_heights=np.array(dfm['lowpoint'])

ialtitudes=[50, 100, 150]
colr=['r', 'g', 'b']
plt.figure()
for i, ialt in enumerate(ialtitudes):

    plt.plot(pointing_heights, np.array(dfm['image_fixalt'].apply(lambda x: x[ialt])), '.', alpha=0.5, label='Signal at '+str(common_heights[ialt]/1000)+' km', color=colr[i])
    plt.xlabel('Pointing Height')
    plt.ylabel('Mean image intensity at '+str(common_heights[ialt]/1000)+' km')
    plt.ylim(0, 230)
    #fit line to data
    m, b = np.polyfit(pointing_heights, np.array(dfm['image_fixalt'].apply(lambda x: x[ialt])), 1)
    plt.plot(pointing_heights, m*pointing_heights + b, color=colr[i], label='Fit at '+str(common_heights[ialt]/1000)+' km , slope:'+str(round(m, 3)))
plt.legend()
plt.title(str(start)+' to '+str(stop)+' channel:'+df['channel'].iloc[0])

# %%
rolling_in_m=np.array(dfm['roll_in_m'])
plt.figure()
for i, ialt in enumerate(ialtitudes):
    
        plt.plot(rolling_in_m, np.array(dfm['image_fixalt'].apply(lambda x: x[ialt])), '.', alpha=0.5, label='Signal at '+str(common_heights[ialt]/1000)+' km', color=colr[i])
        plt.xlabel('Rolling (diff between low left and low right) [m]')  
        plt.ylabel('Mean image intensity at '+str(common_heights[ialt]/1000)+' km')
        plt.ylim(0, 230)
        #fit line to data
        m, b = np.polyfit(rolling_in_m, np.array(dfm['image_fixalt'].apply(lambda x: x[ialt])), 1)
        plt.plot(rolling_in_m, m*rolling_in_m + b, color=colr[i], label='Fit at '+str(common_heights[ialt]/1000)+' km , slope:'+str(round(m, 3)))   
plt.legend()
plt.title(str(start)+' to '+str(stop)+' channel:'+df['channel'].iloc[0])

# %%
#save df as pickle
df.to_pickle(dataadir+'/'+localfolder+'/df_IR3_0.6_5-6march_all_lat.pkl')
             
# %%
#check for how many rows in the df the schdule_name is set to 'CROPD'
df['schedule_name'].value_counts()

# %%
# plot dfm['roll_in_m'] vs latitude (TPlat)

plt.figure()
TPlat=np.array(dfm['TPlat'])
plt.plot(TPlat, rolling_in_m, '.')
plt.xlabel('Latitude')
plt.ylabel('Rolling (diff between low left and low right) [m]')

# %%
dfm['CalibImage'][0].shape()
# %%
