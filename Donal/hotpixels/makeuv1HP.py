# %%
import numpy as np
from numba import jit
import pandas as pd
import datetime as DT
from mats_utils.rawdata.read_data import read_MATS_data
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
import os
import pickle
import sqlite3 as sqlite
from scipy.signal import savgol_filter
from glob import glob
from mats_l1_processing.L1_calibration_functions import *
from mats_utils.geolocation.coordinates import *
import matplotlib.animation as animation
%matplotlib widget
#%%
starttime=DT.datetime(2023,2,11,0,0) 
stoptime=DT.datetime(2023,2,12,0,0)
filter_UV1={'CCDSEL': [5,5]}
#uv1=read_MATS_data(starttime,stoptime,level="1a",version="0.7",filter=filter_UV1)

# %%
#uv1.to_pickle('/Users/donal/projekt/SIW/uv1.pkl')
#uv1 = pd.read_pickle('/Users/donal/projekt/SIW/uv1.pkl')
uv1=pd.read_pickle('/Volumes/Elements/MATSpickles/uv12023-02-11.pkl')
uv1=pd.read_pickle('/Volumes/Elements/MATSpickles/ir42023-02-20.pkl')
ch=uv1

# %%
homecat = os.environ.get('HOME')
#os.remove(homecat + '/Downloads/hpms/uv1.db')
db = sqlite.connect(homecat + '/Downloads/hpms/ir3.db')
cur = db.cursor()
cur.execute("create table if not exists hotpixelmaps ('key' INTEGER not null,'datetime' REAL ,'channel' TEXT,'HPM' blob, PRIMARY KEY (key,channel))")
insertstr = "insert or replace into hotpixelmaps values (?, ?, ?, ? )"
#%%

ch=uv1[uv1.TPsza > 97]
i=100
plt.figure()
image=ch.IMAGE.iloc[i-1]    

plt.plot(ch.IMAGE.iloc[i-1][:,8])
plt.plot(savgol_filter(ch.IMAGE.iloc[i-1][:,8], 11, 3))
# %%
#for each column in the image, calcuate 

# Calculate the difference between the original and filtered values for each column
difference_values = []
#%%
def process_block(ch,threshold=10):
    def calculate_difference_values(image,hp):
        difference_values = np.zeros_like(image, dtype=np.int16)
        for col in range(image.shape[1]):
            column_data = image[:, col]
            filtered_column = savgol_filter(column_data-hp[:,col], 11, 3)
            difference = column_data - filtered_column
            difference_values[:, col] = difference.astype(np.int16)
        return difference_values
    hp_images = []
    se_coordinates = []
    se_values= []
    se_threshold = 400

    for i in range(0,len(ch)):
        image = ch.IMAGE.iloc[i]
        hp = np.zeros_like(image, dtype=np.int16)
        #hp=hp_average_image
        difference_values = calculate_difference_values(image, hp)
        hp[difference_values > threshold] = difference_values[difference_values > threshold]
        difference_values = calculate_difference_values(image, hp)
        hp[difference_values > threshold] = difference_values[difference_values > threshold]
        difference_values = calculate_difference_values(image, hp)
        hp[difference_values > threshold] = difference_values[difference_values > threshold]
        se_coordinates.append([np.argwhere(hp > se_threshold)])
        se_values.append([hp[hp > se_threshold]])
        hp[hp > se_threshold] = 0
        hp_images.append(hp)
        
    presence_count = np.sum([hp > 0 for hp in hp_images], axis=0)
    threshold_count = int(0.6 * len(hp_images))
    present_in_60_percent_mask = presence_count > threshold_count

    # Calculate the average image
    hp_average_image = np.mean([hp for hp in hp_images], axis=0).astype(np.int16)
    hp_average_image[~present_in_60_percent_mask] = 0
    return hp_images, hp_average_image, se_coordinates, se_values
#%%
hpavgs=[]
nimages=len(ch)
nblocks=50
nperblock=int(nimages/nblocks)
for i in range(nblocks):
    print(i)
    hp_images, hp_average_image, se_coordinates, se_values = process_block(ch[i*nperblock:(i+1)*nperblock])
    hpavgs.append(hp_average_image)
    
    
    
    
    
    
    # for j in range(nperblock):
    #     datekey=ch.iloc[i*nperblock+j].EXPDate.year*100000000+ch.iloc[i*nperblock+j].EXPDate.month*1000000+ch.iloc[i*nperblock+j].EXPDate.day*10000+ch.iloc[i*nperblock+j].EXPDate.hour*100+ch.iloc[i*nperblock+j].EXPDate.minute
    #     pic=pickle.dumps(hp_images[j], protocol=pickle.HIGHEST_PROTOCOL)  
    #     print(datekey)              
    #     cur.execute(insertstr,(datekey,ch.iloc[i*nperblock+j].EXPDate.to_pydatetime(),'UV1',pic))
    


# %%len

fig, ax = plt.subplots()
im = ax.imshow(hpavgs[0], cmap='gray', aspect='auto',origin='lower')
im.set_clim(0, 400)
plt.colorbar(im)


def update(frame):
    #print(frame)
    im.set_array(hpavgs[frame])
    ax.set_title(f'Frame {frame}')
    return [im]

ani = animation.FuncAnimation(fig, update, frames=range(50), blit=True, repeat=False)
plt.show()
#%%
# Identify the cells in the series of hp_images that are present in more than 60% of the images
presence_count = np.sum([hp > 0 for hp in hpavgs], axis=0)
threshold_count = int(0.6 * len(hpavgs))
present_in_60_percent_mask = presence_count > threshold_count

# Calculate the average image
hp_average_image = np.mean([hp for hp in hpavgs], axis=0).astype(np.int16)
hp_average_image[~present_in_60_percent_mask] = 0


# Plot the average image
plt.figure()
plt.imshow(hp_average_image, cmap='gray', aspect='auto', origin='lower')
plt.title('Average High Pass Filtered Image')
plt.colorbar()
plt.clim(0, 400)
plt.show()
#%%
# Create a list of images with the high pass average subtracted
images_minus_hp_average = [ch.IMAGE.iloc[i]-hp_average_image for i in range(len(ch.IMAGE))]

fig, ax = plt.subplots()
im = ax.imshow(images_minus_hp_average[0], cmap='viridis', aspect='auto', origin='lower')
im.set_clim(0, 2000)
plt.colorbar(im)

def update(frame):
    im.set_array(images_minus_hp_average[frame])
    return [im]

ani = animation.FuncAnimation(fig, update, frames=range(len(images_minus_hp_average)), blit=True, repeat=False)
plt.show()
#%%pl
def calculate_difference_values(image,hp):
    difference_values = np.zeros_like(image, dtype=np.int16)
    for col in range(image.shape[1]):
        column_data = image[:, col]
        filtered_column = savgol_filter(column_data-hp[:,col], 11, 3)
        difference = column_data - filtered_column
        difference_values[:, col] = difference.astype(np.int16)
    return difference_values
hp=np.zeros_like(image,dtype=np.int16)
difference_values = calculate_difference_values(image,hp)
print(np.mean(difference_values),np.std(difference_values))

threshold = -300
hp[difference_values > threshold] = difference_values[difference_values > threshold]
plt.figure()
plt.hist(hp.flatten(), bins=np.arange (threshold, 300, 1))
plt.title('Histogram of Differences Between Original and Filtered Values')
plt.xlabel('Difference Value')
plt.ylabel('Frequency')
plt.show()
# %%
plt.figure()
plt.imshow(hp, cmap='gray', aspect='auto')
plt.title('High Pass Filtered Image')
plt.colorbar()
#plt.clim(0, 400)
# %%
i=1000
plt.figure()
image=ch.IMAGE.iloc[i-1]  
coldata = image[:,8]
plt.plot(coldata)
plt.plot(savgol_filter(coldata, 11, 3))
plt.plot(savgol_filter(coldata-hp, 11, 3))
# %%
plt.figure()
plt.plot(newcol)
# %%
def process_Dayfile(filename,chan):
    ch=pd.read_pickle(filename)
    ch=ch[ch.TPsza > 97]
    hpavgs=[]
    nimages=len(ch)
    nperblock=50
    nblocks=int(nimages/nperblock)
    if nblocks>0:
        for i in range(nblocks):
            print(i)
            hp_images, hp_average_image, se_coordinates, se_values = process_block(ch[i*nperblock:(i+1)*nperblock])
            hpavgs.append(hp_average_image)
        presence_count = np.sum([hp > 0 for hp in hpavgs], axis=0)
        threshold_count = int(0.6 * len(hpavgs))
        present_in_60_percent_mask = presence_count > threshold_count

        # Calculate the average image
        hp_average_image = np.mean([hp for hp in hpavgs], axis=0).astype(np.int16)
        hp_average_image[~present_in_60_percent_mask] = 0
        datekey=ch.iloc[0].EXPDate.year*100000000+ch.iloc[0].EXPDate.month*1000000+ch.iloc[0].EXPDate.day*10000+ch.iloc[0].EXPDate.hour*100+ch.iloc[0].EXPDate.minute
        pic=pickle.dumps(hp_average_image, protocol=pickle.HIGHEST_PROTOCOL)  
        print(datekey)              
        cur.execute(insertstr,(datekey,ch.iloc[0].EXPDate.to_pydatetime(),chan,pic))
        db.commit()
    else: print('No blocks')
#%%
chan='IR4'
files= glob(f'/Volumes/Elements/MATSpickles/{chan.lower()}*')
for f in files:
    print(f)
    try:
        process_Dayfile(f,chan)
    except Exception as e:
        print(f"An error occurred: {e}")
    
# %%/Volumes/Elements/MATSpickles/uv22023-04-29.pkl