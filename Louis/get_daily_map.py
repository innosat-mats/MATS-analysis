#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import argparse
from datetime import date, timedelta
from mats_utils.plotting.plotCCD import all_channels_plot
#from mats_utils.daily_preview.temp_nadirs import NADIR_geolocation, average_stacking
from temp_nadir import NADIR_geolocation, average_stacking
import numpy as np
import pandas as pd
import multiprocessing
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import ceil
import cartopy.crs as ccrs
import sys
import os
from multiprocessing import Manager

#%%
def parallel_geolocating(part):

    global lat_points, lon_points, sza_points, im_points

    if n > files_per_part:

        if part == 0:
            start_point = 0
        else:
            start_point = part*files_per_part-1

        if (part+1)*files_per_part < n:
            end_point = ((part+1)*files_per_part)
        else:
            end_point = n

        print(f'start: {start_point}; end: {end_point}')

    else:
        start_point=0
        end_point=n

    steps = end_point-start_point

    # defining coordinates arrays
    lat_points = np.zeros((steps,a,b))
    lon_points = np.zeros((steps,a,b))
    sza_points = np.zeros((steps,a,b))
    im_points = np.zeros((steps,a,b))
    
    k = 0
    
    # geolocating images
    for i in tqdm(range(start_point, end_point)):
        ccditem = df1a.iloc[i]
        im = ccditem['IMAGE']
        try:
            lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_sample=6,y_sample=6,interp_method='quintic')
            lat_points[k, :, :] = lat_map
            lon_points[k, :, :] = lon_map   
            sza_points[k, :, :] = sza_map
            im_points[k, :, :] = im 

            np.save(arr=lat_points,file=f'{temp_dir}/lat_points{part}')
            np.save(arr=lon_points,file=f'{temp_dir}/lon_points{part}')
            np.save(arr=sza_points,file=f'{temp_dir}/sza_points{part}')
            np.save(arr=im_points,file=f'{temp_dir}/im_points{part}')

        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            continue

        k = k + 1

# local variables
# value of a saturated pixel 
sat_val = 32880

#%%
for dayvar in range(1):
    # times for start and stop
    # start_time = DT.datetime(2023, 3, dayvar, 0, 0, 0)
    # stop_time = DT.datetime(2023, 3, dayvar+1, 0, 0, 0)

    start_time = DT.datetime(2023, 3, 31, 0, 0, 0)
    stop_time = DT.datetime(2023, 4, 1, 0, 0, 0)

    print("\n\n\n #####################################################")
    print(f"\n Loading images from {start_time} to {stop_time}")

    # filter selecting Nadir chanel
    filter={'CCDSEL': [7,7]}

    # reading mesurements
    df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
    pd.set_option('display.max_rows', 100)
    print(len(df1a_tot))
    df1a = df1a_tot[:]

    # preview of the mapped area
    plt.close('Orbit preview latlon')
    plt.figure('Orbit preview latlon')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    plt.scatter(df1a['satlon'],df1a['satlat'],marker='.')
    plt.show()

    print(f"\n Geolocating images from {start_time} to {stop_time}")

    # geolocating images
    n = int(len(df1a)) # total number 
    a,b = np.shape(df1a.iloc[0]['IMAGE'])

    # parallel processing
    temp_dir=f"/home/louis/MATS/MATS-Data/nadir_animation/{start_time.strftime('%Y_%m_%d')}"
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    files_per_part = 100
    sets = int(np.ceil(len(df1a)/files_per_part))
    parts = list(np.arange(0, sets))
    pool = multiprocessing.Pool(4)
    pool.map(parallel_geolocating, parts)

im_points=np.array(np.shape(df1a.iloc[0]['IMAGE']))


#%%

temp_dir = "/home/louis/MATS/MATS-Data/nadir_animation/2023_12_20" 
#im_points=np.array(np.shape(df1a.iloc[0]['IMAGE']))
sets = 39


for i in range(0,sets):
    try:
        im=np.load(f'{temp_dir}/im_points{i}.npy')
        lat=np.load(f'{temp_dir}/lat_points{i}.npy')
        lon=np.load(f'{temp_dir}/lon_points{i}.npy')
        sza=np.load(f'{temp_dir}/sza_points{i}.npy')

        if i == 0:
            im_points=im
            lat_points=lat
            lon_points=lon
            sza_points=sza
        else:
            im_points=np.append(im_points,im,axis=0)
            lat_points=np.append(lat_points,lat,axis=0)
            lon_points=np.append(lon_points,lon,axis=0)
            sza_points=np.append(sza_points,sza,axis=0)
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        continue

# defining the new coordinate grids for the stacked image
latmin,latmax = np.min(lat_points),np.max(lat_points)
lonmin,lonmax = np.min(lon_points),np.max(lon_points)
nb_lat = ceil((latmax-latmin)/0.05) # number of latitude steps, by default we use a latitude resolution of .05 deg (~ same as NADIR images)
nb_lon = ceil((lonmax-lonmin)/0.05) # number of longitude steps, , by default we use a longitude resolution of .05 deg (~ same as NADIR images)

print(latmin)
print(latmax)

lo = np.linspace(lonmin,lonmax,nb_lon)
la = np.linspace(latmin,latmax,nb_lat)
stacked_lon,stacked_lat = np.meshgrid(lo,la)

# stacking images (surprisingly fast)
stacked_im = average_stacking(im_points,lat_points,lon_points,la,lo)

# Defining several projections
data_crs = ccrs.PlateCarree() # projection for the geolocation data (latitude/longitude)
latlon_projection = ccrs.PlateCarree() # lat/lon projection
# Orthographic projection (point of perspective at infinity ie. an entire hemisphere is visible)
# ortho_projection = ccrs.Orthographic(central_latitude=90,central_longitude=0) # Over the North Pole
ortho_projection = ccrs.Orthographic(central_latitude=-90,central_longitude=0) # Over the South Pole
# Over the image taken in the middle of the path
# h,l = np.shape(stacked_lat)
# ortho_projection = ccrs.Orthographic(central_latitude=stacked_lat[h//2,l//2],central_longitude=stacked_lon[h//2,l//2])



#%% projecting a mosaic of images on a lat/lon projection (fastest)
#im_skip = 5 # gap between images (with a value of 5 there is a full overlap)
#plt.close('Mosaic latlon')
#plt.figure('Mosaic latlon')
#ax = plt.axes(projection=latlon_projection)
#ax.set_global()
#ax.coastlines()
#ax.gridlines()
#for i in tqdm(range(0,n,5)):
#    c = ax.pcolorfast(lon_points[i,:,:],lat_points[i,:,:],im_points[i,:,:], transform=data_crs)
#plt.colorbar(c,ax=ax)
#plt.show()


#%% projecting a mosaic of images on an orthographic projection (slow)
#im_skip = 5 # gap between images (with a value of 5 there is a full overlap)
#plt.close('Mosaic Orthographic')
#plt.figure('Mosaic Orthographic')
#ax = plt.axes(projection=ortho_projection)
#ax.set_global()
#ax.coastlines()
#ax.gridlines()
#or i in tqdm(range(0,n,5)):
#    c = ax.pcolorfast(lon_points[i,:,:],lat_points[i,:,:],im_points[i,:,:], transform=data_crs)
#plt.colorbar(c,ax=ax)
#plt.show()



#%% projecting stacked image on a lat/lon projection (fastest)
plt.close('Stacked latlon')
plt.figure('Stacked latlon')
ax = plt.axes(projection=latlon_projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv',figsize=[12.8,9.6],dpi=500)
plt.colorbar(c,ax=ax)
plt.savefig(f'{temp_dir}/test.png', format='png',dpi=1000)
plt.show()




#%% projecting stacked image on an orthographic projection (slow)
plt.close('Stacked Orthographic')
plt.figure('Stacked Othographic')
ax = plt.axes(projection=ortho_projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
#c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
c = ax.pcolor(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
plt.colorbar(c,ax=ax)
plt.savefig(f'{temp_dir}/test.png', format='png')
#plt.show()





# %%
