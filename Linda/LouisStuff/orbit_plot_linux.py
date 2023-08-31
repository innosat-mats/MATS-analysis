#%%
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import argparse
from datetime import date, timedelta
from mats_utils.geolocation.coordinates import NADIR_geolocation
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
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Louis/geolocation')


from temp_nadir import *

#%%
##### TO BE MODIFIED #####
# directory where the geolocation of each pixel is stored
orb_dir = "/Users/lindamegner/MATS/MATS-retrieval/MATS-Data/NADIR_geolocation/geolocation_test_final"
# directory where the maps are stored
map_dir = "/Users/lindamegner/MATS/MATS-retrieval/Polar_plot_final"

# start and end of timerange
start_time = DT.datetime(2023, 4, 6, 0, 6, 0)
stop_time = DT.datetime(2023, 4, 6, 0, 8, 0)

nb_core = 1 # number of CPU cores to use in multiprocessing (~ half of the cores available)

# simple_weights = np.array(simple_mask,dtype=float)
simple_weights = np.ones((14,56))

# superposition condition. How to handle overlapping images when several orbits cross the same spot
# only usefull when plotting several orbits at the same time
superposition = 'latest' # latest image on top
# superposition = 'oldest' # oldest image on top
# superposition = 'highest' # pixel with highest value on top

##########################

#%%
# creating usefull directories
if not os.path.exists(orb_dir):
        os.makedirs(orb_dir, exist_ok=True)
if not os.path.exists(map_dir):
        os.makedirs(map_dir, exist_ok=True)

#%% Defining usefull function that parallelises geolocation (several CPU cores)
def parallel_geolocating(part,ccditems,temp_dir,geolocation=True):

    global lat_points, lon_points, sza_points, im_points

    n = len(ccditems)

    if ccditems.iloc[0]['DataLevel'] == 'L1B' :
        im_key = 'ImageCalibrated'
    elif ccditems.iloc[0]['DataLevel'] == 'L1A' :
        im_key = 'IMAGE'
    else :
        warnings.warn('DataLevel not recognized (should be L1A or L1B)')

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

    a,b = np.shape(ccditems.iloc[0][im_key])

    # defining coordinates arrays
    lat_points = np.zeros((steps,a,b))
    lon_points = np.zeros((steps,a,b))
    sza_points = np.zeros((steps,a,b))
    im_points = np.zeros((steps,a,b))

    k = 0

    # geolocating images
    for i in tqdm(range(start_point, end_point)):
        ccditem = ccditems.iloc[i]
        im = ccditem[im_key]
        try:
            if geolocation:
                lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_sample=6,y_sample=6,interp_method='quintic')
                lat_points[k, :, :] = lat_map
                lon_points[k, :, :] = lon_map
                sza_points[k, :, :] = sza_map

                np.save(arr=lat_points,file=f'{temp_dir}/lat_points{part}')
                np.save(arr=lon_points,file=f'{temp_dir}/lon_points{part}')
                np.save(arr=sza_points,file=f'{temp_dir}/sza_points{part}')

            im_points[k, :, :] = im

            np.save(arr=im_points,file=f'{temp_dir}/im_points{part}')

        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            continue

        k = k + 1


#%%
# loading data

print("\n\n\n #####################################################")
print(f"\n Loading images from {start_time} to {stop_time}")

# filter selecting Nadir chanel
filter={'CCDSEL':7}

# reading mesurements
df_tot= read_MATS_data(start_time, stop_time,filter,level='1b',version='0.5').sort_values('EXPDate')
df_tot.sort_values('EXPDate')
df = df_tot[:]

if df.iloc[0]['DataLevel'] == 'L1B' :
    im_key = 'ImageCalibrated'
    print('DataLevel L1b')
elif df.iloc[0]['DataLevel'] == 'L1A' :
    im_key = 'IMAGE'
    print('DataLevel L1a')
else :
    warnings.warn('DataLevel not recognized (should be L1A or L1B)')



#%%
# calculating start and end of each orbit
# slice orbite by orbit

orb_ind = [] # orb_ind[j][0] contains the index of the first image from the orbit nb j, orb_ind[j][1] from the last image
orb_times = [] # orb_times[j][0] contains the EXPDate of the first image from the orbit nb j, orb_times[j][1] from the last image
orb_start = df.iloc[0]['EXPDate']
orb_end = df.iloc[0]['EXPDate']
ind_start = 0
ind_end = 0
for i in range(len(df)-1):
    if df.iloc[i+1]['EXPDate']-df.iloc[i]['EXPDate'] > timedelta(seconds = 100):
        orb_end = df.iloc[i]['EXPDate']
        ind_end = i
        orb_times.append([orb_start,orb_end])
        orb_ind.append([ind_start,ind_end])
        orb_start = df.iloc[i+1]['EXPDate']
        ind_start = i+1
orb_times.append([orb_start,df.iloc[-1]['EXPDate']])
orb_ind.append([ind_start,len(df)-1])

#%%
# geolocating each image

print(f"Geolocating images")
geolocation = True
for i in range(len(orb_times)):
    orb_start = orb_times[i][0]
    orb_end = orb_times[i][1]
    ind_start = orb_ind[i][0]
    ind_end = orb_ind[i][1]
    #df = df[(pd.to_datetime(df['EXPDate'])>orb_start) & (pd.to_datetime(df['EXPDate'])<orb_end)]
    df_orb = df[ind_start:ind_end]

    print(f"\n Geolocating images from {orb_start} to {orb_end}")

    # geolocating images
    n = int(len(df_orb)) # total number
    a,b = np.shape(df_orb.iloc[0][im_key])

    # parallel processing
    temp_dir=f"{orb_dir}/{start_time.strftime('%Y_%m_%d')}_orb{i}"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)
    files_per_part = 100
    sets = int(np.ceil(len(df_orb)/files_per_part))
    args = []
    for i in range(sets):
        args.append([i,df_orb,temp_dir,geolocation])
    pool = multiprocessing.Pool(nb_core)
    pool.starmap(parallel_geolocating,args)

#%%

# defining the new coordinate grids for the stacked image
lonmin,lonmax = 0,+360
latmin, latmax = -90,+0

nb_lat = ceil((latmax-latmin)/0.05) # number of latitude steps, by default we use a latitude resolution of .05 deg (~ same as NADIR images)
nb_lon = ceil((lonmax-lonmin)/0.05) # number of longitude steps, , by default we use a longitude resolution of .05 deg (~ same as NADIR images)

lo = np.linspace(lonmin,lonmax,nb_lon)
la = np.linspace(latmin,latmax,nb_lat)
stacked_lon,stacked_lat = np.meshgrid(lo,la)

stacked_im_tot = np.full((nb_lat-1,nb_lon-1),np.nan)

#%%
# aggregate each image in one orbit to create an image mosaic (1 image by orbit)
for j in range(0,len(orb_ind)):

    print('#################')
    print(f"Loading and stacking orbit number {j+1}/{len(orb_ind)} ({orb_times[j][0]} to {orb_times[j][1]})")

    temp_dir=f"{orb_dir}/{start_time.strftime('%Y_%m_%d')}_orb{j}"

    weight_im = simple_weights

    for i in range(0,sets):
        try:
            im=np.load(f'{temp_dir}/im_points{i}.npy')
            lat=np.load(f'{temp_dir}/lat_points{i}.npy')
            lon=np.load(f'{temp_dir}/lon_points{i}.npy')
            sza=np.load(f'{temp_dir}/sza_points{i}.npy')
            n = np.shape(im)[0]

            weight_points = weight_im
            for k in range(n-1):
              weight_points = np.append(weight_points,weight_im,axis=0)

            if i == 0:
                im_points=im
                lat_points=lat
                lon_points=lon
                sza_points=sza
                weights=weight_points

            else:
                im_points=np.append(im_points,im,axis=0)
                lat_points=np.append(lat_points,lat,axis=0)
                lon_points=np.append(lon_points,lon,axis=0)
                sza_points=np.append(sza_points,sza,axis=0)
                weights=np.append(weights,weight_points,axis=0)
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            continue

        lon_points = lon_points%360
        #im_points = sza_points

    print('stacking the images on 1 orbit path')
    # stacking images
    stacked_im = average_stacking(im_points,lat_points,lon_points,la,lo,data_weights=weights)
    np.save(arr=stacked_im,file=f'{temp_dir}/stacked_orbit{j+1}')


#%% plotting the images on the defined projection

# Defining several projections
data_crs = ccrs.PlateCarree() # projection for the geolocation data (latitude/longitude)
latlon_projection = ccrs.PlateCarree() # lat/lon projection
ortho_projectionSP = ccrs.Orthographic(central_latitude=-90,central_longitude=0) # Over the South Pole
ortho_projectionNP = ccrs.Orthographic(central_latitude=+90,central_longitude=0) # Over the North Pole

# plotting each orbit
for j in range(0,len(orb_ind)):

    print('#################')
    print(f"Projecting orbit number {j+1}/{len(orb_ind)} ({orb_times[j][0]} to {orb_times[j][1]})")

    temp_dir=f"{orb_dir}/{start_time.strftime('%Y_%m_%d')}_orb{j}"

    stacked_im=np.load(f'{temp_dir}/stacked_orbit{j+1}.npy')

    orb_start = orb_times[j][0]
    orb_end = orb_times[j][1]
    ind_start = orb_ind[j][0]
    ind_end = orb_ind[j][1]

    # projecting stacked image on an orthographic projection (slow)
    # over the north pole
    plt.close('Stacked Orthographic corrected')
    plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
    plt.title(f"Orbit between {orb_start.strftime('%Y-%m-%d %H:%M:%S')} and {orb_end.strftime('%Y-%m-%d %H:%M:%S')}")
    ax = plt.axes(projection=ortho_projectionNP)
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    # c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
    c = ax.pcolormesh(stacked_lon,stacked_lat,stacked_im_tot, transform=data_crs)
    ax.text(df.iloc[ind_start]['satlon'],df.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df.iloc[ind_end]['satlon'],df.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)
    plt.colorbar(c,ax=ax)
    plt.savefig(f"{map_dir}/nadir_{orb_start.strftime('%Y_%m_%d_%H_%M_%S')}_{orb_end.strftime('%Y_%m_%d_%H_%M_%S')}_NP.png", format='png',dpi=250)
    plt.show()

    # plotting over the south pole
    plt.close('Stacked Orthographic corrected')
    plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
    plt.title(f"{len(orb_times)} orbits between {orb_start.strftime('%Y-%m-%d %H:%M:%S')} and {orb_end.strftime('%Y-%m-%d %H:%M:%S')}")
    ax = plt.axes(projection=ortho_projectionSP)
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    # c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
    c = ax.pcolormesh(stacked_lon,stacked_lat,stacked_im_tot, transform=data_crs)
    ax.text(df.iloc[ind_start]['satlon'],df.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df.iloc[ind_end]['satlon'],df.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)
    plt.colorbar(c,ax=ax)
    plt.savefig(f"{map_dir}/nadir_{orb_start.strftime('%Y_%m_%d_%H_%M_%S')}_{orb_end.strftime('%Y_%m_%d_%H_%M_%S')}_SP.png", format='png',dpi=250)
    plt.show()




#%%
# Plotting the orbits on one single image


for j in range(0,len(orb_ind)):

    print('#################')
    print(f"Loading orbit number {j+1}/{len(orb_ind)} ({orb_times[j][0]} to {orb_times[j][1]})")

    temp_dir=f"{orb_dir}/{start_time.strftime('%Y_%m_%d')}_orb{j}"

    stacked_im=np.load(f'{temp_dir}/stacked_orbit{j+1}.npy')
    if j == 0:
      stacked_im_tot = stacked_im

    print('adding the orbit path to the global image')
    print(f'Superposition parameter : {superposition}')
    if superposition == 'latest':
      stacked_im_tot = np.where(~np.isnan(stacked_im) & np.isnan(stacked_im_tot),stacked_im,stacked_im_tot)
      stacked_im_tot = np.where(~np.isnan(stacked_im) & ~np.isnan(stacked_im_tot),stacked_im,stacked_im_tot)
    elif superposition == 'oldest':
      stacked_im_tot = np.where(~np.isnan(stacked_im) & np.isnan(stacked_im_tot),stacked_im,stacked_im_tot)
    elif superposition == 'highest':
      stacked_im_tot = np.where(~np.isnan(stacked_im) & np.isnan(stacked_im_tot),stacked_im,stacked_im_tot)
      stacked_im_tot = np.where(~np.isnan(stacked_im) & ~np.isnan(stacked_im_tot) & (stacked_im>stacked_im_tot),stacked_im,stacked_im_tot)
    else:
      print("Superposition parameter not defined.")


print("Projecting several orbits")

# projecting over the north pole (Orthographic projection)
plt.close('Stacked Orthographic corrected')
plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
plt.title(f"{len(orb_times)} orbits between {orb_times[0][0].strftime('%Y-%m-%d %H:%M:%S')} and {orb_times[-1][1].strftime('%Y-%m-%d %H:%M:%S')} | Superposition = {superposition}")
ax = plt.axes(projection=ortho_projectionNP)
ax.set_global()
ax.coastlines()
ax.gridlines()
# c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
c = ax.pcolormesh(stacked_lon,stacked_lat,stacked_im_tot, transform=data_crs)
for i in range(len(orb_times)):
    orb_start = orb_times[i][0]
    orb_end = orb_times[i][1]
    ind_start = orb_ind[i][0]
    ind_end = orb_ind[i][1]
    ax.text(df.iloc[ind_start]['satlon'],df.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df.iloc[ind_end]['satlon'],df.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)
plt.colorbar(c,ax=ax)
plt.savefig(f"{map_dir}/nadir_superposition_{superposition}_{orb_times[0][0].strftime('%Y_%m_%d_%H_%M_%S')}_{orb_times[-1][0].strftime('%Y_%m_%d_%H_%M_%S')}_NP.png", format='png',dpi=250)
plt.show()

# projecting over the south pole (Orthographic projection)
plt.close('Stacked Orthographic corrected')
plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
plt.title(f"{len(orb_times)} orbits between {orb_times[0][0].strftime('%Y-%m-%d %H:%M:%S')} and {orb_times[-1][1].strftime('%Y-%m-%d %H:%M:%S')} | Superposition = {superposition}")
ax = plt.axes(projection=ortho_projectionSP)
ax.set_global()
ax.coastlines()
ax.gridlines()
# c = ax.pcolorfast(stacked_lon,stacked_lat,stacked_im, transform=data_crs,cmap='hsv')
c = ax.pcolormesh(stacked_lon,stacked_lat,stacked_im_tot, transform=data_crs)
for i in range(len(orb_times)):
    orb_start = orb_times[i][0]
    orb_end = orb_times[i][1]
    ind_start = orb_ind[i][0]
    ind_end = orb_ind[i][1]
    ax.text(df.iloc[ind_start]['satlon'],df.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df.iloc[ind_end]['satlon'],df.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)
plt.colorbar(c,ax=ax)
plt.savefig(f"{map_dir}/nadir_superposition_{superposition}_{orb_times[0][0].strftime('%Y_%m_%d_%H_%M_%S')}_{orb_times[-1][0].strftime('%Y_%m_%d_%H_%M_%S')}_SP.png", format='png',dpi=250)
plt.show()




# %%
