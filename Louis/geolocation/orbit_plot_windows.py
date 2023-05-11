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
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize_scalar
from skyfield.positionlib import ICRF
from skyfield.api import wgs84 
from skyfield.units import Distance
from skyfield import api as sfapi


#%%

# parameters to change

orb_dir = "/home/louis/MATS/MATS-Data/Geolocation_storage" # directory where the computed geolocation data is stored
map_dir = "/home/louis/MATS/MATS-Data/Polar_plot_test" # directory where the plots are stored

start_time = DT.datetime(2023, 4, 11, 13, 0, 0)
stop_time = DT.datetime(2023, 4, 11, 17, 0, 0)

nb_core = 4 # number of CPU cores to use in multiprocessing


#%%

def parallel_geolocating(part,ccditems,temp_dir):

    global lat_points, lon_points, sza_points, im_points

    n = len(ccditems)

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

    a,b = np.shape(ccditems.iloc[0]['IMAGE'])

    # defining coordinates arrays
    lat_points = np.zeros((steps,a,b))
    lon_points = np.zeros((steps,a,b))
    sza_points = np.zeros((steps,a,b))
    im_points = np.zeros((steps,a,b))
    
    k = 0
    
    

    # geolocating images
    for i in tqdm(range(start_point, end_point)):
        ccditem = ccditems.iloc[i]
        im = ccditem['IMAGE']
        try:
            lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_sample=6,y_sample=6,interp_method='quintic')
            lat_points[k, :, :] = lat_map
            lon_points[k, :, :] = lon_map   
            sza_points[k, :, :] = sza_map
            im_points[k, :, :] = im 

            np.save(arr=lat_points,file=f'{temp_dir}\\lat_points{part}')
            np.save(arr=lon_points,file=f'{temp_dir}\\lon_points{part}')
            np.save(arr=sza_points,file=f'{temp_dir}\\sza_points{part}')
            np.save(arr=im_points,file=f'{temp_dir}\\im_points{part}')

        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            continue

        k = k + 1





def average_stacking(data,lat,lon,stacked_lat,stacked_lon,no_holes=True):
    """
    Function stacking several geolocated images into one image

   
    Arguments
    ----------
    data : array[float]
        array containing the pixel values for each image to be stacked
    lat : array[float]
        array containing the latitude coordinates for each pixel to be stacked
    lon : array[float]
        array containing the longitude coordinates for each pixel to be stacked
    stacked_lat : array[float]
        array containing the latitude coordinates for the pixels in the generated image. It has to be an array of evenly distributed values.
    stacked_lon : array[float]
        array containing the longitude coordinates for the pixels in the generated image. It has to be an array of evenly distributed values. 
    --> there is no required shape for the parameters, but data, lat and lon must have the same shape, aswell as stacked_lat and stacked_lon
    
    
    Returns
    -------
    stacked_im : array[float]
        stacked image, has the same shape as stacked_lon and stacked_lat. If there was no corresponding pixels in the input images, the value of the pixel in stacked_im is None
           
    """
    lat_step = stacked_lat[1]-stacked_lat[0] # an assumed constant latitude step 
    lon_step = stacked_lon[1]-stacked_lon[0] # an assumed constant latitude step
    lon_bins = np.linspace(stacked_lon[0]-lon_step*0.5,stacked_lon[-1]+lon_step*0.5,len(stacked_lon)) # bins for longitude coordinates in the stacked image
    lat_bins = np.linspace(stacked_lat[0]-lat_step*0.5,stacked_lat[-1]+lat_step*0.5,len(stacked_lat)) # bins for latitude coordinates in the stacked image
    data_points = data.ravel()
    lon_ind = np.digitize(lon.ravel(),lon_bins)
    lat_ind = np.digitize(lat.ravel(),lat_bins)    
    stacked_im = np.zeros((len(lat_bins),len(lon_bins)))
    stacked_im_nb = np.zeros((len(lat_bins),len(lon_bins))) # number of pixels used from the input images for each pixel of the output image
    for i in tqdm(range(len(lon.ravel()))):
        stacked_im[lat_ind[i]-1,lon_ind[i]-1] += data_points[i]
        stacked_im_nb[lat_ind[i]-1,lon_ind[i]-1] += 1   
    

    for i in tqdm(range(len(lat_bins))):
        for j in range(len(lon_bins)):
            if stacked_im_nb[i,j] > 0:
                stacked_im[i,j] = stacked_im[i,j]/stacked_im_nb[i,j] # averaging each image 
            else:
                stacked_im[i,j] = None # if no pixel in the input images corresponds to the output coordinates
    
    stacked_im = stacked_im[:-1,:-1]
    stacked_im_nb = stacked_im_nb[:-1,:-1]
    stacked_lat = stacked_lat[:-1]
    stacked_lon = stacked_lon[:-1]
    grid_lon,grid_lat = np.meshgrid(stacked_lon,stacked_lat)
    stacked_coord = np.transpose([grid_lon.ravel(),grid_lat.ravel()])
    stacked_coord = stacked_coord[~np.isnan(stacked_im.ravel())]
    stacked_im_points = stacked_im.ravel()[~np.isnan(stacked_im.ravel())]


    if no_holes: 
        
        sat_path = stacked_im_nb > 0
        tmp = np.copy(sat_path)
        for i in range(1,len(lat_bins)-2):
            for j in range(1,len(lon_bins)-2):
                if not sat_path[i,j]:
                    if sat_path[i-1,j] or sat_path[i+1,j] or sat_path[i,j-1] or sat_path[i,j+1]:
                        tmp[i,j] = True
                        stacked_im[i,j] = np.nanmean(stacked_im[i-1:i+2,j-1:j+2])       
                    
                        
        sat_path = np.copy(tmp)      

        # rerun to remove the last holes  
        for i in range(1,len(lat_bins)-2):
            for j in range(1,len(lon_bins)-2):
                if not sat_path[i,j]:
                    if sat_path[i-1,j] or sat_path[i+1,j] or sat_path[i,j-1] or sat_path[i,j+1]:
                        tmp[i,j] = True
                        stacked_im[i,j] = np.nanmean(stacked_im[i-1:i+2,j-1:j+2]) 
        
        sat_path = np.copy(tmp)      
              
        # interp_stacked_im = LinearNDInterpolator(stacked_coord,stacked_im_points,fill_value=np.nan)
        # stacked_im_no_hole = interp_stacked_im((grid_lon,grid_lat))
        # stacked_im = stacked_im_no_hole * sat_path  
        stacked_im[stacked_im == 0] = np.nan      

    
    return (stacked_im) # removing last row and column which are binning artefacts







#%%


print("\n\n\n #####################################################")
print(f"\n Loading images from {start_time} to {stop_time}")

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}

# reading mesurements
df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5').sort_values('EXPDate')
df1a_tot.sort_values('EXPDate')
df1a = df1a_tot[:]
df1a = df1a[~np.isnan(df1a['satlat'])]


#%% 
# slice orbite by orbit

orb_ind = []
orb_times = []
orb_start = df1a.iloc[0]['EXPDate']
orb_end = df1a.iloc[0]['EXPDate']
ind_start = 0
ind_end = 0
for i in range(len(df1a)-1):
    if df1a.iloc[i+1]['EXPDate']-df1a.iloc[i]['EXPDate'] > timedelta(seconds = 100):
        orb_end = df1a.iloc[i]['EXPDate']
        ind_end = i
        orb_times.append([orb_start,orb_end])
        orb_ind.append([ind_start,ind_end])
        orb_start = df1a.iloc[i+1]['EXPDate']
        ind_start = i+1
orb_times.append([orb_start,df1a.iloc[-1]['EXPDate']])
orb_ind.append([ind_start,len(df1a)-1])

#%% 
# geolocating

print(f"Geolocating images")

for i in range(len(orb_times)):
    orb_start = orb_times[i][0]
    orb_end = orb_times[i][1]
    ind_start = orb_ind[i][0]
    ind_end = orb_ind[i][1]
    #df1a = df1a[(pd.to_datetime(df1a['EXPDate'])>orb_start) & (pd.to_datetime(df1a['EXPDate'])<orb_end)]
    df1a_orb = df1a[ind_start:ind_end]

    print(f"\n Geolocating images from {orb_start} to {orb_end}")

    # geolocating images.to_pydatetime()
    n = int(len(df1a_orb)) # total number 
    a,b = np.shape(df1a_orb.iloc[0]['IMAGE'])

    # parallel processing
    temp_dir=f"{orb_dir}\\{start_time.strftime('%Y_%m_%d')}_orb{i}"
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    files_per_part = 100
    sets = int(np.ceil(len(df1a_orb)/files_per_part))
    args = []
    for i in range(sets):
        args.append([i,df1a_orb,temp_dir])
    pool = multiprocessing.Pool(nb_core)
    pool.starmap(parallel_geolocating,args)

#%% 
# stack orbit by orbit

# defining the new coordinate grids for the stacked image
lonmin,lonmax = 0,+360
latmin, latmax = -90,-10

nb_lat = ceil((latmax-latmin)/0.05) # number of latitude steps, by default we use a latitude resolution of .05 deg (~ same as NADIR images)
nb_lon = ceil((lonmax-lonmin)/0.05) # number of longitude steps, , by default we use a longitude resolution of .05 deg (~ same as NADIR images)

lo = np.linspace(lonmin,lonmax,nb_lon)
la = np.linspace(latmin,latmax,nb_lat)
stacked_lon,stacked_lat = np.meshgrid(lo,la)

stacked_im_tot = np.full((nb_lat-1,nb_lon-1),np.nan)

for j in range(0,len(orb_ind)):
   
    print('#################')
    print(f"Loading and stacking orbit number {j+1}/{len(orb_ind)} ({orb_times[j][0]} to {orb_times[j][1]})")

    sets = 11

    temp_dir=f"{orb_dir}\\{start_time.strftime('%Y_%m_%d')}_orb{j}"
        
    

    for i in range(0,sets):
        try:
            im=np.load(f'{temp_dir}\\im_points{i}.npy')
            lat=np.load(f'{temp_dir}\\lat_points{i}.npy')
            lon=np.load(f'{temp_dir}\\lon_points{i}.npy')
            sza=np.load(f'{temp_dir}\\sza_points{i}.npy')

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
        
        lon_points = lon_points%360

     
     
    print('stacking the images on 1 orbit path')
    # stacking images (surprisingly fast)
    stacked_im = average_stacking(im_points,lat_points,lon_points,la,lo,no_holes = True)

    print('adding the orbit path to the global image')
    stacked_im_tot = np.where(~np.isnan(stacked_im),stacked_im,stacked_im_tot)
    
# stacked_lon = stacked_lon%360






#%%
# Defining several projections
data_crs = ccrs.PlateCarree() # projection for the geolocation data (latitude/longitude)
latlon_projection = ccrs.PlateCarree() # lat/lon projection
ortho_projectionSP = ccrs.Orthographic(central_latitude=-90,central_longitude=0) # Over the South Pole
ortho_projectionNP = ccrs.Orthographic(central_latitude=+90,central_longitude=0) # Over the North Pole





#%% projecting stacked image on an orthographic projection (slow)
plt.close('Stacked Orthographic corrected')
plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
plt.title(f"{len(orb_times)} orbits between {orb_times[0][0].strftime('%Y-%m-%d %H:%M:%S')} and {orb_times[-1][1].strftime('%Y-%m-%d %H:%M:%S')}")
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
    ax.text(df1a.iloc[ind_start]['satlon'],df1a.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df1a.iloc[ind_end]['satlon'],df1a.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)

plt.colorbar(c,ax=ax)
#ax.axis('off')
plt.savefig(f"{map_dir}\\nadir_{orb_times[0][0].strftime('%Y_%m_%d_%H_%M_%S')}_NP.png", format='png',dpi=250)
plt.show()


plt.close('Stacked Orthographic corrected')
plt.figure('Stacked Othographic corrected',figsize=(20,20),dpi=250)
plt.title(f"{len(orb_times)} orbits between {orb_times[0][0].strftime('%Y-%m-%d %H:%M:%S')} and {orb_times[-1][1].strftime('%Y-%m-%d %H:%M:%S')}")
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
    ax.text(df1a.iloc[ind_start]['satlon'],df1a.iloc[ind_start]['satlat'],f"{orb_start.strftime('%H:%M:%S')}",transform=data_crs)
    ax.text(df1a.iloc[ind_end]['satlon'],df1a.iloc[ind_end]['satlat'],f"{orb_end.strftime('%H:%M:%S')}",transform=data_crs)

plt.colorbar(c,ax=ax)
#ax.axis('off')
plt.savefig(f"{map_dir}\\nadir_{orb_times[0][0].strftime('%Y_%m_%d_%H_%M_%S')}_SP.png", format='png',dpi=250)
plt.show()




# %%
