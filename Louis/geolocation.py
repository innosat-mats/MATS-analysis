

#%% Plotting and modules

%matplotlib qt5 
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from math import *
from mats_l1_processing.pointing import pix_deg
from mats_utils.geolocation.coordinates import NADIR_geolocation
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize_scalar
from skyfield.positionlib import ICRF
from skyfield.api import wgs84 
from skyfield.units import Distance
from skyfield import api as sfapi
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from tqdm import tqdm

#%% local variables
# value of a saturated pixel 
sat_val = 32880

# times for start and stop
start_time = DT.datetime(2023, 1, 12, 3, 30, 0)
stop_time = DT.datetime(2023, 1, 12, 4, 0, 0)
start_time = DT.datetime(2023, 3, 13, 3, 27, 0)
stop_time = DT.datetime(2023, 3, 13, 3, 29, 0)
start_time = DT.datetime(2023, 3, 20, 12, 0, 0)
stop_time = DT.datetime(2023, 3, 20, 13, 0, 0)
start_time = DT.datetime(2022, 12, 20, 20, 30, 0)
stop_time = DT.datetime(2022, 12, 20, 22, 0, 0)
# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}

#%% reading mesurements
df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
print(len(df1a_tot))

# displaying keys
pd.set_option('display.max_rows', 100)
df1a_tot.dtypes
df1a = df1a_tot


#%% selecting interesting ccditems in df1a
# df1a = df1a_tot[df1a_tot['nadir_sza']>105.7]
df1a = df1a_tot[df1a_tot['satlon']>179]
df1a = df1a[df1a['satlon']<181]
print(len(df1a))
for i in range(len(df1a)):
    print(df1a.iloc[i]['satlon'])
df1a = df1a_tot[:50]


#%% geolocating images

n = len(df1a)
a,b = np.shape(df1a.iloc[0]['IMAGE'])
lat_points = np.zeros((n,a,b))
lon_points = np.zeros((n,a,b))
sza_points = np.zeros((n,a,b))
im_points = np.zeros((n,a,b))


for i in tqdm(range(n)):
    ccditem = df1a.iloc[i]
    im = ccditem['IMAGE']
    lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_sample=6,y_sample=6,interp_method='quintic')
    lat_points[i,:,:] = lat_map
    lon_points[i,:,:] = lon_map    
    sza_points[i,:,:] = sza_map
    im_points[i,:,:] = im    



#%%

def simple_resampler(data,lat,lon,new_lat,new_lon,method='nearest'):
    points = np.zeros((np.shape(lat.ravel())[0],2))
    points[:,0] = lon.ravel()
    points[:,1] = lat.ravel()
    new_data = griddata(points, data.ravel(), (new_lon, new_lat), method='nearest')
    return new_data

def average_stacking(data,lat,lon,new_lat,new_lon):
    lat_step = new_lat[1]-new_lat[0]
    lon_step = new_lon[1]-new_lon[0]
    lon_bins = np.linspace(new_lon[0]-lon_step*0.5,new_lon[-1]+lon_step*0.5,len(new_lon))
    lat_bins = np.linspace(new_lat[0]-lat_step*0.5,new_lat[-1]+lat_step*0.5,len(new_lat))
    data_points = data.ravel()
    lon_ind = np.digitize(lon.ravel(),lon_bins)
    lat_ind = np.digitize(lat.ravel(),lat_bins)
    new_im = np.zeros((len(lat_bins),len(lon_bins)))
    new_im_nb = np.zeros((len(lat_bins),len(lon_bins)))
    for i in tqdm(range(len(lon.ravel()))):
        new_im[lat_ind[i]-1,lon_ind[i]-1] += data_points[i]
        new_im_nb[lat_ind[i]-1,lon_ind[i]-1] += 1    
    for i in tqdm(range(len(lat_bins))):
        for j in range(len(lon_bins)):
            if new_im_nb[i,j] > 0:
                new_im[i,j] = new_im[i,j]/new_im_nb[i,j]
            else :
                new_im[i,j] = None

    # data_crs = ccrs.PlateCarree()
    # plt.figure('average nb')
    # #projection = ccrs.Robinson()
    # projection = ccrs.PlateCarree()
    # #projection = ccrs.Orthographic(central_latitude=-90,central_longitude=0)
    # ax = plt.axes(projection=projection)
    # ax.set_global()
    # ax.coastlines()
    # ax.gridlines()
    # c= ax.pcolor(new_lon,new_lat,new_im_nb, transform=data_crs)
    # plt.colorbar(c,ax=ax)
    # plt.show()    
    return (new_im[:-1,:-1])


#%%
latmin,latmax = np.min(lat_points),np.max(lat_points)
lonmin,lonmax = np.min(lon_points),np.max(lon_points)
nb_lat = ceil((latmax-latmin)/0.05)
nb_lon = ceil((lonmax-lonmin)/0.05)

lo = np.linspace(lonmin,lonmax,nb_lon)
la = np.linspace(latmin,latmax,nb_lat)
new_lon,new_lat = np.meshgrid(lo,la)


new_im3 = average_stacking(im_points,lat_points,lon_points,la,lo)


#%% lat/lon projection
data_crs = ccrs.PlateCarree()
projection = ccrs.PlateCarree()
plt.close('average')
plt.figure('average')
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
c = ax.pcolor(new_lon,new_lat,new_im3, transform=data_crs)
plt.colorbar(c,ax=ax)
plt.show()


#%% polar projection
data_crs = ccrs.PlateCarree()
projection = ccrs.Orthographic(central_latitude=80,central_longitude=-170)
plt.close('average polar')
plt.figure('average polar')
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
c = ax.pcolor(new_lon,new_lat,new_im3, transform=data_crs)
plt.colorbar(c,ax=ax)
plt.show()

plt.close('mosaic polar')
plt.figure('mosaic polar')
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
for i in tqdm(range(0,n,5)):
    c = ax.pcolor(lon_points[i,:,:],lat_points[i,:,:],im_points[i,:,:], transform=data_crs)
plt.colorbar(c,ax=ax)
plt.show() 
# %%
