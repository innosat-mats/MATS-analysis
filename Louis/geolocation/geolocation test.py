
#%% Plotting and modules

%matplotlib qt5 
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
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
import time
from tqdm import tqdm

#%% local variables
# value of a saturated pixel 
sat_val = 32880

# times for start and stop
start_time = DT.datetime(2023, 3, 20, 10, 59, 57)
stop_time = DT.datetime(2023, 3, 20, 10, 59, 58)
start_time = DT.datetime(2022, 12, 20, 20, 30, 0)
stop_time = DT.datetime(2022, 12, 20, 22, 0, 0)
start_time = DT.datetime(2023, 3, 13, 3, 27, 0)
stop_time = DT.datetime(2023, 3, 13, 3, 29, 0)
start_time = DT.datetime(2023, 1, 12, 0, 0, 0)
stop_time = DT.datetime(2023, 1, 13, 0, 0, 0)
# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}

#%% reading mesurements
df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
pd.set_option('display.max_rows', 100)
df1a_tot.dtypes
print(len(df1a_tot))


ccditem = df1a_tot.iloc[0]
df1a = df1a_tot[:20]
ccditem = df1a.iloc[0]
a,b = np.shape(ccditem['IMAGE'])

#%% geolocating images
df1a = df1a_tot[200:500]
ccditem = df1a.iloc[0]
a,b = np.shape(ccditem['IMAGE'])



# KWARGS = [{'ccditem':ccditem},
    #           {'ccditem':ccditem, 'x_sample':None, 'y_sample':None, 'interp_method':'cubic'}, 
    #           {'ccditem':ccditem, 'x_sample':None, 'y_sample':None, 'interp_method':'pchip'},
    #           {'ccditem':ccditem, 'x_sample':b, 'y_sample':a, 'interp_method':'cubic'},
    #           {'ccditem':ccditem, 'x_sample':4, 'y_sample':4, 'interp_method':'cubic'},
    #           {'ccditem':ccditem, 'x_sample':4, 'y_sample':4, 'interp_method':'pchip'},
    #           {'ccditem':ccditem, 'x_sample':10, 'y_sample':10, 'interp_method':'pchip'},
    #           {'ccditem':ccditem, 'x_sample':6, 'y_sample':6, 'interp_method':'quintic'},
    #           {'ccditem':ccditem, 'x_sample':2, 'y_sample':2, 'interp_method':'linear'},
    #           {'ccditem':ccditem, 'x_sample':100, 'y_sample':100, 'interp_method':'pchip'},
    #           {'ccditem':ccditem, 'x_sample':100, 'y_sample':100, 'interp_method':'quintic'}]

# KWARGS = [{'ccditem':ccditem, 'x_sample':2, 'y_sample':2, 'interp_method':'linear'},
#           {'ccditem':ccditem, 'x_sample':4, 'y_sample':4, 'interp_method':'cubic'},
#           {'ccditem':ccditem, 'x_sample':6, 'y_sample':6, 'interp_method':'quintic'},
#           {'ccditem':ccditem, 'x_sample':10, 'y_sample':10, 'interp_method':'quintic'}]

KWARGS = [{'x_sample':2, 'y_sample':2, 'interp_method':'linear'},
          {'x_sample':4, 'y_sample':4, 'interp_method':'cubic'},
          {'x_sample':6, 'y_sample':6, 'interp_method':'quintic'},
          {'x_sample':10, 'y_sample':10, 'interp_method':'quintic'},
          {'x_sample':100, 'y_sample':100, 'interp_method':'quintic'}]
m = len(KWARGS) # number of interpolation parameters
n = len(df1a) # number of images used to test interpolation
MAX_ERR_LAT = np.zeros((m,n))
MAX_ERR_LON = np.zeros((m,n))
MAX_ERR_SZA = np.zeros((m,n))




for j in tqdm(range(0,n)): 
    ccditem = df1a.iloc[j] 
    lat_map_ref,lon_map_ref,sza_map_ref = NADIR_geolocation(ccditem)
    for i in range(m):
        kwargs = KWARGS[i]        
        lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,**kwargs)
        MAX_ERR_LAT[i,j] = np.max(np.abs(lat_map-lat_map_ref))
        MAX_ERR_LON[i,j] = np.max(np.abs(lon_map-lon_map_ref))
        MAX_ERR_SZA[i,j] = np.max(np.abs(sza_map-sza_map_ref))        
                

#%%
lat_res = (np.max(lat_map_ref)-np.min(lat_map_ref))/b
lon_res = (np.max(lon_map_ref)-np.min(lon_map_ref))/b
sza_res = (np.max(sza_map_ref)-np.min(sza_map_ref))/b
for i in range(0,m):
        kwargs = KWARGS[i]        
        print('\n___________________________')
        print(f"\n x_sample: {kwargs['x_sample']}")
        print(f"\n y_sample: {kwargs['y_sample']}")
        print(f"\n interp_method: {kwargs['interp_method']}")
        print(f"\n latitude : max error {np.max(MAX_ERR_LAT[i,:])} deg ; relative error ~ {np.max(MAX_ERR_LAT[i,:])/lat_res} pixels")
        print(f"\n longitude : max error {np.max(MAX_ERR_LON[i,:])} deg ; relative error ~ {np.max(MAX_ERR_LON[i,:])/lat_res} pixels")
        print(f"\n solar zenith angle : max error {np.max(MAX_ERR_SZA[i,:])} deg ; relative error ~ {np.max(MAX_ERR_SZA[i,:])/lat_res} pixels")
 

#%%

def simple_resampler(data,lat,lon,new_lat,new_lon,method='nearest'):
    points = np.zeros((np.shape(lat.ravel())[0],2))
    points[:,0] = lon.ravel()
    points[:,1] = lat.ravel()
    new_data = griddata(points, data.ravel(), (new_lon, new_lat), method=method)
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



#%% geolocating images
df1a = df1a_tot[:10]
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
projection = ccrs.Orthographic(central_latitude=new_lat[0,0],central_longitude=new_lon[0,0])
projection = ccrs.Orthographic(central_latitude=50,central_longitude=new_lon[0,0])
plt.close('average polar')
plt.figure('average polar')
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()
ax.gridlines()
c = ax.pcolor(new_lon,new_lat,new_im3, transform=data_crs)
plt.colorbar(c,ax=ax)
plt.show()

#%%
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
# Get the l1a photometer data (uncalibrated data)
df = read_MATS_PM_data(start_time,stop_time)

#%%
# Plot temperatures, uncalibrated
plt.plot(df['PMTime'], df['PM1A']/df['PM1ACNTR'], label="PM1_Tphotodiode", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM1B']/df['PM1BCNTR'], label="PM1_Tifilter", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2A']/df['PM2ACNTR'], label="PM2_Tphotodiode", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2B']/df['PM2BCNTR'], label="PM2_Tifilter", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Temperature [bits]')
plt.show()

# Plot photometer signal, uncalibrated
plt.plot(df['PMTime'], df['PM1S']/df['PM1SCNTR'], label="PM1_Signal, Bkg Phot", marker='+', markersize=2.0, linewidth=0)
plt.plot(df['PMTime'], df['PM2S']/df['PM1SCNTR'], label="PM2_Signal, A-band Phot", marker='+', markersize=2.0, linewidth=0)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Signal [bits]')
plt.show()

# %%
# Plot photometer ratio vs satlat and satlon (higher value means higher emission altitude, 0.6 = surface; 13.5 = 10 km cloud top)
plt.figure()
plt.scatter(df['satlon'],df['satlat'],c=df['PM2S']/df['PM1S'], marker='.',clim=[0.6,1.4])
plt.colorbar()
# %%

