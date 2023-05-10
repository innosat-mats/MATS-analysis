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


#%%

def pix_deg(ccditem, xpixel, ypixel):
    """
    Function to get the x and y angle from a pixel relative to the center of the CCD
        
    Arguments
    ----------
    ccditem : CCDitem
        measurement
    xpixel : int or array[int]
        x coordinate of the pixel(s) in the image
    ypixel : int or array[int]
        y coordinate of the pixel(s) in the image
        
    Returns
    -------
    xdeg : float or array[float]
        angular deviation along the x axis in degrees (relative to the center of the CCD)
    ydeg : float or array[float]
        angular deviation along the y axis in degrees (relative to the center of the CCD) 
    """
    h = 6.9 # height of the CCD in mm
    d = 27.6 # width of the CCD in mm
    # selecting effective focal length
    if (ccditem['CCDSEL']) == 7: # NADIR channel
        f = 50.6 # effective focal length in mm
    else: # LIMB channels
        f = 261    
    
    ncskip = ccditem['NCSKIP']
    try:
        ncbin = ccditem['NCBIN CCDColumns']
    except:
        ncbin = ccditem['NCBINCCDColumns']
    nrskip = ccditem['NRSKIP']
    nrbin = ccditem['NRBIN']
    ncol = ccditem['NCOL'] # number of columns in the image MINUS 1

    y_disp = (h/(f*511))
    x_disp = (d/(f*2048))
  
    if (ccditem['CCDSEL']) in [1, 3, 5, 6, 7]:
        xdeg = np.rad2deg(np.arctan(x_disp*((2048-ncskip - (ncol+1)*ncbin + ncbin*(xpixel+0.5)) - 2047./2)))
    else:
        xdeg = np.rad2deg(np.arctan(x_disp*(ncskip + ncbin * (xpixel+0.5) - 2047./2)))
        
    ydeg = np.rad2deg(np.arctan(y_disp*(nrskip + nrbin * (ypixel+0.5) - 510./2)))

    return xdeg, ydeg

def deg_map(ccditem):
    """
    Function to get the x and y angular deviation map for each pixel of the image. 
    The deviation is given in degrees relative to the center of the CCD
    
    
    Arguments
    ----------
    ccditem : CCDitem
        measurement
            
    Returns
    -------
    xmap : array[float]
        angular deviation map along the x axis in degrees (relative to the center of the CCD)
    ymap : array[float]
        angular deviation map along the y axis in degrees (relative to the center of the CCD) 
    """    
    im = ccditem['IMAGE']

    a,b = np.shape(im)
    X = range(b)
    Y = range(a)
    xpixel, ypixel = np.meshgrid(X,Y)
    xmap,ymap = pix_deg(ccditem, xpixel, ypixel)
    return xmap,ymap


def funheight_square(s, t, pos, FOV):
    """
    Function to get the distance between a point at position pos + s*FOV and the surface of the Geoid (wgs84 model),
     at time t.
    
    
    Arguments
    ----------
    s : float
        length along the straight line
    t : skyfield.timelib.Time
        time
    pos : array[float]
        position in space where the line starts (~position of MATS). Array of 3 position coordinates in m in the ICRF reference frame
    FOV : array[float]
        angle of the line (direction of the line), array of 3 elements in the IFRC reference frame
            
    Returns
    -------
    elevation**2 : float
        elevation squared of the point pos+s*FOV in m**2
    """
    newp = pos + s * FOV
    newp = ICRF(Distance(m=newp).au, t=t, center=399)
    return wgs84.subpoint(newp).elevation.m**2


def findsurface(t, pos, FOV):
    """
    Function to get the distance between a point at position pos and the surface of the Geoid (wgs84 model),
     at time t, along the line oriented along the FOV direction and starting at position pos
    
    
    Arguments
    ----------
    t : skyfield.timelib.Time
        time
    pos : array[float]
        position in space where the line starts (~position of MATS). Array of 3 position coordinates in m in the ICRF reference frame
    FOV : array[float]
        angle of the line (direction of the line), array of 3 elements in the IFRC reference frame
            
    Returns
    -------
    res : OptimizeResult object
        res.x is the distance found in m   
    """
    res = minimize_scalar(funheight_square, args=(t, pos, FOV), bracket=(3e5, 8e5))
    return res


def NADIR_geolocation(ccditem,x_sample=None,y_sample=None,interp_method='quintic'):
    """
    Function to get the latitude, longitude and solar zenith angle map for each pixel of the image.
    The values are calculated for some points and then interpolated for each pixel.
    WARNING : no images are flipped
    
    Arguments
    ----------
    ccditem : CCDitem
        measurement
    x_sample : int
        number of geolocated points along the x axis used for the interpolation. Default value is None, which means that there is no interpolation along the x-axis (each value is computed)
    y_step : int
        number of geolocated points along the y axis used for the interpolation. Default value is None, which means that there is no interpolation along the y-axis (each value is computed)
    interp_method :
        interpolation method : 'linear', 'nearest', 'slinear', 'cubic', 'quintic' and 'pchip'
        WARNING : choose the minimum x and y sampling according to the interpolation method
            
    Returns
    -------
    lat_map : array[float]
        map giving the latitude for each pixel in the image
    lon_map : array[float]
        map giving the longitude for each pixel in the image
    sza_map : array[float]
        map giving the solar zenith angle for each pixel in the image
    """
    im = ccditem['IMAGE']
    x_deg_map, y_deg_map = deg_map(ccditem) # creating angle deviation map for each pixel (degrees)
    a,b = np.shape(im)

    metoOHB  = R.from_matrix([[0,0,-1],[0,-1,0],[-1,0,0]])
    ts=sfapi.load.timescale()
    t=ts.from_datetime((ccditem['EXPDate']+timedelta(seconds=ccditem['TEXPMS']/(2*1000))).replace(tzinfo=sfapi.utc)) # exposure time (middle of the exposure timespan)  
    q=ccditem.afsAttitudeState
    quat=R.from_quat(np.roll(q,-1)) # quaternion of MATS attitude (for the satellite frame) 
    pos=ccditem.afsGnssStateJ2000[0:3] # position of MATS
    
    if x_sample == None or x_sample >= b: # no upsampling
        x_sample = b
    if y_sample == None or y_sample >= a: # no upsampling
        y_sample = a
    
    interpolation = True
    if x_sample == b and y_sample == a: # if both axis have enough sampling points, there is no interpolation
        interpolation = False

    xd = np.linspace(np.min(x_deg_map),np.max(x_deg_map),x_sample) # sampled angles on the x axis
    yd = np.linspace(np.min(y_deg_map),np.max(y_deg_map),y_sample) # sampled angles on the y axis
    x_deg_sample,y_deg_sample = np.meshgrid(xd,yd)

    if not interpolation:
        y_deg_sample,x_deg_sample = y_deg_map,x_deg_map # the sampled angles are the calculated angles for each pixel

    # sampled latitude, longitude and solar zenith angle values
    LAT = np.zeros((y_sample,x_sample))
    LON = np.zeros((y_sample,x_sample))
    SZA = np.zeros((y_sample,x_sample))

    # computing the latitude, longitude and solar zenith angles at the intersection of the line of sight and the earth surface
    # only the line of sights from some sampled pixels are computed
    for i in range(y_sample):
        for j in range(x_sample):
                # angular transformations
                # rotation from the line of sight of the LIMB imager to the line of sight of the NADIR pixel
                angle = R.from_euler('XYZ', [x_deg_sample[i,j],-(90-23)+y_deg_sample[i,j],0] , degrees=True).apply([1, 0, 0])
                FOV = quat.apply(metoOHB.apply(angle)) # attitude state for the line of sight of the NADIR pixel    
                # finding the distance between the point pos and the Geoid along the line of sight
                res = findsurface(t,pos,FOV)
                newp = pos + res.x * FOV 
                newp = ICRF(Distance(m=newp).au, t=t, center=399) # point at the intersection between the line of sight at the pixel and the Geoid surface
                LAT[i,j]=wgs84.subpoint(newp).latitude.degrees # latitude of the point
                LON[i,j]=wgs84.subpoint(newp).longitude.degrees # longitude of the point E [-180,+180] 

                # finding the solar zenith angle of the point
                planets = sfapi.load('de421.bsp')
                earth=planets['Earth']
                sun=planets['Sun']
                SZA[i,j]=90-((earth+wgs84.subpoint(newp)).at(t).observe(sun).apparent().altaz())[0].degrees
    
    # to get a continuous longitudinal field
    if np.max(LON)-np.min(LON) > 300: # this condition is met if points are on both sides of the -180/+180 deg line
        LON = np.where(LON<0,LON+360,LON)

    if interpolation: # interpolating the results along all the pixels
        # each interpolator object takes as argument an y and x angular deviation and gives a lat/lon/sza value
        interp_lat = RegularGridInterpolator((yd,xd),LAT,interp_method,bounds_error=False,fill_value=None) 
        interp_lon = RegularGridInterpolator((yd,xd),LON,interp_method,bounds_error=False,fill_value=None)
        interp_sza = RegularGridInterpolator((yd,xd),SZA,interp_method,bounds_error=False,fill_value=None)
        # interpolating on the real angular deviations for each pixel
        lat_map = interp_lat((y_deg_map,x_deg_map))
        lon_map = interp_lon((y_deg_map,x_deg_map))
        sza_map = interp_sza((y_deg_map,x_deg_map))
    else: # no interpolation       
        lat_map = LAT
        lon_map = LON
        sza_map = SZA

    return(lat_map,lon_map,sza_map)


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
    pool = multiprocessing.Pool(4)
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
plt.savefig(f"{map_dir}\\nadir_{orb_times[0][0].strftime('%Y_%m_%d')}_NP.png", format='png',dpi=250)
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
plt.savefig(f"{map_dir}\\nadir_{orb_times[0][0].strftime('%Y_%m_%d')}_SP.png", format='png',dpi=250)
plt.show()




# %%
