#%% Import modules
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from math import *
from mats_l1_processing.pointing import pix_deg
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize_scalar
from skyfield.positionlib import ICRF
from skyfield.api import wgs84 
from skyfield.units import Distance
from skyfield import api as sfapi
from scipy.interpolate import RegularGridInterpolator


#%% Functions
def pix_deg2(ccditem, xpixel, ypixel):
    """
    Function to get the x and y angle from a pixel relative to the center of the CCD
    WARNING : no images are flipped in this function
    
    Parameters
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
        
    yCCDpix = (nrskip + nrbin * (ypixel+0.5)) # y position of the pixel on the CCD (0.5 at the bottom, 510.5 on top)
    xCCDpix = (ncskip + ncbin * (xpixel+0.5)) # x position of the pixel on the CCD (0.5 on the left, 2047.5 on the right)
    
    xdeg = (180/pi)*np.arctan(d*(xCCDpix/2048-0.5)/f) # angular deviation along the x axis in degrees
    ydeg = (180/pi)*np.arctan(h*(yCCDpix/511-0.5)/f) # angular deviation along the y axis in degrees
    return xdeg, ydeg

def deg_map(ccditem):
    """
    Function to get the x and y angular deviation map for each pixel of the image. 
    The deviation is given in degrees relative to the center of the CCD
    WARNING : no images are flipped before calculating the angular deviation
    
    Parameters
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
    xmap,ymap = pix_deg2(ccditem, xpixel, ypixel)
    return xmap,ymap


def funheight(s, t, pos, FOV):
    """
    Function to get the distance between a point at position pos + s*FOV and the surface of the Geoid (wgs84 model),
     at time t.
    
    
    Parameters
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
    
    
    Parameters
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
    res = minimize_scalar(funheight, args=(t, pos, FOV), bracket=(3e5, 8e5))
    return res


def NADIR_geolocation(ccditem,x_step=2,y_step=2):
    """
    Function to get the latitude, longitude and solar zenith angle map for each pixel of the image.
    The values are calculated for some points and then interpolated for each pixel.
    WARNING : no images are flipped
    
    Parameters
    ----------
    ccditem : CCDitem
        measurement
    x_step : int
        step along the x-axis in the image between 2 sampled points used for interpolation. The default value is 2.
    y_step : int
        step along the y-axis in the image between 2 sampled points used for interpolation. The default value is 2.
            
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
    x_deg_map, y_deg_map = deg_map(ccditem) # creating angle deviation map for each pixel (degress)
    a,b = np.shape(im)

    metoOHB  = R.from_matrix([[0,0,-1],[0,-1,0],[-1,0,0]])
    ts=sfapi.load.timescale()
    t=ts.from_datetime(ccditem['EXPDate'].replace(tzinfo=sfapi.utc)) # exposure time  
    q=ccditem.afsAttitudeState
    quat=R.from_quat(np.roll(q,-1)) # quaternion of MATS attitude (for the LIMB imager) 
    pos=ccditem.afsGnssStateJ2000[0:3] # position of MATS
    
    xd = range(0,b,x_step) # sampled pixels on the x axis
    yd = range(0,a,y_step) # sampled pixels on the y axis
    LAT = np.zeros((len(yd),len(xd)))
    LON = np.zeros((len(yd),len(xd)))
    SZA = np.zeros((len(yd),len(xd)))

    # computing the latitude, longitude and solar zenith angles at the intersection of the line of sight and the earth surface
    # only the line of sights from some sampled pixels are computed
    for i in range(len(yd)):
        for j in range(len(xd)):
                x = xd[j]
                y = yd[i]
                # angular transformations
                # rotation from the line of sight of the LIMB imager to the line of sight of the NADIR pixel
                angle = R.from_euler('XYZ', [x_deg_map[y,x],-(90-23)+y_deg_map[y,x],0] , degrees=True).apply([1, 0, 0])
                FOV = quat.apply(metoOHB.apply(angle)) # attitude state for the line of sight of the NADIR pixel    
                # finding the distance between the point pos and the Geoid along the line of sight
                res = findsurface(t,pos,FOV)
                newp = pos + res.x * FOV 
                newp = ICRF(Distance(m=newp).au, t=t, center=399) # point at the intersection between the line of sight at the pixel and the Geoid surface
                LAT[i,j]=wgs84.subpoint(newp).latitude.degrees # latitude of the point
                LON[i,j]=wgs84.subpoint(newp).longitude.degrees # longitude of the point    

                # finding the solar zenith angle of the point
                planets = sfapi.load('de421.bsp')
                earth=planets['Earth']
                sun=planets['Sun']
                SZA[i,j]=90-((earth+wgs84.subpoint(newp)).at(t).observe(sun).apparent().altaz())[0].degrees
    
    # interpolating the results along all the pixels
    interp_lat = RegularGridInterpolator((yd,xd),LAT,method="quintic",bounds_error=False,fill_value=None) 
    interp_lon = RegularGridInterpolator((yd,xd),LON,method="quintic",bounds_error=False,fill_value=None)
    interp_sza = RegularGridInterpolator((yd,xd),SZA,method="quintic",bounds_error=False,fill_value=None)

    X_map,Y_map = np.meshgrid(range(b),range(a))
    lat_map = interp_lat((Y_map,X_map))
    lon_map = interp_lon((Y_map,X_map))
    sza_map = interp_sza((Y_map,X_map))

    return(lat_map,lon_map,sza_map)



#%% local variables
# value of a saturated pixel 
sat_val = 32880

# times for start and stop
start_time = DT.datetime(2023, 1, 12, 5, 0, 0)
stop_time = DT.datetime(2023, 1, 12, 7, 0, 0)

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}


#%% reading mesurements
df1a= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
print(len(df1a))

# displaying keys
pd.set_option('display.max_rows', 100)
df1a.dtypes


# %% 
ccditem = df1a.iloc[1]
ccditem['IMAGE'] = np.fliplr(ccditem['IMAGE'])
im = ccditem['IMAGE']
a,b = np.shape(im)
X = range(b)
Y = range(a)
xpixel, ypixel = np.meshgrid(X,Y)

#%%
XDEG2,YDEG2 = deg_map(ccditem)

# %%
lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_step=4,y_step=2)

# %% plotting

plt.figure()
plt.title('Latitude')
plt.imshow(lat_map)
plt.colorbar()

plt.figure()
plt.title('Longitude')
plt.imshow(lon_map)
plt.colorbar()

plt.figure()
plt.title('Solar Zenith Angle')
plt.imshow(sza_map)
plt.colorbar()

#%%
n = len(df1a)
a,b = np.shape(df1a.iloc[0]['IMAGE'])
lat_points = np.zeros((n,a,b))
lon_points = np.zeros((n,a,b))
sza_points = np.zeros((n,a,b))
im_points = np.zeros((n,a,b))


for i in range(n):
    ccditem = df1a.iloc[i]
    im = ccditem['IMAGE']
    lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_step=4,y_step=2)
    lat_points[i,:,:] = lat_map
    lon_points[i,:,:] = lon_map
    sza_points[i,:,:] = sza_map
    im_points[i,:,:] = im
    print(f"image n# {i}/{n}")

#%%

##%matplotlib qt
plt.figure('map')
ax = plt.axes(projection=ccrs.PlateCarree())
plt.figure('map')
plt.pcolormesh(np.reshape(lon_points,(a,b*n)), np.reshape(lat_points,(a,b*n)), np.reshape(im_points,(a,b*n)))
ax.coastlines()
plt.show()   




# %%
