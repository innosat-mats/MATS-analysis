#%%
#%matplotlib qt5
import warnings
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
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull,Delaunay
from datetime import timedelta
from tqdm import tqdm


#%% Importing several functions
# TEMPORARY to be used until branch merging in MATS-utility-functions and MATS-L1-processing





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
    # im = ccditem['IMAGE']
    if ccditem['DataLevel'] == 'L1B' :
        im_key = 'ImageCalibrated'
    elif ccditem['DataLevel'] == 'L1A' :
        im_key = 'IMAGE'
    else :
        warnings.warn('DataLevel not recognized (should be L1A or L1B)')
    im = ccditem[im_key]
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
    if ccditem['DataLevel'] == 'L1B' :
        im_key = 'ImageCalibrated'
    elif ccditem['DataLevel'] == 'L1A' :
        im_key = 'IMAGE'
    else :
        warnings.warn('DataLevel not recognized (should be L1A or L1B)')
    im = ccditem[im_key] # image
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








def average_stacking(data,lat,lon,stacked_lat,stacked_lon,data_weights=None):
    """
    Function stacking several geolocated images into one image. All the images are binned following the stacked_lat and stacked_lon coordinates.
    For each pixel in the new stacked grid, the value is the weighted average of the corresponding pixels in the input images. The weights are given by data_weights.
    The imaging path is then determined, and pixels inside the path for which no corresponding pixels have been found in the input images are interpolated linearly


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
    data_weights : array[float]
        array containing the weights for each pixels to be stacked. If None, all images are weighted equally. Default value is None

    --> there is no required shape for the parameters, but data, data_weights, lat and lon must have the same shape, aswell as stacked_lat and stacked_lon


    Returns
    -------
    mosaic_im : array[float]
        stacked image, has the same shape as stacked_lon and stacked_lat. If there was no corresponding pixels in the input images, the value of the pixel in stacked_im is None

    """
    lat_step = stacked_lat[1]-stacked_lat[0] # an assumed constant latitude step
    lon_step = stacked_lon[1]-stacked_lon[0] # an assumed constant latitude step
    lon_bins = np.linspace(stacked_lon[0]-lon_step*0.5,stacked_lon[-1]+lon_step*0.5,len(stacked_lon)) # bins for longitude coordinates in the stacked image
    lat_bins = np.linspace(stacked_lat[0]-lat_step*0.5,stacked_lat[-1]+lat_step*0.5,len(stacked_lat)) # bins for latitude coordinates in the stacked image
    data_points = data.ravel()

    n,a,b = np.shape(data)


    if type(data_weights) == type(None):
        data_weights = np.ones_like(data_points)
    weights = data_weights.ravel()
    lon_ind = np.digitize(lon.ravel(),lon_bins)
    lat_ind = np.digitize(lat.ravel(),lat_bins)
    stacked_im = np.zeros((len(lat_bins),len(lon_bins)))
    stacked_im_nb = np.zeros((len(lat_bins),len(lon_bins))) # number of pixels used from the input images for each pixel of the output image
    for i in tqdm(range(len(lon.ravel())),desc='stacking images part 1/2'):
        stacked_im[lat_ind[i]-1,lon_ind[i]-1] += data_points[i]*weights[i]
        stacked_im_nb[lat_ind[i]-1,lon_ind[i]-1] += weights[i]


    for i in tqdm(range(len(lat_bins)),desc='stacking images part 2/2'):
        for j in range(len(lon_bins)):
            if stacked_im_nb[i,j] > 0:
                stacked_im[i,j] = stacked_im[i,j]/stacked_im_nb[i,j] # averaging each image
            else:
                stacked_im[i,j] = None # if no pixel in the input images corresponds to the output coordinates

    stacked_im = stacked_im[:-1,:-1] # removing last row and column which are binning artefacts
    stacked_im_nb = stacked_im_nb[:-1,:-1] # removing last row and column which are binning artefacts
    stacked_lat = stacked_lat[:-1] # removing last row
    stacked_lon = stacked_lon[:-1] # removing last column
    grid_lon,grid_lat = np.meshgrid(stacked_lon,stacked_lat)
    mosaic_coord = np.transpose([grid_lon.ravel(),grid_lat.ravel()]) # list of coordinates in the mosaic image
    stacked_coord = mosaic_coord[~np.isnan(stacked_im.ravel())] # list of coordinates containing stacked data
    stacked_im_points = stacked_im.ravel()[~np.isnan(stacked_im.ravel())] # list of the stacked data value

    print('interpolating images part 1/2')
    interp_stacked_im = LinearNDInterpolator(stacked_coord,stacked_im_points,fill_value=np.nan)
    im_mosaic_points = interp_stacked_im((mosaic_coord[:,0],mosaic_coord[:,1])) # the stacked data is interpolated on the output image grid


    inside_coord = np.full(np.shape(mosaic_coord)[0], False, dtype=bool) # array with value True if the corresponding point in mosaic_coord is inside the satellite imaging path, False otherwise
    for i in tqdm(range(n),desc='interpolating images part 2/2'):
        corner_list = np.array([[lon[i,0,0],lat[i,0,0]],[lon[i,-1,0],lat[i,-1,0]],[lon[i,-1,-1],lat[i,-1,-1]],[lon[i,0,-1],lat[i,0,-1]]])
        if np.abs(max(corner_list[:,0])-min(corner_list[:,0])) > 180: # dealing with images crossing the +/- 180 deg in longitude
            corner_list[:,0] = np.where(corner_list[:,0]>180,corner_list[:,0]-360,corner_list[:,0])
            triangulation = Delaunay(corner_list)
            inside_coord_im = triangulation.find_simplex(np.transpose([np.where(mosaic_coord[:,0]>180,mosaic_coord[:,0]-360,mosaic_coord[:,0]),mosaic_coord[:,1]])) >= 0 # list of points inside one image
        else:
            triangulation = Delaunay(corner_list)
            inside_coord_im = triangulation.find_simplex(mosaic_coord) >= 0 # list of points inside one image
        inside_coord = inside_coord + inside_coord_im

    sat_path = np.array(inside_coord).reshape(np.shape(stacked_im)) # reshaping satellite path
    im_mosaic_tmp = im_mosaic_points.reshape(np.shape(stacked_im)) # reshaping interpolated image
    im_mosaic=np.where(sat_path,im_mosaic_tmp,np.nan) # only keeping interpolated pixels that are in the satellite imaging path

    return (im_mosaic)
# %%
