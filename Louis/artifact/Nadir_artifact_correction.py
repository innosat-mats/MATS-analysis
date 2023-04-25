# script for studying the number of saturated pixels in NADIR images as a function of solare zenith angle 

#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from mats_utils.statistiscs.images_functions import create_imagecube
from tqdm import tqdm
from mats_utils.geolocation.coordinates import NADIR_geolocation 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from skyfield.api import wgs84 
from skyfield.api import load
from numpy.linalg import norm
from scipy.optimize import curve_fit
from datetime import datetime, timedelta
import warnings
from skyfield.positionlib import Geocentric, ICRF
from skyfield.units import Distance




#%%


def nadir_az(ccditem):
    """
    Function giving the solar azinuth angle  

   
    Arguments:
        ccditem 
    Returns:
        nadir_az: solar azimuth angle at nadir imager (degrees)       
        
    """
    planets=load('de421.bsp')
    earth,sun,moon= planets['earth'], planets['sun'],planets['moon']
   
     
    d = ccditem['EXPDate']
    ts =load.timescale()
    t = ts.from_datetime(d)
    satlat, satlon, satheight = coordinates.satpos(ccditem)
    TPlat, TPlon, TPheight = coordinates.TPpos(ccditem)
    
    sat_pos=earth + wgs84.latlon(satlat, satlon, elevation_m=satheight)
    TP_pos=earth + wgs84.latlon(TPlat, TPlon, elevation_m=TPheight)
    sundir=sat_pos.at(t).observe(sun).apparent()
    limbdir = TP_pos.at(t) - sat_pos.at(t)
    obs_limb = limbdir.altaz()
    obs_sun=sundir.altaz()
    nadir_az = (obs_sun[1].degrees - obs_limb[1].degrees) #nadir solar azimuth angle    
    return nadir_az



def nadir_mask(ccditems,bias_threshold):
    """
    Function calculating linear regressions for each pixel in the nadir images in the dataframe. The 
    regression is performed on the correlation between a pixel in the corrected zone and the corresponding
    pixel in the reference zone. As the ground is moving under the satellite, a same object below is seen 
    by the corrected pixels and some exposures later by the reference pixel.
    Bias values under the threshhold value are set to zero.
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the images
        bias_threshhold : float
            all bias values smaller than this value are not taken into account (set to zero)

    Returns:
        bias_mask : np.array[float]
            array of bias values between the "corrected zone" and the "reference zone" in the image. The array
            has the same shape as the images (in the referecence zone)
        R2_mask : np.array[float]
            array of R2 values in the regression between the "corrected zone" and the "reference zone" in the 
            image. The array has the same shape as the images (in the referecence zone)
        
    """

    # only using NADIR images to set the mask values
    ccditems = ccditems[ccditems['CCDSEL'] == 7]
    pixel_shift = 2.6 # pixel shift along the y axis between 2 consecutive images. For a sampling time of 2s it is equivalent to 2.6 pixels

    n = len(ccditems)
    if n == 0:
        return None, None
    
    a,b = np.shape(ccditems.iloc[0]['IMAGE'])            
    
    im_R2 = np.ones((a,b)) # mask of R2 values
    im_bias = np.zeros((a,b)) # mask of bias values    
    

    def func(X,b):
        return X+b

    for x_art in range(0,b):
        for y_art in range(3,a): # skipping the first 3 "reference" rows

            # calculating the position of the reference pixel and the number of images offset
            step = int(y_art//pixel_shift) # images offset between the reference and corrected pixel
            x_dark = x_art
            y_dark = int(round(y_art%pixel_shift))

            X = []
            Y = []

            # Data points for the regression
            for i in range(n-step):
                ccditem_art = ccditems.iloc[i]
                ccditem_ref = ccditems.iloc[i+step]

                #check if the image taken as reference is indeed taken the right amount of time after the other image
                if (ccditem_ref['EXPDate'] - ccditem_art['EXPDate'] - step*timedelta(seconds=2)) < timedelta(seconds= 0.01):
                    X.append(ccditem_ref['IMAGE'][y_dark,x_dark])
                    Y.append(ccditem_art['IMAGE'][y_art,x_art])
            

            if len(X) + len(Y) > 0:
                # linear regression, the slope is set to 1 
                fit_param, cov = curve_fit(func,X,Y)
                abs_err = Y-func(X,fit_param[0])
                rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
                intercept = fit_param[0]
                slope = 1.0
            else :                 
                return None,None 


            im_R2[y_art,x_art] = rsquare
            im_bias[y_art,x_art] = intercept

    # thresholding the correction
    mask = im_bias>bias_threshold       
    bias_mask = im_bias * mask
    R2_mask = im_R2 * mask

    return bias_mask,R2_mask



def mask_correction(ccditems,mask):
    """
    Function applying a correction mask to all the images in a dataframe. The values in 
    the mask are substracted to the pixel values in the ccditems 
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the ccditems 
        mask : np.array[float]
            array with the same shape as the images in the dataframe. 

    Returns:
        corrected_ccditems : Panda dataframe
            dataframe containing the corrected ccditems      
    """
    corrected_ccditems = ccditems.copy()
    corrected_ccditems['IMAGE'] = [im - mask for im in ccditems['IMAGE']]
    return corrected_ccditems

    


def azimuth_bias_mask(ccditems,bias_threshold,az_list):
    """
    Function creating correction masks dependant on solar azimuth angles. The mask
    creation follows the same rules as in the function nadir_mask. The bias and R2 
    masks are created for each azimuth angle intervall.     
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the images
        bias_threshhold : float
            all bias values smaller than this value are not taken into account (set to zero)
        az_list : list of float
            list of azimuth value. The regression is made on azimuth angle intervalls with the 
            given angles as center points 

    Returns:
        azimuth_masks : Pandas dataframe
            'bias_mask' : np.array[float]
                thresholded bias for each pixel, has the same shape as the images
            'R2_mask' : np.array[float]
                R2 value for each pixel, has the same shape as the images
            'azimuth' : float
                center value of each azimuth intervall. This might not be the case 
            
    """
    # only using NADIR images to set the mask values
    ccditems = ccditems[ccditems['CCDSEL'] == 7]
    n = len(ccditems)

    NADIR_AZ = [] # list of nadir solar azimuth angles
    
    IM_BIAS = [] # list of bias masks
    IM_R2 = []  # list of R2 masks
   
    for i in tqdm(range(n)):
        ccditem = ccditems.iloc[i]
        NADIR_AZ.append(nadir_az(ccditem))      

    NADIR_AZ = np.array(NADIR_AZ)

    
    az_max = max(NADIR_AZ)

    

    for i in tqdm(range(len(az_list))):
        
        if i == 0:
            az_min = min(NADIR_AZ) # lower interval limit
        else :
            az_min = (az_list[i-1] + az_list[i])/2.0 # lower interval limit

        if i == len(az_list)-1:
            az_max = max(NADIR_AZ) # higher interval limit
        else :
            az_max = (az_list[i+1] + az_list[i])/2.0 # higher interval limit
        
        
        # indexes of ccditems with an angle between ang_min and ang_max
        indexes = (NADIR_AZ>az_min) & (NADIR_AZ<az_max)

        if np.any(indexes): # if ccditems exist in this angle interval
            # regression on the selected images
            bias_mask,R2_mask = nadir_mask(ccditems[indexes],bias_threshold)
        else: # if no image is found
            bias_mask = None
            R2_mask = None

        IM_BIAS.append(bias_mask)
        IM_R2.append(R2_mask)
        
    
    # creating dataframe
    azimuth_masks = pd.DataFrame({'bias_mask': IM_BIAS,
                                 'R2_mask': IM_R2,
                                 'azimuth': az_list})
    return azimuth_masks




def azimuth_bias_mask2(ccditems,bias_threshold,az_list):
    """
    Function creating correction masks dependant on solar azimuth angles. The mask
    creation follows the same rules as in the function nadir_mask. The bias and R2 
    masks are created for each azimuth angle intervall.     
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the images
        bias_threshhold : float
            all bias values smaller than this value are not taken into account (set to zero)
        az_list : list of float
            list of azimuth value. The regression is made on azimuth angle intervalls with the 
            given angles as center points 

    Returns:
        azimuth_masks : Pandas dataframe
            'bias_mask' : np.array[float]
                thresholded bias for each pixel, has the same shape as the images
            'R2_mask' : np.array[float]
                R2 value for each pixel, has the same shape as the images
            'azimuth' : float
                center value of each azimuth intervall. This might not be the case 
            
    """
    # only using NADIR images to set the mask values
    ccditems = ccditems[ccditems['CCDSEL'] == 7]
    n = len(ccditems)
    a,b = np.shape(ccditems.iloc[0]['IMAGE'])

    pixel_shift = 2.6 # pixel shift along the y axis between 2 consecutive images. For a sampling time of 2s it is equivalent to 2.6 pixels
    
    im_points = np.zeros((n,a,b)) # array containing all the pixel values
    NADIR_AZ = [] # list of nadir solar azimuth angles
    EXP_DATE = [] # list of exposition dates
    
    IM_BIAS = [] # list of bias masks
    IM_R2 = []  # list of R2 masks
   
    # filling the several arrays
    for i in tqdm(range(n)):
        ccditem = ccditems.iloc[i]
        NADIR_AZ.append(nadir_az(ccditem)) 
        im = ccditem['IMAGE']
        im_points [i,:,:] = im
        EXP_DATE.append(ccditem['EXPDate'])

    NADIR_AZ = np.array(NADIR_AZ)
    EXP_DATE = np.array(EXP_DATE)

    # intermediate function for regression (only the bias is optimized)
    def func(X,b):
        return X+b
    

    for i in tqdm(range(len(az_list))): # looping on the azimuth angle intervals
        
        # lower interval limit
        if i == 0:
            az_min = min(NADIR_AZ) 
        else :
            az_min = (az_list[i-1] + az_list[i])/2.0

        # higher interval limit
        if i == len(az_list)-1:
            az_max = max(NADIR_AZ) 
        else :
            az_max = (az_list[i+1] + az_list[i])/2.0 
        
        
        # regression on the selected images
        im_R2 = np.ones_like(im_points[0,:,:]) # R2 values for each pixel
        im_bias = np.zeros_like(im_points[0,:,:]) # bias value for each pixel

        for x_art in range(0,56):
            for y_art in range(3,14): # the first 3 rows are taken as reference

                # calculating the position of the reference pixel and the number of images offset
                step = int(y_art//pixel_shift) # images offset between the reference and corrected pixel
                x_ref= x_art
                y_ref = int(round(y_art%2.6))

                X = [] # list of reference pixel values
                Y = [] # list of artifact pixel values

                az_index_art = (az_min < NADIR_AZ) & (NADIR_AZ < az_max) # indices for the images with in the given azimuth angle interval (the pixels to be corrected are taken from these images)
                art_indices = np.arange(n)[az_index_art][:-step] # indices of the corresponding reference images where the reference pixel is taken
                
                if len(art_indices) > 0 : 
                    for art_ind in art_indices: # iterating over the images (image index of the artifact pixel)
                        ref_ind = art_ind + step # image index of the reference 
                        #check if the image taken as reference is indeed taken the right amount of time after the other image
                        if (EXP_DATE[art_ind] - EXP_DATE[ref_ind] - step*timedelta(seconds=2)) < timedelta(seconds= 0.01):
                            X.append(im_points[ref_ind,y_ref,x_ref])
                            Y.append(im_points[art_ind,y_art,x_art])

                if len(X) + len(Y) > 0:
                # linear regression, the slope is set to 1 
                    fit_param, cov = curve_fit(func,X,Y)
                    abs_err = Y-func(X,fit_param[0])
                    rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
                    intercept = fit_param[0]
                    im_R2[y_art,x_art] = rsquare
                    im_bias[y_art,x_art] = intercept
                else : 
                    im_R2[y_art,x_art] = None
                    im_bias[y_art,x_art] = None               
        
        # thresholding the correction
        mask = im_bias>bias_threshold       
        bias_mask = im_bias * mask
        R2_mask = im_R2 * mask
            
        IM_BIAS.append(bias_mask)
        IM_R2.append(R2_mask)
        
    
    # creating dataframe
    azimuth_masks = pd.DataFrame({'bias_mask': IM_BIAS,
                                 'R2_mask': IM_R2,
                                 'azimuth': az_list})
    return azimuth_masks


      
def azimuth_corr_mask(ccditems, azimuth_masks):
    """
    Function applying a correction mask to all the images in a dataframe. The values in 
    the mask are substracted to the pixel values in the ccditems. The correction mask is
    a function of the solar azimuth angle. 
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the ccditems 
        azimuth_masks : Panda dataframe
            dataframe containing the masks to apply, sorted by azimuth angle. The structure 
            of the dataframe is the same as the one described for the output of the azimuth_bias_mask function

    Returns:
        corrected_ccditems : Panda dataframe
            dataframe containing the corrected ccditems      
    """

    m = len(azimuth_masks)
    corrected_ccditems = ccditems.copy()

    IMAGES = [] # list of corrected images

    for index in tqdm(ccditems.index):
        azimuth = nadir_az(ccditems.loc[index])
        im = ccditems.loc[index,'IMAGE'] # image to be corrected
        
        # finding the mask which corresponding azimuth angle interval is the closest to the image's azimuth angle
        distance = 360.0
        best_ind = 0
        for j in range(m):
            if abs(azimuth_masks.iloc[j]['azimuth']-azimuth) < distance and (type(azimuth_masks.iloc[j]['bias_mask']) != type(None)):
                best_ind = j
                distance = abs(azimuth_masks.iloc[j]['azimuth']-azimuth)
        mask = azimuth_masks['bias_mask'][best_ind]

        # substracting the mask
        IMAGES.append(im - mask)

    # creating the corrected dataframe
    corrected_ccditems['IMAGE'] = IMAGES

    return corrected_ccditems

def artifact_correction(ccditems,ccditems_correction,bias_threshold=600.0,az_step=5.0):
    """
    Function computing and applying a correction mask on the nadir images. The correction masks are computed 
    by assuming a constant bias between the expected pixel value and the measured one in the artifact. Several 
    azimuth angles intervals are defined and a mask is computed for each interval.
   
    Arguments:
        ccditems : Panda dataframe
            dataframe containing the ccditems to be corrected
        ccditems_correction : Panda dataframe
            dataframe containing the ccditems used to compute the masks
        bias_threshold : float
            all bias values smaller than this value are not taken into account (set to zero)
        az_step : float
            number of degrees in each azimuth angle interval in the masks

    Returns:
        corrected_ccditems : Panda dataframe
            dataframe containing the corrected ccditems           
    """


    print(f"\n Calculating the correction mask")
    azimuth_masks = azimuth_bias_mask(ccditems_correction,500,az_step=5)

    print(f"\n Correcting the images")
    corrected_ccditems = azimuth_corr_mask(ccditems,azimuth_masks)

    return corrected_ccditems


#%%
def azimuth_masks_plot(azimuth_masks,angles):
    m = len(azimuth_masks)
    
    for ang in angles:
        distance = 360.0
        best_ind = 0
        for j in range(m):
            if abs(azimuth_masks.iloc[j]['azimuth']-ang) < distance and (type(azimuth_masks.iloc[j]['bias_mask']) != type(None)):
                best_ind = j
                distance = abs(azimuth_masks.iloc[j]['azimuth']-ang)
        mask = azimuth_masks['bias_mask'][best_ind]
        R2 = azimuth_masks['R2_mask'][best_ind]
        az = azimuth_masks['azimuth'][best_ind]
        

        
        fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
        fig.suptitle(f" solar azimuth angle : {az:.2f} deg")
        fig = ax1.imshow(mask,origin='lower')
        ax1.set_title('bias mask')
        plt.colorbar(fig,ax=ax1,fraction=0.02)
            
        fig = ax2.imshow(R2,vmin=0.5,vmax=1.0,origin='lower')
        plt.colorbar(fig,ax=ax2,fraction=0.02)
        ax2.set_title('R2 values')
        plt.show()
        










# %%

