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


#%% local variables
# value of a saturated pixel 
sat_val = 32880

# times for start and stop
start_time = DT.datetime(2023, 3, 20, 12, 0, 0)
stop_time = DT.datetime(2023, 3, 21, 13, 0, 0)

start_time = DT.datetime(2023, 3, 22, 3, 0, 0)
stop_time = DT.datetime(2023, 3, 22, 3, 15, 0)

start_time = DT.datetime(2022, 12, 20, 0, 0, 0)
stop_time = DT.datetime(2022, 12, 22, 0, 0, 0)

start_time = DT.datetime(2022, 12, 11, 0, 0, 0)
stop_time = DT.datetime(2022, 12, 13, 0, 0, 0)

start_time = DT.datetime(2022, 12, 21, 4, 0, 0)
stop_time = DT.datetime(2022, 12, 21, 8, 0, 0)

start_time = DT.datetime(2023, 4, 1, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 2, 0, 0, 0)

# start_time = DT.datetime(2022, 12, 21, 0, 0, 0)
# stop_time = DT.datetime(2022, 12, 22, 0, 0, 0)

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}


#%% reading mesurements
df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
print(len(df1a_tot))
df1a = df1a_tot

# displaying keys
pd.set_option('display.max_rows', 100)
df1a.dtypes

plt.plot(range(len(df1a)),df1a['EXPDate'])



#%% computing solar zenith angles for each pixel
# df1a = df1a_tot[:]

#df1a = df1a_tot[:20]
n = len(df1a)
a,b = np.shape(df1a.iloc[0]['IMAGE'])
im_points = np.zeros((n,a,b))

for i in tqdm(range(n)):
    ccditem = df1a.iloc[i]
    im = ccditem['IMAGE']
    im_points [i,:,:] = im




#%%


def nadir_az(ccditem):
    """
    Function giving various angles.. 

   
    Arguments:
        ccditem or dataframe with the 'EXPDate'

    Returns:
        nadir_sza: solar zenith angle at satelite position (degrees)
        TPsza: solar zenith angle at TP position (degrees)
        TPssa: solar scattering angle at TP position (degrees), 
        tpLT: Local time at the TP (string)
        
    """
    planets=load('de421.bsp')
    earth,sun,moon= planets['earth'], planets['sun'],planets['moon']
   
    d = ccditem['EXPDate']
    ts =load.timescale()
    t = ts.from_datetime(d)
    satlat, satlon, satheight = coordinates.satpos(ccditem)
    sat_pos=earth + wgs84.latlon(satlat, satlon, elevation_m=satheight)
    sundir=sat_pos.at(t).observe(sun).apparent()
    obs=sundir.altaz()
    nadir_az = (obs[1].degrees) #nadir solar zenith angle    
    return nadir_az



def nadir_mask(ccditems,bias_threshold):
    """
    Function giving various angles.. 

   
    Arguments:
        ccditem or dataframe with the 'EXPDate'

    Returns:
        nadir_sza: solar zenith angle at satelite position (degrees)
        TPsza: solar zenith angle at TP position (degrees)
        TPssa: solar scattering angle at TP position (degrees), 
        tpLT: Local time at the TP (string)
        
    """

    # only using NADIR images to set the mask values
    ccditems = ccditems[ccditems['CCDSEL'] == 7]
    pixel_shift = 2.6 # pixel shift along the y axis between 2 consecutive images

    n = len(ccditems)
    if n == 0:
        return None, None
    
    a,b = np.shape(ccditems.iloc[0]['IMAGE'])
    im_points = np.zeros((n,a,b))
    NADIR_AZ = []

    for i in range(n):
        ccditem = df1a.iloc[i]
        im = ccditem['IMAGE']
        im_points [i,:,:] = im
        
    
    im_R2 = np.ones_like(im_points[0,:,:])
    im_bias = np.zeros_like(im_points[0,:,:])
    

    def func(X,b):
        return X+b

    for x_art in range(0,b):
        for y_art in range(3,a):

            step = int(y_art//pixel_shift)
            x_dark = x_art
            y_dark = int(round(y_art%pixel_shift))

            X = im_points[step:-1,y_dark,x_dark]
            Y = im_points[0:-(step+1),y_art,x_art]                

            

            fit_param, cov = curve_fit(func,X,Y)
            abs_err = Y-func(X,fit_param[0])
            rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
            intercept = fit_param[0]
            slope = 1.0

            im_R2[y_art,x_art] = rsquare
            im_bias[y_art,x_art] = intercept

    mask = im_bias>bias_threshold       
    bias_mask = im_bias * mask
    R2_mask = im_R2 * mask

    return bias_mask,R2_mask


def mask_correction(ccditems,mask):
    n = len(ccditems)
    corrected_ccditems = ccditems.copy()
    corrected_ccditems['IMAGE'] = [im - mask for im in ccditems['IMAGE']]
    return corrected_ccditems

    


def azimuth_bias_mask(ccditems,bias_threshold,az_min=0,az_max=360,az_step=5):
    
    # only using NADIR images to set the mask values
    ccditems = ccditems[ccditems['CCDSEL'] == 7]

    NADIR_AZ = [] # list of nadir solar azimuth angles
    
    n = len(ccditems)
    a,b = np.shape(ccditems.iloc[0]['IMAGE'])

    IM_BIAS = []
    IM_R2 = []
    AZ_MIN = []
    AZ_MAX = []

    for i in tqdm(range(n)):
        ccditem = df1a.iloc[i]
        NADIR_AZ.append(nadir_az(ccditem))      

    NADIR_AZ = np.array(NADIR_AZ)

    for ang in range(az_min,az_max,az_step):
        ang_min = ang
        ang_max = ang + az_step
        
        indexes = (NADIR_AZ>ang_min) & (NADIR_AZ<ang_max)

        bias_mask,R2_mask = nadir_mask(ccditems[indexes],bias_threshold)

        IM_BIAS.append(bias_mask)
        IM_R2.append(R2_mask)
        AZ_MIN.append(ang_min)
        AZ_MAX.append(ang_max)
    
    azimuth_masks = pd.DataFrame({'bias_mask': IM_BIAS,
                                 'R2_mask': IM_R2,
                                 'azimuth_min': AZ_MIN,
                                 'azimuth_max': AZ_MAX})
    return azimuth_masks



# def azimuth_corr_mask(ccditems, azimuth_masks):
    
#     m = len(azimuth_masks)
#     n = len(ccditems)
#     corrected_ccditems = ccditems.copy()

#     NADIR_AZ = [] # list of nadir solar azimuth angles
    
#     n = len(ccditems)
#     a,b = np.shape(ccditems.iloc[0]['IMAGE'])

#     for index in tqdm(ccditems.index):
#         azimuth = nadir_az(ccditems.loc[index])
#         im = ccditems.loc[index,'IMAGE']
#         for j in range(m):
#             indexes = (azimuth_masks['azimuth_min']<azimuth) & (azimuth_masks['azimuth_max']>azimuth)
#             mask = azimuth_masks['bias_mask'][indexes].iloc[0]
#         corrected_ccditems.loc[index,'IMAGE'] = im - mask

#     return corrected_ccditems

        
def azimuth_corr_mask(ccditems, azimuth_masks):
    
    m = len(azimuth_masks)
    n = len(ccditems)
    corrected_ccditems = ccditems.copy()

    NADIR_AZ = [] # list of nadir solar azimuth angles
    
    n = len(ccditems)
    a,b = np.shape(ccditems.iloc[0]['IMAGE'])

    IMAGES = []

    for index in tqdm(ccditems.index):
        azimuth = nadir_az(ccditems.loc[index])
        im = ccditems.loc[index,'IMAGE']
        for j in range(m):
            indexes = (azimuth_masks['azimuth_min']<azimuth) & (azimuth_masks['azimuth_max']>azimuth)
            mask = azimuth_masks['bias_mask'][indexes].iloc[0]
            
        IMAGES.append(im - mask)

    corrected_ccditems['IMAGE'] = IMAGES

    return corrected_ccditems

def artifact_correction(ccd_items,ccd_items_correction,az_min=0,az_max=360,az_step=5):
    return 4


azimuth_masks = azimuth_bias_mask(df1a,500)

df_corr = azimuth_corr_mask(df1a,azimuth_masks)



