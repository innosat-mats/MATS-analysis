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

start_time = DT.datetime(2023, 4, 2, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 3, 0, 0, 0)

start_time = DT.datetime(2023, 4, 7, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 8, 0, 0, 0)

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

# df1a = df1a_tot

#df1a = df1a_tot[:20]
# df1a = corrected_ccd
# df1a = df_corr

# df1a.dropna(inplace=True)

df1a = df1a_mask_1day

df1a = df1a[~np.isnan(df1a['satlat'])]

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


#%%


NADIR_AZ = [] # list of nadir solar azimuth angles
n = len(df1a)
for i in tqdm(range(n)):
    ccditem = df1a.iloc[i]
    #print(i)
    NADIR_AZ.append(nadir_az(ccditem))      

NADIR_AZ = np.array(NADIR_AZ)
plt.plot(df1a['EXPDate'],NADIR_AZ,'.')

# %% 1 pixel, 1 angle range

ang_min = -90
ang_max = ang_min+0.5
step = 3
x_art = 17
y_art = 10
x_dark = x_art
y_dark = y_art - 8


X = im_points[step:-1,y_dark,x_dark]
Y = im_points[0:-(step+1),y_art,x_art]



# indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['nadir_sza']>ang_min) & (df1a.iloc[:-(step+1)]['nadir_sza']<ang_max)
indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (NADIR_AZ[:-(step+1)]>ang_min) & (NADIR_AZ[:-(step+1)]<ang_max)

X = X[indexes]
Y = Y[indexes]
# C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']
C = NADIR_AZ[np.append(indexes,(step+1)*[False])]
# C = df1a.iloc[np.append(indexes,(step+1)*[False])]['nadir_sza']

from scipy.optimize import curve_fit






# def func(X,b,a):
#     return X*a + b*np.ones_like(X)

# fit_param, cov = curve_fit(func,X,Y)
# abs_err = Y-func(X,fit_param[0],fit_param[1])
# rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
# intercept = fit_param[0]
# slope = fit_param[1]


def func(X,b):
    return X+b

fit_param, cov = curve_fit(func,X,Y)
abs_err = Y-func(X,fit_param[0])
rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
intercept = fit_param[0]
slope = 1.0


plt.figure()
plt.scatter(X,Y,c=C)
plt.plot(X,intercept + slope*X,color='red')
plt.title(f"Offset of {step} images. Slope = {slope:.3f}, intercept = {intercept:.1f}, R**2 = {rsquare:.3f}")
# plt.title(f"Artifact correlation, offset of {step} images. Slope = {slope:.3f}, intercept = {intercept:.1f}, R**2 = {rsquare:.3f}")
plt.xlabel(f'Pixel ({x_dark},{y_dark}) (correct value)')
plt.ylabel(f'Pixel ({x_art},{y_art}) (artifact)')
plt.colorbar(label='nadir azimuth')
plt.show()



# %% all the image, 1 angle range
ang_min = -80.17
ang_max = -80.07
im_R = np.ones_like(im_points[0,:,:])
im_bias = np.zeros_like(im_points[0,:,:])
im_slope = np.ones_like(im_points[0,:,:])

for x_art in range(0,56):
    for y_art in range(3,14):

        step = int(y_art//2.6)
        x_dark = x_art
        y_dark = int(round(y_art%2.6))

        X = im_points[step:-1,y_dark,x_dark]
        Y = im_points[0:-(step+1),y_art,x_art]
     

        # indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['TPssa']>ang_min) & (df1a.iloc[:-(step+1)]['TPssa']<ang_max)
        indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (NADIR_AZ[:-(step+1)]>ang_min) & (NADIR_AZ[:-(step+1)]<ang_max)


        X = X[indexes]
        Y = Y[indexes]
        # C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']
        C = NADIR_AZ[np.append(indexes,(step+1)*[False])]
        # C = df1a.iloc[np.append(indexes,(step+1)*[False])]['nadir_sza']

        from scipy.optimize import curve_fit

        # plt.figure()
        # plt.scatter(X,Y,c=C)

        # def func(X,b,a):
        #     return X*a + b*np.ones_like(X)

        # fit_param, cov = curve_fit(func,X,Y)
        # abs_err = Y-func(X,fit_param[0],fit_param[1])
        # rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
        # intercept = fit_param[0]
        # slope = fit_param[1]


        def func(X,b):
            return X+b

        fit_param, cov = curve_fit(func,X,Y)
        abs_err = Y-func(X,fit_param[0])
        rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
        intercept = fit_param[0]
        slope = 1.0

        im_R[y_art,x_art] = rsquare
        im_bias[y_art,x_art] = intercept
        im_slope[y_art,x_art] = slope

plt.figure()
plt.title(f'Bias value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.title(f'Bias value')
plt.imshow(im_bias[:,:],origin='lower')
plt.colorbar()
plt.show()

plt.figure()
plt.title(f'Slope value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.imshow(im_slope[:,:],origin='lower')
plt.colorbar()
plt.show()

plt.figure()
plt.title(f'R squared value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.title(f'R squared value')
plt.imshow(im_R[:,:],origin='lower',vmin = 0.8,vmax=1)
plt.colorbar()
plt.show()

plt.figure()
plt.title(f'Bias histogram ')
plt.hist(im_bias.ravel(),30)
plt.show()

# %% 1 pixel, several angles
ang_min = -105
ang_max = -71
ang_step = 0.1
ANG = np.linspace(ang_min,ang_max,int((ang_max-ang_min)//ang_step))
BIAS = []
R2 = []
step = 3
x_art = 17
y_art = 10
x_dark = x_art
y_dark = y_art - 8
for ang in ANG:
    ang_min = ang
    ang_max = ang + ang_step
    X = im_points[step:-1,y_dark,x_dark]
    Y = im_points[0:-(step+1),y_art,x_art]

    # indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['nadir_sza']>ang_min) & (df1a.iloc[:-(step+1)]['nadir_sza']<ang_max)
    indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (NADIR_AZ[:-(step+1)]>ang_min) & (NADIR_AZ[:-(step+1)]<ang_max)

    X = X[indexes]
    Y = Y[indexes]
    # C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']
    C = NADIR_AZ[np.append(indexes,(step+1)*[False])]


    from scipy.optimize import curve_fit




    # def func(X,b,a):
    #     return X*a + b*np.ones_like(X)

    # fit_param, cov = curve_fit(func,X,Y)
    # abs_err = Y-func(X,fit_param[0],fit_param[1])
    # rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
    # intercept = fit_param[0]
    # slope = fit_param[1]


    def func(X,b):
        return X+b

    fit_param, cov = curve_fit(func,X,Y)
    abs_err = Y-func(X,fit_param[0])
    rsquare = 1.0 - (np.var(abs_err)/np.var(Y))
    intercept = fit_param[0]
    slope = 1.0

    BIAS.append(intercept)
    R2.append(rsquare)


plt.figure()
plt.title(f'Bias value (Pixel ({x_art},{y_art}))')
plt.title(f'Bias value')
plt.plot(ANG,BIAS,'.')
plt.xlabel('Azimuth angle (deg)')
plt.ylabel('Bias value')
plt.show()

plt.figure()
plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.title(f'R squared value')
plt.plot(ANG,R2,'.')
plt.xlabel('Azimuth angle (deg)')
plt.ylabel('R squared value')
plt.ylim([0.5,1.05])
plt.show()


# %%  dependency 


ang_min = -360
ang_max = ang_min + 500
step = 3
x_art = 17
y_art = 10
x_dark = x_art
y_dark = y_art - 8


X = im_points[step:-1,y_dark,x_dark]
Y = im_points[0:-(step+1),y_art,x_art]


indexes = (~np.isnan(X)) & (~np.isnan(Y))
indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['nadir_sza']>ang_min) & (df1a.iloc[:-(step+1)]['nadir_sza']<ang_max)
indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (NADIR_AZ[:-(step+1)]>ang_min) & (NADIR_AZ[:-(step+1)]<ang_max)

X = X[indexes]
Y = Y[indexes]
TPSSA = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']
AZ = NADIR_AZ[np.append(indexes,(step+1)*[False])]
SZA = df1a.iloc[np.append(indexes,(step+1)*[False])]['nadir_sza']
EXP = df1a.iloc[np.append(indexes,(step+1)*[False])]['EXPDate']


# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(TPSSA,Y-X,'.')
plt.xlabel('TPSSA angle (deg)')
plt.ylabel('Pixel value difference (artifact-ref)')
plt.show()

# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(AZ,Y-X,'.')
plt.xlabel('Azimuth angle (deg)')
plt.ylabel('Pixel value difference (artifact-ref)')
plt.show()

# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(SZA,Y-X,'.')
plt.xlabel('Solar zenith angle (deg)')
plt.ylabel('Pixel value difference (artifact-ref)')
plt.show()

# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(AZ,SZA,'.')
plt.xlabel('Azimuth angle (deg)')
plt.ylabel('Solar zenith angle (deg)')
plt.show()

# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(EXP,Y-X,'.')
plt.xlabel('Exposure time')
plt.ylabel('Pixel value difference (artifact-ref)')
plt.show()

# plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.figure()
plt.plot(AZ,Y-X,'.')
plt.xlabel('Solar azimuth angle')
plt.ylabel('Pixel value difference (artifact-ref)')
plt.show()
# %%
