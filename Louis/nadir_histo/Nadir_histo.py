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

start_time = DT.datetime(2023, 3, 31, 0, 0, 0)
stop_time = DT.datetime(2023, 4, 1, 0, 0, 0)

start_time = DT.datetime(2022, 12, 21, 0, 0, 0)
stop_time = DT.datetime(2022, 12, 22, 0, 0, 0)

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


#%% computing the saturation level of a ccd item
def saturation_level(ccditem, sat_val = 32880):
    """
    Function to get the proportion of saturated pixels in an image

    Parameters
    ----------
    ccditem : CCDitem
        measurement
    sat_val : int
        value of a saturated pixel
        
    Returns
    -------
    sat_prop : float
        proportion (between 0 and 1) of saturated pixels
    """
    im = ccditem['IMAGE']
    a,b = np.shape(im)
    saturated_im = sat_val*np.ones_like(im) # completely saturated image
    SAT_LEVEL = im == saturated_im # boolean array : True value for saturated pixels, False otherwise
    sat_prop = SAT_LEVEL.sum()/(a*b) # proportion of saturated pixels
    return(sat_prop)

#%% computing solar zenith angles for each image 
NADIR_SZAS = [] # list of nadir solar zenih angles
SAT_LEV = [] # list of saturation level (proportion of saturated pixels)
LAT = [] # list of latitude points
n = len(df1a)
for i in tqdm(range(n)):
    ccditem = df1a.iloc[i]
    nadir_sza, TPsza, TPssa, TPlt = coordinates.angles(ccditem)   
    NADIR_SZAS.append(nadir_sza) 
    SAT_LEV.append(saturation_level(ccditem)) 
    LAT.append(ccditem['satlat'])

print(f"{np.average(SAT_LEV)*100:.1f} % of saturated pixels over the whole image set")
print(f"{np.sum(SAT_LEV==np.ones_like(SAT_LEV))/len(SAT_LEV)*100:.1f} % of completely saturated images over the whole image set")


#%% plotting
# plotting NADIR SZA
plt.figure()
plt.plot(df1a['EXPDate'],NADIR_SZAS,linestyle=' ',marker='.')
plt.ylabel('NADIR SZA (deg)')

# plotting Latitude and SZA
plt.figure()
plt.plot(NADIR_SZAS,LAT,linestyle=' ',marker='.')
m = len(LAT)
plt.arrow(NADIR_SZAS[m//2-1],LAT[m//2-1],NADIR_SZAS[m//2]-NADIR_SZAS[m//2-1],LAT[m//2]-LAT[m//2-1], shape='full', lw=0, length_includes_head=False, head_width=.5, head_length=5)
plt.ylabel('NADIR latitude')
plt.xlabel('NADIR SZA (deg)')

plt.show()

#%% computing solar zenith angles for each pixel
win_mod = '13..2'
df1a = df1a_tot[df1a_tot['WDWInputDataWindow'] == win_mod]
df1a = df1a_tot[:]

#df1a = df1a_tot[:20]
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
    im_points [i,:,:] = im


#%% studying only the pixels with sza > 105 deg 

ang_lim = 95
im_points_night = np.zeros_like(im_points)
sza_points_night = np.zeros_like(sza_points)
time_points_night = np.zeros_like(im_points)
n,a,b = np.shape(im_points)

for j in tqdm(range(n)):
    for k in range(a):
        for l in range(b):
            if sza_points[j,k,l] <= ang_lim :
                im_points_night[j,k,l] = None
                sza_points_night[j,k,l] = None  
            else :
                im_points_night[j,k,l] = im_points[j,k,l]
                sza_points_night[j,k,l] = sza_points[j,k,l]
                         
                
stacked_im = np.nanmean(im_points_night,axis=0)
mask = stacked_im < np.nanmean(stacked_im) + 1000 # arbitrary
correction1 = -(stacked_im-np.nanmean(im_points_night))*(1-mask)
correction2 = -(stacked_im-np.nanmean(im_points_night[:,1,17]))*(1-mask)

plt.close('Stacked image')
plt.figure('Stacked image')
plt.imshow(stacked_im)
plt.colorbar()
plt.show()      

plt.close('Stacked image - mean')
plt.figure('Stacked image - mean')
plt.imshow(stacked_im-np.nanmean(im_points_night))
plt.colorbar()
plt.show() 

plt.close('Correction 1')
plt.figure('Correction 1')
plt.imshow(correction1)
plt.colorbar()
plt.show()  

plt.close('Correction 2')
plt.figure('Correction 2')
plt.imshow(correction2)
plt.colorbar()
plt.show()  

plt.close('Correction 3')
plt.figure('Correction 3')
plt.imshow(mask)
plt.colorbar()
plt.show()  

plt.close('Pixel value')
plt.figure('Pixel value')
plt.plot(im_points_night[:,10,17], label='pixel(10,17)')
plt.plot(im_points_night[:,1,17], label='pixel(1,17)')
#plt.plot(im_points_night[:,10,17]-np.nanmean(im_points_night,axis=(1,2)), label='pixel(10,17) - mean')
#plt.plot(im_points_night[:,10,17]-np.nanmean(im_points_night[:,10,17]), label='pixel(10,17) - mean(10,17)')
#plt.plot(im_points_night[:,10,17]-np.nanmean(im_points_night), label='pixel(10,17) - mean (over all the pictures)')
plt.plot(im_points_night[:,10,17]-(np.nanmean(im_points_night[:,10,17])-np.nanmean(im_points_night)), label='pixel(10,17) - [mean(10,17)-mean]')
plt.plot(im_points_night[:,10,17]-(np.nanmean(im_points_night[:,10,17])-np.nanmean(im_points_night[:,1,17])), label='pixel(10,17) - [mean(10,17)-mean]')
plt.plot(np.nanmean(im_points_night,axis=(1,2)), label='mean pixel value')
plt.ylabel('Pixel value')
plt.legend()


plt.show()

# %%
# DataWindow = ['13..2','14..3','15..4']


# im_points_cor1 = np.zeros_like(im_points) # substracting (mean image - mean pixel value)
# im_points_cor2 = np.zeros_like(im_points) # substracting (mean image - mean dark pixel value)
# im_points_cor3 = np.zeros_like(im_points) # removing pixels with artefacts
# for j in tqdm(range(n)):
#     for k in range(a):
#         for l in range(b):
#             if im_points_night[j,k,l] != None:
#                 im_points_cor1[j,k,l] = im_points_night[j,k,l] + correction1[k,l]
#                 im_points_cor2[j,k,l] = im_points_night[j,k,l] + correction2[k,l]
#                 im_points_cor3[j,k,l] = im_points_night[j,k,l] * mask[k,l]
#             else :
#                 im_points_cor1[j,k,l] = None
#                 im_points_cor2[j,k,l] = None
#                 im_points_cor3[j,k,l] = None

# MAX_NO_COR = [[],[],[]]
# MAX_COR1 = [[],[],[]]
# MAX_COR2 = [[],[],[]]
# MAX_COR3 = [[],[],[]]

# VAL_NO_COR = [[],[],[]]
# VAL_COR1 = [[],[],[]]
# VAL_COR2 = [[],[],[]]
# VAL_COR3 = [[],[],[]]
# for i in tqdm(range(len(DataWindow))):
#         MAX_NO_COR[i] = np.nanmax(im_points_night[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:],axis=(1,2))
#         MAX_COR1[i] = np.nanmax(im_points_cor1[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:],axis=(1,2))
#         MAX_COR2[i] = np.nanmax(im_points_cor2[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:],axis=(1,2))
#         MAX_COR3[i] = np.nanmax(im_points_cor3[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:],axis=(1,2))

#         VAL_NO_COR[i] = im_points_night[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:].ravel()
#         VAL_COR1[i] = im_points_cor1[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:].ravel()
#         VAL_COR2[i] = im_points_cor2[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:].ravel()
#         VAL_COR3[i] = im_points_cor3[df1a.iloc[:]['WDWInputDataWindow'] == DataWindow[i],:,:].ravel()

#         # removing nan values
#         MAX_NO_COR[i] = MAX_NO_COR[i][~np.isnan(MAX_NO_COR[i])]
#         MAX_COR1[i] = MAX_COR1[i][~np.isnan(MAX_COR1[i])]
#         MAX_COR2[i] = MAX_COR2[i][~np.isnan(MAX_COR2[i])]
#         MAX_COR3[i] = MAX_COR3[i][~np.isnan(MAX_COR3[i])]
        
#         VAL_NO_COR[i] = VAL_NO_COR[i][~np.isnan(VAL_NO_COR[i])]
#         VAL_COR1[i] = VAL_COR1[i][~np.isnan(VAL_COR1[i])]
#         VAL_COR2[i] = VAL_COR2[i][~np.isnan(VAL_COR2[i])]
#         VAL_COR3[i] = VAL_COR3[i][~np.isnan(VAL_COR3[i])]

#         # computing number of pixels that could have been counted in the previous window
#         WIN_VAL = [2**13,2**14,2**15,2**16]
#         print('#################')
#         print(f"Mode {DataWindow[i]} :")
#         for lim in WIN_VAL:
#             if len(VAL_NO_COR[i]) > 0:
#                 print(f"{(sum(VAL_NO_COR[i]<lim)/len(VAL_NO_COR[i]))*100.:2f} % of pixels below {lim}")
                    
#         for lim in WIN_VAL:
#             if len(VAL_COR3[i]) > 0:
#                 print(f"{(sum(VAL_COR3[i]<lim)/len(VAL_COR3[i]))*100.:2f} % of corrected pixels below {lim}")


        


# xlabels = ["no correction","correction 1","correction 2","correction 3"]
# for i in range(len(DataWindow)):
#     win = DataWindow[i] 
#     DATA = [VAL_NO_COR[i],VAL_COR1[i],VAL_COR2[i],VAL_COR3[i]]   
#     plt.close(f"Window mode {win} pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.figure(f"Window mode {win} pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.boxplot(DATA,whis=[0,100],showfliers=False)
#     plt.title(f"Window mode {win} pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.xlabel('Correction')
#     plt.xticks([1,2,3,4],xlabels)
#     plt.ylabel('Pixel value')
#     plt.show()

#     DATA = [MAX_NO_COR[i],MAX_COR1[i],MAX_COR2[i],MAX_COR3[i]]   
#     plt.close(f"Window mode {win} MAX pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.figure(f"Window mode {win} MAX pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.boxplot(DATA,whis=[0,100],showfliers=False)
#     plt.title(f"Window mode {win} maximum pixel value distribution (SZA > {ang_lim:.1f} deg)")
#     plt.xlabel('Correction')
#     plt.xticks([1,2,3,4],xlabels)
#     plt.ylabel('Maximum pixel value')
#     plt.show()

        
        
            







# %%

ang_min = 100
ang_max = ang_min + 1
step = 3
x_art = 17
y_art = 10
x_dark = x_art
y_dark = y_art - 8

# X = im_points_night[step:-1,y_dark,x_dark]
# Y = im_points_night[0:-(step+1),y_art,x_art]
# C = sza_points_night[0:-(step+1),y_art,x_art]

X = im_points[step:-1,y_dark,x_dark]
Y = im_points[0:-(step+1),y_art,x_art]
C = sza_points[0:-(step+1),y_art,x_art]



indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['nadir_sza']>ang_min) & (df1a.iloc[:-(step+1)]['nadir_sza']<ang_max)

X = X[indexes]
Y = Y[indexes]
C = C[indexes]
C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']

from scipy.optimize import curve_fit




plt.figure()
plt.scatter(X,Y,c=C)

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

plt.plot(X,intercept + slope*X,color='red')
plt.title(f"Artifact correlation, offset of {step} images. Slope = {slope:.3f}, intercept = {intercept:.1f}, R**2 = {rsquare:.3f}")
plt.xlabel(f'Pixel ({x_dark},{y_dark}) (correct value)')
plt.ylabel(f'Pixel ({x_art},{y_art}) (artifact)')
plt.colorbar(label='TPssa')
plt.show()
# %%


# %%
ang_min = 80
ang_max = ang_min +15
im_R = np.ones_like(im_points[0,:,:])
im_bias = np.zeros_like(im_points[0,:,:])
im_slope = np.ones_like(im_points[0,:,:])

for x_art in range(0,56):
    for y_art in range(3,14):

        step = int(y_art//2.6)
        x_dark = x_art
        y_dark = int(round(y_art%2.6))

        # X = im_points_night[step:-1,y_dark,x_dark]
        # Y = im_points_night[0:-(step+1),y_art,x_art]
        # C = sza_points_night[0:-(step+1),y_art,x_art]

        X = im_points[step:-1,y_dark,x_dark]
        Y = im_points[0:-(step+1),y_art,x_art]
        C = sza_points[0:-(step+1),y_art,x_art]

        indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['TPssa']>ang_min) & (df1a.iloc[:-(step+1)]['TPssa']<ang_max)

        X = X[indexes]
        Y = Y[indexes]
        C = C[indexes]
        C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']

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


plt.title(f'Bias value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.imshow(im_bias[:,:],origin='lower')
plt.colorbar()
plt.show()

plt.title(f'Slope value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.imshow(im_slope[:,:],origin='lower')
plt.colorbar()
plt.show()

plt.title(f'R squared value ({ang_min} deg < TPssa < {ang_max} deg)')
plt.imshow(im_R[:,:],origin='lower')
plt.colorbar()
plt.show()

# %%
ang_min = 80
ang_max = 95
ang_step = 1
ANG = range(ang_min,ang_max,ang_step)
BIAS = []
R2 = []
step = 3
x_art = 17
y_art = 9
x_dark = x_art
y_dark = y_art - 8
for ang in ANG:
    X = im_points[step:-1,y_dark,x_dark]
    Y = im_points[0:-(step+1),y_art,x_art]
    C = sza_points[0:-(step+1),y_art,x_art]

    indexes = (~np.isnan(X)) & (~np.isnan(Y)) & (df1a.iloc[:-(step+1)]['TPssa']>ang) & (df1a.iloc[:-(step+1)]['TPssa']<ang+ang_step)
    X = X[indexes]
    Y = Y[indexes]
    C = C[indexes]
    C = df1a.iloc[np.append(indexes,(step+1)*[False])]['TPssa']

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

    BIAS.append(intercept)
    R2.append(rsquare)


plt.title(f'Bias value (Pixel ({x_art},{y_art}))')
plt.plot(ANG,BIAS)
plt.xlabel('TPSSA angle (deg)')
plt.ylabel('Bias value')
plt.show()

plt.title(f'R squared value (Pixel ({x_art},{y_art}))')
plt.plot(ANG,R2)
plt.xlabel('TPSSA angle (deg)')
plt.ylabel('R squared value')
plt.show()
# %%
