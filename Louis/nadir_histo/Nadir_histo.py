# script for studying the number of saturated pixels in NADIR images as a function of solare zenith angle 

#%% Import modules
%matplotlib qt5
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

start_time = DT.datetime(2022, 12, 20, 20, 30, 0)
stop_time = DT.datetime(2022, 12, 20, 22, 0, 0)

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}


#%% reading mesurements
df1a_tot= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
print(len(df1a_tot))
df1a = df1a_tot

# displaying keys
pd.set_option('display.max_rows', 100)
df1a.dtypes


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

#%% computing the saturation level of an array of images
def saturation_level_imagecube(imagecube, sat_val = 32880):
    """
    Function to get the proportion of saturated pixels in an imagecube

    Parameters
    ----------
    imagecube : np.array
        array of 2D images (axes order doesn't matter)
    sat_val : int
        value of a saturated pixel
        
    Returns
    -------
    sat_prop : float
        proportion (between 0 and 1) of saturated pixels
    """
    a,b,c = np.shape(imagecube)
    saturated_im = sat_val*np.ones_like(imagecube) # completely saturated image
    SAT_LEVEL = imagecube == saturated_im # boolean array : True value for saturated pixels, False otherwise
    sat_prop = SAT_LEVEL.sum()/(a*b*c) # proportion of saturated pixels
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

# plotting the histogram of saturation levels
plt.figure()
plt.hist(SAT_LEV)
plt.title(f"Proportion of saturated pixels in each NADIR image ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
plt.xlabel('Proportion of saturated pixels')
plt.ylabel('Number of images')

# saturation level as a function of SZA
plt.figure()
plt.plot(NADIR_SZAS,SAT_LEV,linestyle=' ',marker='.')
plt.title(f"Proportion of saturated pixels in each NADIR image ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
plt.ylabel('Proportion of saturated pixels')
plt.xlabel('Zenith angle (deg)')

plt.show()

#%% computing solar zenith angles for each pixel
df1a = df1a_tot[::10]
n = len(df1a)
a,b = np.shape(df1a.iloc[0]['IMAGE'])
lat_points = np.zeros((n,a,b))
lon_points = np.zeros((n,a,b))
sza_points = np.zeros((n,a,b))
im_points = np.zeros((n,a,b))

for i in tqdm(range(n)):
    ccditem = df1a.iloc[i]
    im = ccditem['IMAGE']
    lat_map,lon_map,sza_map = NADIR_geolocation(ccditem,x_sample=4,y_sample=4,interp_method='cubic')
    lat_points[i,:,:] = lat_map
    lon_points[i,:,:] = lon_map
    sza_points[i,:,:] = sza_map
    im_points [i,:,:] = im


#%%

nb_bin = 30
sza_bin_lim = np.linspace(np.min(sza_points),np.max(sza_points),nb_bin+1)
step = sza_bin_lim[1]-sza_bin_lim[0]
sza_bin = np.linspace(sza_bin_lim[0]+0.5*step,sza_bin_lim[-1]-0.5*step,nb_bin)
VALUES = []
NB_SAT = np.zeros(nb_bin)
NB_ZERO = np.zeros(nb_bin)
data = []

for i in tqdm(range(nb_bin)):
    sza_min = sza_bin_lim[i]
    sza_max = sza_bin_lim[i+1]
    SZAS = []
    VALUES = []
    for j in range(n):
        for k in range(a):
            for l in range(b):
                if sza_min <= sza_points[j,k,l] and sza_points[j,k,l] <= sza_max :
                    SZAS.append(sza_points[j,k,l])
                    VALUES.append(im_points[j,k,l])
                    if im_points[j,k,l] == sat_val:
                        NB_SAT[i] += 1
                    
    NB_SAT[i] = NB_SAT[i]/len(VALUES)
    data.append(VALUES)
    # plt.figure()
    # plt.hist(VALUES)
    # plt.title(f"Histogram of pixel values for  {sza_min:.2f}deg < sza < {sza_max:.2f}deg  ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
    # plt.xlabel('Proportion of pixels')
    # plt.ylabel('Pixel values')
    # plt.show()
    

plt.figure()
plt.plot(sza_bin,NB_SAT,linestyle=' ',marker='o')
plt.title(f"Proportion of saturated pixels ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
plt.xlabel('SZA')
plt.ylabel('Proportion of saturated pixels')
plt.show()


xlabels = []
for lab in sza_bin:
    xlabels.append(f"{lab:.1f}")
plt.figure()
plt.boxplot(data,whis=[0,100],showfliers=False)
plt.title(f"Pixel values ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
plt.xlabel('SZA')
plt.xticks(range(1,nb_bin+1),xlabels)
plt.ylabel('Pixel value')
plt.show()



# %%
DataWindow = ['13..2','14..3','15..2']
MAX_VAL = [[],[],[]]
ang_lim = 105
for i in range(len(DataWindow)):
    for j in range(n):
        max_val = 0
        if df1a.iloc[j]['WDWInputDataWindow']:
            for k in range(a):
                for l in range(b):
                    if im_points[j,k,l] > max_val and sza_points[j,k,l] > ang_lim :
                            max_val = im_points[j,k,l]
                            maxk = k
                            maxl = l
            print(maxk,maxl,max_val)
            if max_val != 0:
                MAX_VAL[i].append(max_val)
        

AVG = np.average(im_points,axis=0)

# plt.figure()
# plt.boxplot(MAX_VAL,whis=[0,100],showfliers=False)
# plt.title(f"Max pixel values ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms ; {len(df1a)} images)")
# plt.ylabel('Pixel value')
# plt.show()

xlabels = []
for i in range(len(DataWindow)):
    xlabels.append(f"{DataWindow[i]}|{len(MAX_VAL[i])} images")
plt.figure()
plt.boxplot(MAX_VAL,whis=[0,100],showfliers=False)
plt.title(f"Maximum pixel values ('TEXPMS'={df1a.iloc[0]['TEXPMS']}ms)")
plt.xlabel('Window mode')
plt.xticks([1,2,3],xlabels)
plt.ylabel('Maximum pixel value')
plt.show()

# %%
