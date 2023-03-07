# script for studying the number of saturated pixels in NADIR images as a function of solare zenith angle 

#%% Import modules
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from mats_utils.statistiscs.images_functions import create_imagecube


#%% local variables
# value of a saturated pixel 
sat_val = 32880

# times for start and stop
start_time = DT.datetime(2023, 3, 3, 9, 30, 0)
stop_time = DT.datetime(2023, 3, 3, 11, 30, 0)

# filter selecting Nadir chanel
filter={'CCDSEL': [7,7]}


#%% reading mesurements
df1a= read_MATS_data(start_time, stop_time,filter,level='1a',version='0.5')
print(len(df1a))


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





#%%
NADIR_SZAS = [] # list of nadir solar zenih angles
SAT_LEV = [] # list of saturation level (proportion of saturated pixels)
for i in range(len(df1a)):
    ccditem = df1a.iloc[i]
    nadir_sza, TPsza, TPssa, TPlt = coordinates.angles(ccditem)   
    NADIR_SZAS.append(nadir_sza) 
    SAT_LEV.append(saturation_level(ccditem)) 

print(f"{np.average(SAT_LEV)*100:.1f} % of saturated pixels over the whole image set")
print(f"{np.sum(SAT_LEV==np.ones_like(SAT_LEV))/len(SAT_LEV)*100:.1f} % of completely saturated images over the whole image set")

# plotting NADIR SZA
plt.plot(df1a['EXPDate'],NADIR_SZAS,linestyle=' ',marker='.')
plt.ylabel('NADIR SZA (deg)')

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


#%% filtering out completely saturated images

# a,b = np.shape(df1a.iloc[0]['IMAGE'])
# saturated_im = sat_val*np.ones_like(df1a.iloc[0]['IMAGE']) # completely saturated image

# SATURATED = []
# for i in range(len(df1a)):
#     ccditem = df1a.iloc[i]
#     SATURATED.append(np.array_equal(df1a.iloc[i]['IMAGE'],saturated_im))

# df1a_nonsat = df1a[[not sat for sat in SATURATED]]
# df1a_sat = df1a[SATURATED]
