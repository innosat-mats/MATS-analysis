# %% [markdown]
# This code will generate 3x1 monthly image means for IR1

# %%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
import argparse
from datetime import date, timedelta
from mats_utils.plotting.plotCCD import all_channels_plot
#from mats_utils.daily_preview.temp_nadirs import NADIR_geolocation, average_stacking
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
import matplotlib
import xarray as xr
import scipy
import proplot as pplt
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
from database_generation.experimental_utils import plot_CCDimage
def stack(image):
    stacked_im=np.stack(image.values)
    return stacked_im

def generate_movies(df, data_folder, outputfolder,movienameprefix='movie', field_of_choise='ImageCalibrated'):
    """
    Generate a movie from a dataframe with a certain channel

    Parameters
    ----------
    df : pandas dataframe
        Dataframe with the data
    channel : str
        Channel to be used
    data_folder : str
        Folder where the images are stored
    outputfolder : str
        Folder where the movie is to be stored
    moviename : str, optional
        Name of the movie. The default is 'movie.gif'.
    field_of_choise : str, optional
        Field to be used. The default is 'ImageCalibrated'.

    """
    from mats_utils.plotting.plotCCD import orbit_plot
    from mats_utils.plotting.animate import generate_gif
    import os
    

    imagedir=data_folder+'/movies/'
    # Create directory if it does not exist
    if not os.path.exists(imagedir):
        os.makedirs(imagedir)

    orbit_plot(df,imagedir,nbins=7, cmap='YlGnBu_r', plothistogram=False, field_of_choise='ImageCalibrated', useplotCCDimage=True)
    
    for CCDno in range(0, 8):
        CCDs = df[df['CCDSEL'] == CCDno]
        if len(CCDs) > 0:
            imagedir_CCDno=f"{imagedir}CCDSEL{str(CCDno)}"   
    
            generate_gif(imagedir_CCDno, outputfolder+'/'+movienameprefix+'_CCDSEL'+str(CCDno)+'.gif')

    print('Gifs are to be found in ',outputfolder+'/'+movienameprefix+'_CCDSEL*.gif')
    



# %%
df=pd.read_pickle(f'/Users/lindamegner/MATS/MATS-retrieval/data/BjornsData/022023_2.pk1')



# %%
# ASCENDING MAP PLOT /// UPDATED TO READ FROM L1b_v05

fig = pplt.figure(figwidth='12cm')
axs = fig.subplots(ncols=1, nrows=3, proj='cyl')
fig.format(suptitle='Dayglow observations from IR1 (~17:30 LT)',abc='a.',abcborder=False)
axs.format(coast=True,landzorder=4,latlines=30, lonlines=60)
#axs.format(lonlabels='b')
#axs[0].format(latlabels='l')
#axs[3].format(latlabels='l')


#monthstrings=["022023_2","032023_1","032023_2","042023_1"]
monthstrings=["022023_1","032023_1","042023_1","022023_2","032023_2",False]
monthstrings=["022023","032023","042023"]
monthstrings=["022023"]

#monthstrings=["022023_1","032023_1",False,False,False,False]
#monthstrings=[False,"032023_1",False,False,False,False]
# dayglow
dmin,dmax = 0, 95
ascending = True

def stack(image):
    stacked_im=np.stack(image.values)
    return stacked_im

#subtitles=['Feb (1st-14th)']#, 'Feb (15st-28th)','March (1st-14th)', 'March (15st-31st)','April (1st-14th)', 'April (15st-31st)']
#for monthstring in monthstrings:
monthstring="022023"
week=2
#    for offset in [1]:
#        for j in [1]:
#            week = j + offset
#            print(f'Processing {monthstring}_{week}')
df=pd.read_pickle(f'/Users/lindamegner/MATS/MATS-retrieval/data/BjornsData/{monthstring}_{week}.pk1')
#df=df[:100]

#df=pd.read_pickle(f'/media/waves/AVAGO/data/MATS/pandas_csv/L1b_v05/7070/CCDSEL1/alts_80to90/finished/{monthstring}_{week}.pk1')

midday=DT.time(12, 0, 0)

df = df[(df['CCDSEL'] == 1)]

# ascending
if ascending:
    df = df[(df['TPlocaltime'].dt.time > midday)]
else:
    df = df[(df['TPlocaltime'].dt.time < midday)]

heights=np.array([80,81,82,83,84,85,86,87,88,89,90])
df.ImageCalibrated = df.ImageCalibrated.div(df.TEXPMS/1000)
#df.ImageCalibrated= df.ImageCalibrated.apply(np.mean, axis=1)

#add_field_with_subtracted_rolling_mean(df, 'ImageCalibrated', 'ImageCalibrated_minus_rolling_mean', window_before=20, window_after=20, skipbefore=5, skipafter=5)
#plot_CCDimage(df['ImageCalibrated_minus_rolling_mean'].iloc[0])
#df.ImageCalibrated_minus_rolling_mean_std= df.ImageCalibrated_minus_rolling_mean.apply(np.std, axis=1)

# Run all ImageCalibrated through median filter
#df.ImageCalibratedFiltered = df.ImageCalibrated.apply(lambda x: scipy.ndimage.median_filter(x, size=(3, 3)))


stacked_im=stack(df.ImageCalibrated)
#stacked_im_minus_rolling_mean=stack(df.ImageCalibrated_minus_rolling_mean)
CCDs_allheights =xr.Dataset(data_vars=dict(
            im=(["time", "height","xpix"], stacked_im[:,:,:]),
            lon=(["time"], df["TPlon"]),
            lat=(["time"], df["TPlat"])),
coords=dict(
    time=pd.to_datetime(df["EXPDate"]),
    height=heights),

attrs=dict(description="CCD data."),
)


# Calculate the rolling mean with a window size of your choice (e.g., 3)
CCDs_allheights['rolling_mean'] = CCDs_allheights.im.rolling(time=9, center=True, min_periods=1).mean()
# Subtract the rolling mean from the original mean_im
CCDs_allheights['im_detrended'] = CCDs_allheights.im - CCDs_allheights.rolling_mean
#Apply median filter to detrended image
from scipy.signal import medfilt
CCDs_allheights['im_detrended_filtered'] = (["time", "height","xpix"], medfilt(CCDs_allheights.im_detrended, kernel_size=[1,1,7]))

CCDs_allheights['std_horizontal']=CCDs_allheights['im_detrended_filtered'].std(dim='xpix').mean(dim='height')
CCDs_allheights['std_vertical']=CCDs_allheights['im_detrended_filtered'].std(dim='height').mean(dim='xpix')




CCDs=CCDs_allheights.sel(height=86).where(CCDs_allheights.sel(height=86))

m=axs[0].scatter(CCDs.lon.values,CCDs.lat.values,c=CCDs.rolling_mean.values,s=0.3,cmap='Thermal',vmin=50,vmax=130)
axs[0].format(title='Airglow at 86 km', titleloc='l')
m=axs[1].scatter(CCDs.lon.values,CCDs.lat.values,c=CCDs.std_horizontal.values,s=0.5,cmap='Thermal',vmin=0,vmax=1)
axs[1].format(title='horizontal variability', titleloc='l')
m=axs[2].scatter(CCDs.lon.values,CCDs.lat.values,c=CCDs.std_vertical.values,s=0.5,cmap='Thermal',vmin=0,vmax=1)    
axs[2].format(title='vertical variability', titleloc='l')
        
        
#axs[1].colorbar(m)
fig.colorbar(m, label="Radiance [units]",loc='r',length=0.8)




#fig.save('/home/waves/projects/MATS/global_patterns/notebooks/paper_figures/output/test.png',format='png')

# %%

#convert CCDs_allheights to dataframe
df_new=CCDs_allheights.to_dataframe()
generate_movies(df_new[10:20], '/Users/lindamegner/MATS/MATS-retrieval/data/temp', '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output', movienameprefix='mymovie',field_of_choise='ImageCalibrated')

#%%
# Assuming CCDs is already defined as in your example
# Calculate the rolling mean with a window size of your choice (e.g., 3)
rolling_mean = CCDs_allheights.im.rolling(time=3, center=True).mean()
# Subtract the rolling mean from the original mean_im
CCDs_allheights['im_detrended'] = CCDs_allheights.im - rolling_mean
CCDs_allheights['im_detrended'].shape
plot_CCDimage(CCDs_allheights['im_detrended'][40,:,:])
hej=CCDs_allheights['im_detrended'][40,:,:].std(dim='xpix', dim='height')
