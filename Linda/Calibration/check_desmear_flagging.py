#%%
import pandas as pd

from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate, calibrate_all_items
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_l1_processing.read_in_functions import read_CCDitems 
import numpy as np
import pickle
import matplotlib.pyplot as plt
from database_generation.experimental_utils import plot_CCDimage

#%%

def make_calibration_plots(images,datechannel, savefig=False):

    print(len(images))
    (image_lsb,
     image_se_corrected, 
     image_hot_pixel_corrected, 
     image_bias_sub, 
     image_linear,
     image_desmeared, 
     image_dark_sub, 
     image_flatfielded, 
     image_flipped, 
     image_calibrated, 
     errors,
    ) = images


    plot_calib_step(image_lsb,image_se_corrected,'SE_correction_' + datechannel,errors)
    plot_calib_step(image_se_corrected,image_hot_pixel_corrected,'hot-pixel_correction' + datechannel,errors)
    plot_calib_step(image_hot_pixel_corrected,image_bias_sub,'bias_subtraction' + datechannel,errors)
    plot_calib_step(image_bias_sub,image_linear,'linearization' + datechannel,errors,divide=True)#, clim3=[0.99,1.002])
    change=image_linear/image_bias_sub
    print('datechannel', datechannel)
    print('mean change: ', change)
    plot_calib_step(image_linear,image_desmeared,'desmear' + datechannel,errors)
    plot_calib_step(image_desmeared,image_dark_sub,'dark_subtraction' + datechannel,errors)
    plot_calib_step(image_dark_sub,image_flatfielded,'flatfielding' + datechannel,errors,divide=True)
    plot_calib_step(image_flatfielded,image_flipped,'flipping' + datechannel,errors)
    #plot_calib_step(image_flipped,image_calibrated,'image_calibrated',errors)

    return

def plot_calib_step(step1,step2,title,error,divide=False, clim1=None, clim2=None,clim3=None,fig=None, ax=None):

    error = np.zeros(np.shape(step1))

    if fig is None:
        fig, ax = plt.subplots(3,1)
        fig.suptitle(title, fontsize=16)
        ax0 = ax[0]
        ax1 = ax[1]
        ax2 = ax[2]
    if clim1 is None:
        clim1=[np.mean(step1)-np.std(step1), np.mean(step1)+np.std(step1)]
    sc = ax0.imshow(step1, clim=clim1,origin='lower')
    cbar = fig.colorbar(sc, ax=ax0, orientation='vertical')
    if clim2 is None:
        clim2=[np.mean(step2)-np.std(step2), np.mean(step2)+np.std(step2)]   
    sc = ax1.imshow(step2, clim=clim2,origin='lower')
    cbar = fig.colorbar(sc, ax=ax1, orientation='vertical')
    if divide:
        change=step2/step1
        if clim3 is None:
            clim3=[np.mean(change)-np.std(change), np.mean(change)+np.std(change)]
        sc = ax2.imshow(change,clim=clim3,origin='lower')
    else:
        change=step2-step1
        sc = ax2.imshow(change,origin='lower')
    cbar = fig.colorbar(sc, ax=ax2, orientation='vertical')
    # axs[3].imshow(error,origin='lower')
    # cbar = fig.colorbar(sc, ax=axs[3], orientation='vertical')
    #set aspect of the plots to be 4:1
    ax0.set_aspect(1/12)
    ax1.set_aspect(1/12)
    ax2.set_aspect(1/12)   


    plt.tight_layout()
    plt.savefig('../output/'+ title +'.png')
    plt.show()

    #plot_CCDimage(change,title='change'+title, borders=True, nrsig=3)

    #print('mean change: ', np.mean(change))
    #print('max change: ', np.max(change))
    #print('min change: ', np.min(change))



    return



#%%
# #%% Select on explicit time n 2023 03 22 09:15:00 UTC - 2023 03 22 09:25:00 UTC
#start_time = DT.datetime(2023, 3, 22, 9, 15)
#stop_time = DT.datetime(2023, 3, 22, 9, 25)
start_time = DT.datetime(2023, 4, 22, 0, 7, 5)
stop_time = DT.datetime(2023, 4, 22, 0, 7, 16)


df = read_MATS_data(start_time,stop_time,version='1.0',level='1a',dev=False)

#%%
pickle_filename = f"../output/df_{start_time.strftime('%Y%m%d_%H%M%S')}_to_{stop_time.strftime('%Y%m%d_%H%M%S')}.pkl"
df.to_pickle(pickle_filename) 

#%%
#load df from pickle
df = pd.read_pickle(pickle_filename)


df_IR2 = df[df['channel']=='IR2'] 



#%%
instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')


# #IR1§§§§§
# i = 5
# images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
# datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
# make_calibration_plots(images,datechannel, savefig=True)
 #%%


#2023-03-22 09:15:00.529052734+0000
#target_time = pd.Timestamp('2023-03-22 09:15:00.529052734+0000', tz='UTC')
#2023-04-22 00:07:14.378448+00:00
#target_time = pd.Timestamp('2023-03-22 09:20:12.528030396+0000', tz='UTC')
target_time = pd.Timestamp('2023-04-22 00:07:14.378448+00:00', tz='UTC')
exp_dates = pd.to_datetime(df_IR2.EXPDate, utc=True)
nearest_idx = (exp_dates - target_time).abs().idxmin()
df_IR2_selected = df_IR2.loc[[nearest_idx]]

CCDitems = dataframe_to_ccd_items(df_IR2_selected)
#Print timestamp of the selected item
print('Selected timestamp: ', df_IR2_selected.EXPDate)
# print tplat and tplon of the selected item
print('Selected tplat: ', df_IR2_selected.TPlat)
print('Selected tplon: ', df_IR2_selected.TPlon)


#%%
for i in range(0,1):
    images= L1_calibrate(CCDitems[i], instrument,return_steps=True)
    (image_lsb,
        image_se_corrected, 
        image_hot_pixel_corrected, 
        image_bias_sub, 
        image_linear,
        image_desmeared, 
        image_dark_sub, 
        image_flatfielded, 
        image_flipped, 
        image_calibrated, 
        errors,
    ) = images
    datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
    plot_calib_step(image_lsb,image_se_corrected,'SE_correction_' + datechannel,errors)
    plot_calib_step(image_se_corrected,image_hot_pixel_corrected,'hot-pixel_correction' + datechannel,errors)
    plot_calib_step(image_hot_pixel_corrected,image_bias_sub,'bias_subtraction' + datechannel,errors)
    plot_calib_step(image_bias_sub,image_linear,'linearization' + datechannel,errors,divide=True)
    plot_calib_step(image_linear,image_desmeared,'desmear' + datechannel,errors)
    plot_calib_step(image_desmeared,image_dark_sub,'dark_subtraction' + datechannel,errors)
    plot_calib_step(image_dark_sub,image_flatfielded,'flatfielding' + datechannel,errors,divide=True)
    #Print timestamps of each step

 

#%%
calibrate_all_items(CCDitems, instrument, plot=False)

# %%
