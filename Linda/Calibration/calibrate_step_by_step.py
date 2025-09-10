#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_l1_processing.read_in_functions import read_CCDitems 
import numpy as np
import pickle
import matplotlib.pyplot as plt
from database_generation.experimental_utils import plot_CCDimage



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
        fig, ax = plt.subplots(1,3)
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


    plt.tight_layout()
    plt.savefig('../output/'+ title +'.png')
    plt.show()

    #plot_CCDimage(change,title='change'+title, borders=True, nrsig=3)

    #print('mean change: ', np.mean(change))
    #print('max change: ', np.max(change))
    #print('min change: ', np.min(change))



    return




# #%% Select on explicit time
start_time = DT.datetime(2023, 5, 5, 20, 10)
stop_time = DT.datetime(2023, 5, 5, 20, 15)



# # #%% Select on explicit time NLC
start_time = DT.datetime(2023, 2, 9, 18, 54, 39) 
stop_time = DT.datetime(2023, 2, 9, 19, 3, 44)

# # #%% Select on explicit time nadir
start_time = DT.datetime(2023, 4, 6, 0, 6, 0) 
stop_time = DT.datetime(2023, 4, 6, 0, 7, 0)


#start_time = DT.datetime(2023, 2, 12, 1, 10)
#stop_time = DT.datetime(2023, 2, 12, 1, 15)
df = read_MATS_data(start_time,stop_time,version='1.0',level='1a',dev=False)

#%%
#select only nadir
dfnadir = df[df['channel']=='NADIR'] 
CCDitems = dataframe_to_ccd_items(dfnadir)


#%%

# # #%% Select on explicit time NLC

# # # #%% Select on explicit time NLC
# start_time = DT.datetime(2023, 2, 3, 0, 27)
# stop_time = DT.datetime(2023, 2, 3, 0, 31)

# df = read_MATS_data(start_time,stop_time,version='0.6',level='1a',dev=False)
# uv1=df[df.channel=='UV1']
# uv2=df[df.channel=='UV2']
# ir1=df[df.channel=='IR1']
# ir2=df[df.channel=='IR2']
# ir3=df[df.channel=='IR3']
# ir4=df[df.channel=='IR4']


# CCDitemsuv1 = dataframe_to_ccd_items(uv1)
# CCDitemsuv2 = dataframe_to_ccd_items(uv2)

# pickle.dump(CCDitemsuv1, open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'wb'))
# pickle.dump(CCDitemsuv2, open('testdata/CCD_items_in_orbit_NLCuv2.pkl', 'wb'))

# #%%
# #with open('testdata/CCD_items_in_orbit_UVIR.pkl', 'rb') as f:
# with open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'rb') as f:
#     CCDitems = pickle.load(f)
#%%
instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')


#%%
#IR1§§§§§
i = 5
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)
 #%%
#IR2
i = 2
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)

#%%
#IR3
i = 1
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)
#%%
#IR4
i = 6
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)
#%%
#UV1
i = 3
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)

#UV2
i = 4
images = L1_calibrate(CCDitems[i], instrument,return_steps=True)
datechannel = str(CCDitems[i]["TMHeaderTime"])[0:10] + '_' + CCDitems[i]["channel"]
make_calibration_plots(images,datechannel, savefig=True)


#%%
for i in range(0,7):
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
    #plot_calib_step(image_desmeared,image_dark_sub,'dark_subtraction' + datechannel,errors)
    #plot_calib_step(image_dark_sub,image_flatfielded,'flatfielding' + datechannel,errors,divide=True)
    plot_calib_step(image_se_corrected,image_hot_pixel_corrected,'hot-pixel_correction' + datechannel,errors)



#%%
