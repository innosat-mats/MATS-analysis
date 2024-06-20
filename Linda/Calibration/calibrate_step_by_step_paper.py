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
from mats_utils.rawdata.calibration import calibrate_dataframe



def plot_calib_step(dfentry, step1name,step2name,title,divide=False, clim1=None, clim2=None,clim3=None,fig=None, ax=None):

    step1=dfentry[step1name]
    step2=dfentry[step2name]

    if fig is None:
        fig, ax = plt.subplots(3,1)
        fig.suptitle(title+' '+dfentry['channel'], fontsize=16)
    if clim1 is None:
        clim1=[np.mean(step1)-np.std(step1), np.mean(step1)+np.std(step1)]
    if clim2 is None:
        clim2=[np.mean(step2)-np.std(step2), np.mean(step2)+np.std(step2)]   
    if divide:
        change=step2/step1
    else:
        change=step2-step1

    plot_CCDimage(step1, fig=fig, axis=ax[0], title=step1name, clim=clim1)
    ax[0].text(10, 60, 'max: '+str(np.max(step1)), color='white')
    ax[0].text(10, 80, 'min: '+str(np.min(step1)), color='white')
    plot_CCDimage(step2, fig=fig, axis=ax[1], title=step2name, clim=clim2)
    ax[1].text(10, 60, 'max: '+str(np.max(step2)), color='white')
    ax[1].text(10, 80, 'min: '+str(np.min(step2)), color='white')

    plot_CCDimage(change, fig=fig, axis=ax[2], title='change')
    ax[2].text(10, 60, 'max: '+str(np.max(change)), color='white')
    ax[2].text(10, 80, 'min: '+str(np.min(change)), color='white')

    plt.tight_layout()
    plt.savefig('../output/'+ title +'.png')
    plt.show()





    return



# # #%% Select on explicit time NLC
start_time = DT.datetime(2023, 2, 2, 19, 38, 0) 
stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)

# Dayglow
#start_time = DT.datetime(2023, 2, 20, 19, 20, 0)
#stop_time = DT.datetime(2023, 2, 20, 19, 21, 0)
#Nightglow
#start_time = DT.datetime(2023, 2, 20, 20, 48, 0)
#stop_time = DT.datetime(2023, 2, 20, 20, 49, 0)

df = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)

#%%
n=10
uv1=df[df.channel=='UV1'][:n].reset_index()
uv2=df[df.channel=='UV2'][:n].reset_index()
ir1=df[df.channel=='IR1'][:n].reset_index()
ir2=df[df.channel=='IR2'][:n].reset_index()
ir3=df[df.channel=='IR3'][:n].reset_index()
ir4=df[df.channel=='IR4'][:n].reset_index()


# pickle.dump(CCDitemsuv1, open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'wb'))
# pickle.dump(CCDitemsuv2, open('testdata/CCD_items_in_orbit_NLCuv2.pkl', 'wb'))

# #%%
# #with open('testdata/CCD_items_in_orbit_UVIR.pkl', 'rb') as f:
# with open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'rb') as f:
#     CCDitems = pickle.load(f)
#%%
#instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

uv1cal=calibrate_dataframe(uv1, instrument, debug_outputs=True)
uv2cal=calibrate_dataframe(uv2, instrument, debug_outputs=True)
ir1cal=calibrate_dataframe(ir1, instrument, debug_outputs=True)
ir2cal=calibrate_dataframe(ir2, instrument, debug_outputs=True)
ir3cal=calibrate_dataframe(ir3, instrument, debug_outputs=True)
ir4cal=calibrate_dataframe(ir4, instrument, debug_outputs=True)

#%%
# plot all the different calibration steps
#
dfentry=ir3cal.iloc[0]

plot_calib_step(dfentry, 'image_lsb','image_se_corrected','SE correction')
plot_calib_step(dfentry, 'image_se_corrected','image_hot_pixel_corrected','hot-pixel correction')
plot_calib_step(dfentry, 'image_hot_pixel_corrected','image_bias_sub','bias subtraction')
plot_calib_step(dfentry, 'image_bias_sub','image_linear','linearization',divide=True)
plot_calib_step(dfentry, 'image_linear','image_desmeared','desmear')
plot_calib_step(dfentry, 'image_desmeared','image_dark_sub','dark subtraction')
plot_calib_step(dfentry, 'image_dark_sub','image_flatfielded','flatfielding',divide=True)
plot_calib_step(dfentry, 'image_flatfielded','image_flipped','flipping')





# %%

substeps=['image_lsb',
 #'image_se_corrected',
 #'image_hot_pixel_corrected',
 'image_bias_sub',
 'image_linear',
 'image_desmeared',
 #'image_dark_sub',
 'image_flatfielded',
# 'image_flipped',
 'image_calibrated']

nmax=len(substeps)-1
stepnames=['a) '+ dfentry.channel+' '+str(dfentry.TMHeaderTime)[0:10]+' Level 1a [Counts]', 
          #'SE corrected', 
          #'hot-pixel corrected', 
           'c) Bias, hotpixel and single events subtracted [Counts]', 
           'e) Linearized [Counts]', 
           'g) Desmeared [Counts]', 
           #'Dark current subtracted [Counts]', 
           'i) Calibrated [$10^{12}$ photons/nm/s/m$^{2}$/str]', 
          # 'flipped', 
           'calibrated']

processnames=[
                # 'SE correction',
                # 'Hot-pixel correction',
                'b) Bias, hot pixel and single event subtraction [Counts]',
                'd) Linearization factor [Unitless]',
                'f) Desmearing [counts]',
                #'Dark current subtraction [Counts]',
                'h) Calib. factor [$10^{12}$ photons/nm/s/m$^{2}$/str/counts]',
                # 'Flipping',
                'Calibration']


#All substeps below:

substeps=['image_lsb',
 'image_se_corrected',
 'image_hot_pixel_corrected',
 'image_bias_sub',
 'image_linear',
 'image_desmeared',
 'image_dark_sub',
 'image_flatfielded',
 'image_flipped',
 'image_calibrated']

nmax=len(substeps)-1
stepnames=[dfentry.channel+' '+str(dfentry.TMHeaderTime)[0:10]+' Level 1a [Counts]', 
          'SE corrected', 
          'hot-pixel corrected', 
           'Bias, hotpixel and single events subtracted [Counts]', 
           'Linearized [Counts]', 
           'Desmeared [Counts]', 
           'Dark current subtracted [Counts]', 
           ' Calibrated [$10^{12}$ photons/nm/s/m$^{2}$/str]', 
           'flipped', 
           'calibrated']

processnames=[
                'SE correction',
                'Hot-pixel correction',
                'Bias, hot pixel and single event subtraction [Counts]',
                'Linearization factor [Unitless]',
                'Desmearing [counts]',
                'Dark current subtraction [Counts]',
                'Calib. factor [$10^{12}$ photons/nm/s/m$^{2}$/str/counts]',
                'Flipping',
                'Calibration']




fig, ax = plt.subplots(nmax, 1, figsize=(10,18))
for i in range(nmax):
    plot_CCDimage(dfentry[substeps[i]], fig=fig, axis=ax[i], title=stepnames[i])
    ax[i].text(5, 60,'min: '+str(np.min(dfentry[substeps[i]])), color='white')
    ax[i].text(5, 90,'max: '+str(np.max(dfentry[substeps[i]])), color='white')
plt.tight_layout()







# %%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create a figure with a custom layout using GridSpec
fig = plt.figure(figsize=(10, 9))
gs = gridspec.GridSpec(nmax*2, 2)

# Display the images in the left column and their differences in the right column
for i in range(nmax-1):
    # Display dfentry[substeps[i]] in the left column
    ax1 = fig.add_subplot(gs[2*i:2*i+2, 0])
    sc1=ax1.imshow(dfentry[substeps[i]], origin='lower', interpolation="none")
    cbar = fig.colorbar(sc1, ax=ax1, orientation='vertical')
    ax1.set_title(stepnames[i])
    ax1.set_aspect('auto')

    # Display dfentry[substeps[i]] - dfentry[substeps[i+1]] in the right column
    ax2 = fig.add_subplot(gs[2*i+1:2*i+3, 1])
    if substeps[i+1] in ['image_linear', 'image_flatfielded']:
        field=dfentry[substeps[i+1]]/dfentry[substeps[i]]
    else:
        field=dfentry[substeps[i]]-dfentry[substeps[i+1]]

    if substeps[i+1] in ['image_flatfielded']:
        clim=[0.050, 0.060]
        sc2 = ax2.imshow(field, origin='lower', interpolation="none", clim=clim)
    else:
        sc2 = ax2.imshow(field, origin='lower', interpolation="none")
    ax2.text(3, 40, 'max: '+str(np.max(field)), color='white')
    ax2.text(3, 60, 'min: '+str(np.min(field)), color='white')

    ax2.set_title(processnames[i])    
    ax2.set_aspect('auto')
    cbar = fig.colorbar(sc2, ax=ax2, orientation='vertical')

# Display dfentry[substeps[8]] in the last row of the left column
ax3 = fig.add_subplot(gs[-2:, 0])
clim=[np.mean(dfentry[substeps[nmax-1]])-2*np.std(dfentry[substeps[nmax-1]]), np.mean(dfentry[substeps[nmax-1]])+2*np.std(dfentry[substeps[nmax-1]])]
sc3 = ax3.imshow(dfentry[substeps[nmax-1]], origin='lower', interpolation="none", clim=clim)
ax3.set_title(stepnames[nmax-1])
ax3.set_aspect('auto')
cbar = fig.colorbar(sc3, ax=ax3, orientation='vertical')

# Adjust the space between the subplots
plt.subplots_adjust(hspace=0.5)
plt.tight_layout()
# Display the figure
plt.show()

#fig.savefig('../output/calibration_steps.png')

# %%
#Error Analysis
import mats_utils.error_estimate.error_estimate as ee
CCDitems = dataframe_to_ccd_items(ir1)

