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

from mats_l1_processing.L1_calibration_functions import meanbin_image_with_BC, bin_image_with_BC, absolute_calibration
import datetime as DT
import toml
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
#from lindas_own_functions import plot_CCDimage_transp

def savefield(field,filename,directory='output', plot=False, mat=False):
    from database_generation.experimental_utils import plot_CCDimage
    from pathlib import Path 


    Path(directory).mkdir(parents=True, exist_ok=True)
    np.savetxt(directory+'/'+filename+'.csv', field)
    np.save(directory+'/'+filename+'.npy', field)
    if mat:
        from scipy.io import savemat
        savemat(directory+'/'+filename+'.mat', {'field': field})
    if plot:
        fig, ax=plt.subplots(1,1)
        plot_CCDimage(field, fig=fig, axis=ax, title=filename)
        fig.savefig(directory+'/'+filename + ".jpg")
    return




def mark_current_cropping(CCDitem, flatfield_scalefield, flatfield_binned):

    fig1, ax1 = plt.subplots(2, 1)
    plot_CCDimage(flatfield_scalefield,fig1, ax1[0], title='Flatfield scalefield'+CCDitem['channel'])
    plot_CCDimage(flatfield_binned,fig1, ax1[1], title='Flatfield binned'+CCDitem['channel'])
 
    ncshift=CCDitem['NCSKIP'] 
    nrshift=CCDitem['NRSKIP']
    ncol=CCDitem['NCOL']*CCDitem['NCBIN CCDColumns']*CCDitem['NCBIN FPGAColumns']
    nrow=CCDitem['NROW']*CCDitem['NRBIN']
    rectangle = plt.Rectangle((ncshift, nrshift), ncol, nrow, facecolor='none', ec='white')
    ax1[0].add_patch(rectangle)


def get_flatfield_comparison_area_parameters():
    # The values below are taken from scale_field in flatfield.py in the database_generation folder
    # # The values of the below should give an area not affected by the baffle and is the area that the 
    # calibration factors from the lab is defined by
    FirstRow = 350 
    LastRow = 400 
    FirstCol =524
    LastCol =1523
    return FirstRow, LastRow, FirstCol, LastCol



def printsettings(CCDitem):
    print("################################################################################")
    print('TMHeaderTime:', CCDitem['TMHeaderTime'])
    print('channel:', CCDitem['channel'])
    print('NCBIN CCDColumns:', CCDitem['NCBIN CCDColumns'])
    print('NCBIN FPGAColumns:', CCDitem['NCBIN FPGAColumns'])
    print('NRBIN:', CCDitem['NRBIN'])
    print('NCSKIP:', CCDitem['NCSKIP'])
    print('NRSKIP:', CCDitem['NRSKIP'])
    print('NCOL:', CCDitem['NCOL'])
    print('NROW:', CCDitem['NROW'])
    print('TEXPMS:', CCDitem['TEXPMS'])
    



#%%
# # # #%% Select on explicit time
start_time = DT.datetime(2023, 4, 11, 6, 10)
stop_time = DT.datetime(2023, 4, 11, 6, 15)
df = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)
CCDitems = dataframe_to_ccd_items(df)
#pickle.dump(CCDitems, open('testdata/CCD_items_in_orbit_l1a.pkl', 'wb'))

# #%% Select on explicit time
start_time = DT.datetime(2023, 5, 5, 20, 10)
stop_time = DT.datetime(2023, 5, 5, 20, 15)
df = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)
CCDitems = dataframe_to_ccd_items(df)
#pickle.dump(CCDitems, open('testdata/CCD_items_in_orbit_UVIR?.pkl', 'wb'))

#%%


#%%

# # #
with open('testdata/CCD_items_in_orbit_UVIR.pkl', 'rb') as f:
# with open('testdata/CCD_items_in_orbit_l1b.pkl', 'rb') as f:
    CCDitems = pickle.load(f)
#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'



#%%


calibration_data = toml.load(calibration_file)
instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

nccditems=7 #len(CCDitems)
#%%
for i in range(0,nccditems):


    CCDitem=CCDitems[i]
    channel = CCDitem['channel']
    #totbin = int(CCDitem["NRBIN"])*int(CCDitem["NCBIN CCDColumns"]) * \
    #    int(CCDitem["NCBIN FPGAColumns"])
    
    


    # The baffle scalefield, where 1 means no effect, ie the flatfield is the 
    # same as the one without baffle. 
    flatfield = np.load(
                calibration_data["flatfield"]["flatfieldfolder"]
                + "flatfield_"
                + channel
                + "_HSM.npy")

    flatfield_binned=meanbin_image_with_BC(CCDitem, flatfield)
    
    # For Nickolay
    savefield(flatfield_binned,'flatfield_binned_CROPD_'+channel, directory='ToNickolay',mat=True, plot=True)
    savefield(flatfield,'flatfield_unbinned_'+channel, directory='ToNickolay', mat=True, plot=True)

    #Compare flatfields to Nickolay's stars
    FirstRow, LastRow, FirstCol, LastCol = get_flatfield_comparison_area_parameters()
    flatfield_compare_area = flatfield[FirstRow:LastRow, FirstCol:LastCol]
    figstar, axstar = plt.subplots(2, 1)
    plot_CCDimage(flatfield_compare_area,figstar, axstar[0], title='Flatfield compare area '+channel)
    # plot the horizonal mean
    axstar[1].plot(np.mean(flatfield_compare_area, axis=0))
    # fit a line to the mean, and print the slope and the error of the slope
    x = np.arange(flatfield_compare_area.shape[1])
    y = np.mean(flatfield_compare_area, axis=0)
    V=np.polyfit(x, y, 1, cov=True)
    axstar[1].plot(x, V[0][0]*x+V[0][1])
    axstar[1].set_title('Mean of flatfield compare area')
    # print the slope on the plot with 2 decimals
    

    axstar[1].text(0.1, 0.9, 'slope: '+str(V[0][0]), transform=axstar[1].transAxes)

    # print the constant on the plot
    axstar[1].text(0.1, 0.8, 'constant: '+str(V[0][1]), transform=axstar[1].transAxes)
    # print the error estimate of the constant on the plot
    #axstar[1].text(0.1, 0.7, 'constant error: '+str(np.sqrt(V[1][1][1])), transform=axstar[1].transAxes)   

    figstar.savefig('output/flatfield_compare_area_'+channel + ".jpg")



    #printsettings(CCDitems[i])
# %%

# %%
