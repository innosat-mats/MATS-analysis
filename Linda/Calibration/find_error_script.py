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
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
#from lindas_own_functions import plot_CCDimage_transp

def add_error_est_to_dataframe(df, calibration_file):
    """
    docstring
    df: dataframe for which to add error
    """
    CCDitems = dataframe_to_ccd_items(df, calibration_file)
    
    add_error_est_to_CCDitems(CCDitems, calibration_file)
    

def add_error_est_to_CCDitems(CCDitems, calibration_file):
    """
    Takes a CCDitem and adds an error estimate to it based on the calibration errors
    for the different steps. The error estimate is added as a new key in the CCDitem.

    CCDitem: CCDitem for which to add error
    """
    import toml
    calibration_data = toml.load(calibration_file)
    instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

    
    
    for i in range(0,len(CCDitems)):
        CCDitems[i]['CCDunit'] =instrument.get_CCD(CCDitems[i]['channel'])
        calib_denominator=CCDitems[i]['CCDunit'].calib_denominator(CCDitems[i]['GAIN Mode'])
        #this can also be done using L1_calibration_functions.absolute_calibration
        CCDitems[i]['FlatfieldError']=get_flatfield_error(CCDitems[i], calibration_data)/calib_denominator

        CCDitems[i]['DarkCurrentError']=get_darkcurrent_error(CCDitems[i])
        CCDitems[i]['TotSysError'] = np.sqrt(CCDitems[i]['FlatfieldError']**2+CCDitems[i]['DarkCurrentError']**2)
        

def bin_abs_error(CCDitem, error_nonbinned):
    """
    This is a function to bin an error estimate. 
    Bins according to binning and NSKIP settings in CCDitem.

    Args:
        CCDitem:  dictonary containing CCD image and information
        erroe_nonbinned (optional): numpy array with error estimate for each pixel

    Returns:
        binned_error: binned error (currently assumes that the error is the same for all subpixels in a superpixel)

    """
    totbin = int(CCDitem["NRBIN"])*int(CCDitem["NCBIN CCDColumns"]) * \
        int(CCDitem["NCBIN FPGAColumns"])

    #assumes that the error is the same for all subpixels in a superpixel, 
    # could be imoproved by taking the root of the squared sum of the error of all subpixels
    # NOTE: This is the error of the superpixel, ie of the sum, not the mean of all the pixels
    # htis is why it is multiplied by sqrt(totbin) rather than devided by it
    binned_error=meanbin_image_with_BC(CCDitem, error_nonbinned)*np.sqrt(totbin)


    return binned_error

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

def get_darkcurrent_error(CCDitem):
    """
    Takes a CCDitem and returns the error estimate for the darkvcurrent.

    CCDitem: CCDitem for which to add error
    """
    CCDunit=CCDitem['CCDunit']
    T=CCDitem["temperature"]

    if CCDitem["GAIN Mode"] == 'High':
        log_a_img_avr=CCDunit.log_a_img_avr_HSM
        log_b_img_avr=CCDunit.log_b_img_avr_HSM
        log_a_img_std=CCDunit.log_a_img_err_HSM
        log_b_img_std=CCDunit.log_b_img_err_HSM
    elif CCDitem["GAIN Mode"] == 'Low':
        log_a_img_avr=CCDunit.log_a_img_avr_LSM
        log_b_img_avr=CCDunit.log_b_img_avr_LSM 
        log_a_img_std=CCDunit.log_a_img_err_LSM
        log_b_img_std=CCDunit.log_b_img_err_LSM           
    else:
        raise Exception("Undefined mode")
    

    rawdark=CCDunit.getrawdark(log_a_img_avr, log_b_img_avr, T)

    #Add errors from log_a and log_b since they are correlated - this may be an overestimate - in fact they may be anticorrelated /LM 240513 
    errorab=CCDunit.getrawdark(log_a_img_avr+log_a_img_std, log_b_img_avr+log_b_img_std, T)-rawdark
    deltaT=3
    errorT=CCDunit.getrawdark(log_a_img_avr, log_b_img_avr, T+deltaT)-rawdark
    error=np.sqrt(errorT**2+errorab**2)


    totdarkcurrent2Derr=error* int(CCDitem["TEXPMS"])/ 1000.0
    




    dark_calc_err_image = (
        CCDunit.ampcorrection
        * totdarkcurrent2Derr
        / CCDunit.alpha_avr(CCDitem["GAIN Mode"])
    )

    darkcurrent_err_binned = bin_abs_error(CCDitem, dark_calc_err_image)
    return darkcurrent_err_binned 



def get_flatfield_error(CCDitem, calibration_data):
    """
    Takes a CCDitem and returns the error estimate for the flatfielding step.

    CCDitem: CCDitem for which to add error
    """
    channel = CCDitem['channel']
    #totbin = int(CCDitem["NRBIN"])*int(CCDitem["NCBIN CCDColumns"]) * \
    #    int(CCDitem["NCBIN FPGAColumns"])
    
    

    # The error measured as the standard deviation of three images divided by square root of 3
    flatfield_wo_baffle_err = np.load(
                calibration_data["flatfield"]["flatfieldfolder"]
                + "flatfield_wo_baffle_err_"
                + channel
                + "_HSM.npy")
    # The baffle scalefield, where 1 means no effect, ie the flatfield is the 
    # same as the one without baffle. 
    flatfield_scalefield = np.load(
                calibration_data["flatfield"]["flatfieldfolder"]
                + "baffle_scalefield_"
                + channel
                + "_HSM.npy")

    # Another way to estimate the error is to take the error of the flatfield
    #flatfield_err_rel=flatfield_wo_baffle_err/flatfield_scalefield
    #flatfield_err_rel_binned=meanbin_image_with_BC(CCDitem, flatfield_err_rel) 
    #the above is not correct since it does not reduce the error when binning
    flatfield_err_binned = bin_abs_error(CCDitem, flatfield_wo_baffle_err)
    flatfield_binned=bin_image_with_BC(CCDitem, flatfield_scalefield)


    mark_current_cropping(CCDitem, flatfield_scalefield, flatfield_binned)

    if "ImageCalibrated" in CCDitem.keys(): 
        flatfield_err_nonflipped =CCDitem['ImageCalibrated']/ flatfield_binned * flatfield_err_binned 
    else:
        flatfield_err_nonflipped =absolute_calibration(CCDitem,image=CCDitem['IMAGE'] )/ flatfield_binned * flatfield_err_binned 
        Warning("Image is not calibrated, which will result in a wrong error estimate")

    #relative error=flatfield_wo_baffle_err/flatfield_scalefield
    return flatfield_err_nonflipped    




# #%%
# # # # #%% Select on explicit time
# start_time = DT.datetime(2023, 5, 11, 6, 10)
# stop_time = DT.datetime(2023, 5, 11, 6, 15)
# df = read_MATS_data(start_time,stop_time,version='0.7',level='1b',dev=False)

# CCDitems = dataframe_to_ccd_items(df)
# # pickle.dump(CCDitems, open('testdata/CCD_items_in_orbit_l1b.pkl', 'wb'))

#%%
#
#with open('testdata/CCD_items_in_orbit_UVIR.pkl', 'rb') as f:
with open('testdata/CCD_items_in_orbit_l1b.pkl', 'rb') as f:
    CCDitems = pickle.load(f)
#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
nccditems=7
add_error_est_to_CCDitems(CCDitems[:nccditems], calibration_file)
#%%

for i in range(0,nccditems):
    fig, ax = plt.subplots(4, 1)
    plot_CCDimage(CCDitems[i]['ImageCalibrated'],fig, ax[0], title=CCDitems[i]['channel']+' Original image')
    plot_CCDimage(CCDitems[i]['DarkCurrentError'],fig, ax[1], title='Dark current error estimate')
    plot_CCDimage(CCDitems[i]['FlatfieldError'],fig, ax[2], title='Flatfield error estimate')
    plot_CCDimage(CCDitems[i]['TotSysError'],fig, ax[3], title='Total systematic error estimate')

#%%


#%%
