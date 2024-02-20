# Author: Linda Megner
# Created: October 26, 2023


#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
import datetime as DT
from mats_utils.rawdata.calibration import  calibrate_dataframe
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from database_generation import flatfield as flatfield
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items, CCDitems_to_dataframe
import pickle


from mats_l1_processing.L1_calibration_functions import (
    get_true_image,
    desmear_true_image,
    subtract_dark,
    flatfield_calibration,
    get_linearized_image,
    combine_flags,
    make_binary,
    flip_image,
    handle_bad_columns,
    artifact_correction,
    correct_hotpixels,
    correct_single_events
)
 

def wrap_get_CCD(CCDitem, instrument):
    CCDunit=instrument.get_CCD(CCDitem["channel"])
    return CCDunit  
def wrap_correct_single_events(CCDitem):
    image_se_corrected, error_flags_se  = correct_single_events(CCDitem,CCDitem['IMAGE'])
    return image_se_corrected
def wrap_correct_hotpixels(CCDitem):
    image_hot_pixel_corrected, error_flags_hp  = correct_hotpixels(CCDitem,CCDitem["image_se_corrected"])
    return image_hot_pixel_corrected
def wrap_get_true_image(CCDitem):
    image_bias_sub, error=get_true_image(CCDitem, CCDitem["image_hot_pixel_corrected"])
    return image_bias_sub
def wrap_get_linearized_image(CCDitem):    
    image_linear, error = get_linearized_image(CCDitem, CCDitem["image_bias_sub"], force_table=True)
    return image_linear
def wrap_desmear_true_image(CCDitem):
    image_desmeared, error=desmear_true_image(CCDitem, CCDitem["image_linear"])
    return image_desmeared
def wrap_subtract_dark(CCDitem):
    image_dark_sub, error=subtract_dark(CCDitem, CCDitem["image_desmeared"])
    return image_dark_sub
def wrap_flatfield_calibration(CCDitem):
    image_flatfielded, error=flatfield_calibration(CCDitem, CCDitem["image_dark_sub"])
    return image_flatfielded
def wrap_flip_image(CCDitem):
    image_flipped=flip_image(CCDitem, CCDitem["image_flatfielded"])
    return image_flipped



# Calibrate the images so that all steps are saved
def calibration_in_steps(CCDitem, instrument):
    CCDitem["CCDunit"] = wrap_get_CCD(CCDitem, instrument)
    CCDitem["image_lsb"] = CCDitem['IMAGE']
    CCDitem["image_se_corrected"] = wrap_correct_single_events(CCDitem)
    CCDitem["image_hot_pixel_corrected"] = wrap_correct_hotpixels(CCDitem)
    CCDitem["image_bias_sub"] = wrap_get_true_image(CCDitem)
    CCDitem["image_linear"] = wrap_get_linearized_image(CCDitem)
    CCDitem["image_desmeared"] = wrap_desmear_true_image(CCDitem)
    CCDitem["image_dark_sub"] = wrap_subtract_dark(CCDitem)
    CCDitem["image_flatfielded"] = wrap_flatfield_calibration(CCDitem)
    CCDitem["image_flipped"] = wrap_flip_image(CCDitem)
            #ccd["ImageCalibrated"] = image_calibrated
            #ccd["CalibrationErrors"] = errors
    return CCDitem


def calibration_of_df_in_steps(df, instrument):
    CCDitems = dataframe_to_ccd_items(df)
    for CCDitem in CCDitems:
        CCDitem = calibration_in_steps(CCDitem, instrument)
    df_with_calib_steps=CCDitems_to_dataframe(CCDitems)
    return df_with_calib_steps





#%%




# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'
starttime=DT.datetime(2023,5,7,1,0,0)
endtime=DT.datetime(2023,5,7,2,0,0)


filter_UV2={'CCDSEL': [6,6]} 
df=read_MATS_data(starttime, endtime,filter_UV2, level='1a',version='0.6')
print(len(df))
print(df.columns)

pickle.dump(df, open('testdata/df_in_orbit_NLCuv2.pkl', 'wb'))

#%%

with open('testdata/df_in_orbit_NLCuv2.pkl', 'rb') as f:
    df = pickle.load(f)


calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)

#%%

#df= calibration_of_df_in_steps(df[:5], instrument)
df=calibrate_dataframe(df[:5], instrument,debug_outputs=True)

for index, CCD in df.iterrows():
    plot_image(CCD, save=False, image_field='image_flatfielded')
    

# %%
