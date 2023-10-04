
#%%    
from mats_utils.rawdata.read_data import read_MATS_data
from database_generation.experimental_utils import plot_CCDimage
from database_generation.flatfield import read_flatfield
from mats_l1_processing.instrument import Instrument
import datetime as DT
from mats_utils.rawdata.calibration import  calibrate_dataframe
from mats_l1_processing.L1_calibration_functions import calculate_flatfield, bin_image_with_BC, meanbin_image_with_BC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from lindas_own_functions import rename_CCDitem_entries
from database_generation import flatfield as flatfield
from mats_l1_processing.instrument import CCD
import toml

#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'

CCDlist= ['IR1', 'IR2','IR3','IR4','UV1','UV2']#
for i, channel in enumerate(CCDlist):
    CCDunit=CCD(channel,calibration_file)
    #instrument=Instrument(calibration_file)
    #CCD =instrument.get_CCD(channel)

    calibration_data=toml.load(calibration_file)
    flatfield_wo_baffle_HSM= read_flatfield(CCDunit, 'HSM', calibration_data["flatfield"]["flatfieldfolder_cold_unprocessed"])
    flatfield_wo_baffle_LSM= read_flatfield(CCDunit, 'LSM', calibration_data["flatfield"]["flatfieldfolder_cold_unprocessed"])
    
    flatfield_wo_baffle_HSMm=flatfield_wo_baffle_HSM/flatfield_wo_baffle_HSM.mean()
    flatfield_wo_baffle_LSMm=flatfield_wo_baffle_LSM/flatfield_wo_baffle_LSM.mean()


    fig, ax=plt.subplots(2,1)
    HSM=plot_CCDimage(flatfield_wo_baffle_HSMm, fig, ax[0],  title=channel+' HSM')
    climHSM=HSM.get_clim()
    LSM=plot_CCDimage(flatfield_wo_baffle_LSMm, fig, ax[1], clim=climHSM, title=channel+' LSM')



# %%
