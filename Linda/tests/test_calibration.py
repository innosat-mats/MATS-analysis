#%%
import pytest

from mats_l1_processing.read_and_calibrate_all_files_parallel import main
from mats_l1_processing.instrument import Instrument, CCD, Photometer
from mats_l1_processing import photometer
import pandas as pd
from mats_l1_processing.L1_calibration_functions import inverse_model_real,inverse_model_table,make_binary,combine_flags,desmear,artifact_correction, correct_single_events,correct_hotpixels, padlastrowsofimage
from mats_l1_processing.L1_calibrate import L1_calibrate

import pickle
import numpy as np
import matplotlib.pyplot as plt
import time

__author__ = "Ole Martin Christensen"
__copyright__ = "Ole Martin Christensen"
__license__ = "MIT"

def test_image_padding():
    with open('/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/testdata/flatfield_IR2_HSM.pkl', 'rb') as f:
        image = pickle.load(f)
    #image = CCDitems[0]['IMAGE']
    [rows, colums]=image.shape
    print('shape of image',image.shape)
    nrow=3
    image_padded=padlastrowsofimage(image,nrow)
    assert np.shape(image_padded)==(rows+nrow,colums), "Image padding shape mismatch"
    assert np.all(image_padded[512:,:]==image_padded[511,:]), "Image padding values mismatch"
    print('shape of padded image',image_padded.shape)
if __name__ == "__main__":
    

    #test_calibrate()
    #test_calibrate_plots()
    #test_error_algebra()
    #test_channel_quaterion()
    #test_photometer()
    #test_hp_correction()
    #test_se_correction()
    test_image_padding()

# %%
