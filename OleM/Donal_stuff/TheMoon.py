#%%
from pyarrow import fs, schema, string
from pyarrow.dataset import FilenamePartitioning
import pyarrow.dataset as ds
import boto3
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from mats_utils.rawdata.read_data import read_MATS_data
import matplotlib.pylab as plt
import matplotlib.colors as colors
from scipy import optimize
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation as R
import io
from PIL import Image, ImageFilter
from mats_l1_processing.read_parquet_functions import *
from mats_l1_processing.L1_calibration_functions import *
from mats_l1_processing import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_utils.geolocation.coordinates import *
from mats_utils.geolocation import satellite as satellite
from skyfield import api as sfapi
import skyfield.sgp4lib as sgp4lib
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from pandas import concat, DataFrame

#%%
starttime = DT.datetime(2023,2, 13, 0, 30, 0)
endtime = DT.datetime(2023, 2, 13, 12, 45, 0)
ccd_data = read_MATS_data(starttime, endtime,filter=None,level='1a',version='0.5')


#%%
calibration_file ="calibration_data.toml"    
instrument = Instrument(calibration_file)

#%%
ccd_items = dataframe_to_ccd_items(
        ccd_data,
        remove_empty=False,
        remove_errors=False,
        remove_warnings=False,
    )

#%%
for ccd in ccd_items:
    if ccd["IMAGE"] is None:
        image_calibrated = None
        errors = None
    else:
        (
            _,
            _,
            _,
            _,
            _,
            _,
            image_calibrated,
            errors,
        ) = L1_calibrate.L1_calibrate(ccd, instrument, force_table=False)
    ccd["ImageCalibrated"] = image_calibrated
    ccd["CalibrationErrors"] = errors

calibrated = DataFrame.from_records(
    ccd_items,
    columns=[
        "ImageCalibrated",
        "CalibrationErrors",
        "qprime",
        "channel",
        "flipped",
        "temperature",
        "temperature_HTR",
        "temperature_ADC",
    ],
)
l1b_data = concat([
    ccd_data,
    calibrated,
], axis=1).set_index("EXPDate").sort_index()
l1b_data.drop(["ImageData", "Errors", "Warnings"], axis=1, inplace=True)
l1b_data = l1b_data[l1b_data.ImageCalibrated != None]  # noqa: E711
l1b_data["ImageCalibrated"] = [
    ic.tolist() for ic in l1b_data["ImageCalibrated"]
]
l1b_data["CalibrationErrors"] = [
    ce.tolist() for ce in l1b_data["CalibrationErrors"]
]