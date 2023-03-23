#%% Donal script to calibrate
import datetime as DT
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
import matplotlib.pyplot as plt
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage

calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/scripts/calibration_data_linda.toml'

start_time=DT.datetime(2023,2,13,9,5,0)
stop_time=DT.datetime(2023,2,13,9,15,0)

df = read_MATS_data(start_time,stop_time,version='0.4')
CCDitems = dataframe_to_ccd_items(df)
#%%
ccdnames=('IR1','IR4','IR3','IR2','UV1','UV2','NADIR')
ir1=dataframe_to_ccd_items(df[df.CCDSEL==1])
ir4=dataframe_to_ccd_items(df[df.CCDSEL==2])
ir3=dataframe_to_ccd_items(df[df.CCDSEL==3])
ir2=dataframe_to_ccd_items(df[df.CCDSEL==4])
uv1=dataframe_to_ccd_items(df[df.CCDSEL==5])
uv2=dataframe_to_ccd_items(df[df.CCDSEL==6])

instrument=Instrument(calibration_file)
def calibrate(CCDitem, instrument):
    (image_lsb, image_bias_sub, 
    image_desmeared, image_dark_sub, 
    image_calib_nonflipped, image_calibrated, errors) = L1_calibrate(CCDitem,instrument)
    return image_calibrated
n=0
plt.close('all')
ir1cal=calibrate(ir1[n],instrument)
ir2cal=calibrate(ir2[n],instrument)
ir3cal=calibrate(ir3[n],instrument)
ir4cal=calibrate(ir4[n],instrument)
uv1cal=calibrate(uv1[n],instrument)
uv2cal=calibrate(uv2[n],instrument)
# %%
fig, ax=plt.subplots(1)

plot_CCDimage(ir1cal, fig, ax)
# %%
