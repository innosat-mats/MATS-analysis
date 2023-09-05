#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from database_generation.flatfield import read_pic_and_picd
from mats_l1_processing.L1_calibrate import L1_calibrate, flatfield_calibration
from database_generation.experimental_utils import plot_CCDimage, calibrate_CCDitems
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
import datetime as DT
from mats_utils.plotting.plotCCD import orbit_plot, simple_plot, plot_image
import matplotlib.pyplot as plt
from mats_utils.plotting.plotCCD import all_channels_plot
from mats_utils.plotting.animate import generate_gif
import numpy as np
import copy
from mats_utils.geolocation import satellite as satellite

def calibrate(CCDitem, instrument):
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated,
        errors
    ) = L1_calibrate(CCDitem, instrument)
#%%

calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/scripts/calibration_data_linda.toml'
instrument=Instrument(calibration_file)


#%%
# Read in anything just to get a CCDitem structure for IR1
start_time=DT.datetime(2023,2,13,9,7,0)
stop_time=DT.datetime(2023,2,13,9,10,0)

df = read_MATS_data(start_time,stop_time,version='0.4')
CCDitems = dataframe_to_ccd_items(df)
ccdnames=('IR1','IR4','IR3','IR2','UV1','UV2','NADIR')
ir1=dataframe_to_ccd_items(df[df.CCDSEL==1])
ir4=dataframe_to_ccd_items(df[df.CCDSEL==2])
ir3=dataframe_to_ccd_items(df[df.CCDSEL==3])
ir2=dataframe_to_ccd_items(df[df.CCDSEL==4])
uv1=dataframe_to_ccd_items(df[df.CCDSEL==5])
uv2=dataframe_to_ccd_items(df[df.CCDSEL==6])



#%%




CCDlist=[ir1[0],ir2[0],ir3[0],ir4[0],uv1[0],uv2[0]]


directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/Flatfield20210421/PayloadImages/'

for CCD in CCDlist:

    CCDitem=copy.deepcopy(CCD)
    pic, picd=read_pic_and_picd(CCDitem['channel'],directory)
    image_before=pic-picd
    CCDitem['IMAGE']=image_before

    CCDitem["CCDunit"] =instrument.get_CCD(CCDitem["channel"])

    image_calib_nonflipped, error_flags_flatfield = flatfield_calibration(CCDitem)

    fig, ax=plt.subplots(2,1)

    image_calib_nonflipped_scaled=image_calib_nonflipped/image_calib_nonflipped[300:400, 800:1200].mean()
    image_before_scaled=image_before/image_before[300:400,800:1200].mean()
    sp_before=plot_CCDimage(image_before_scaled, fig, ax[0],title=CCDitem["channel"]+' uncalibrated image')
    clim1=sp_before.get_clim()
    plot_CCDimage(image_calib_nonflipped_scaled, fig, ax[1], title='calibrated image', clim=clim1)





"""
(
    image_lsb,
    image_bias_sub,
    image_desmeared,
    image_dark_sub,
    image_calib_nonflipped,
    image_calibrated,
    errors
) = L1_calibrate(CCDitem, instrument)

fig, ax=plt.subplots(7,1)

plot_CCDimage(CCDitem['IMAGE'], fig, ax[0],title='uncalibrated image')
plot_CCDimage(image_lsb, fig, ax[1],title='image_lsb')
plot_CCDimage(image_bias_sub, fig, ax[2],title='image_bias_sub')
plot_CCDimage(image_desmeared, fig, ax[3],title='image_desmeared')
plot_CCDimage(image_dark_sub, fig, ax[4],title='image_dark_sub')
plot_CCDimage(image_calib_nonflipped, fig, ax[5],title='image_calib_nonflipped')
plot_CCDimage(image_calibrated, fig, ax[6],title='image_calibrated')
#plot_CCDimage(CCDitem['image_calibrated'], fig, ax[7], title='calibrated image')
plt.tight_layout
"""

# #%%
# fig, ax=plt.subplots(3,1)
# plot_CCDimage(CCDitem_orig['IMAGE'], fig, ax[0], title='org image')
# plot_CCDimage(CCDitem['IMAGE'], fig, ax[1],title='uncalibrated image')
# plot_CCDimage(CCDitem['image_calibrated'], fig, ax[2], title='calibrated image')




# %%
