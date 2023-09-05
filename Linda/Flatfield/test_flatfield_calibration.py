#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from database_generation.flatfield import read_pic_and_picd
from mats_l1_processing.L1_calibrate import L1_calibrate
from database_generation.experimental_utils import plot_CCDimage, calibrate_CCDitems
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
import datetime as DT
from mats_utils.plotting.plotCCD import orbit_plot, simple_plot, plot_image
import matplotlib.pyplot as plt
from mats_utils.plotting.plotCCD import all_channels_plot
from mats_utils.plotting.animate import generate_gif
import copy
from mpldatacursor import datacursor

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


CCDitem_orig=ir1[0]

#%%


CCDitem=copy.deepcopy(CCDitem_orig)






directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/Flatfield20210421/PayloadImages/'
pic, picd=read_pic_and_picd('IR1',directory)
#CCDitem['IMAGE']=pic



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
datacursor(lines)
plt.show(block=True)

# #%%
# fig, ax=plt.subplots(3,1)
# plot_CCDimage(CCDitem_orig['IMAGE'], fig, ax[0], title='org image')
# plot_CCDimage(CCDitem['IMAGE'], fig, ax[1],title='uncalibrated image')
# plot_CCDimage(CCDitem['image_calibrated'], fig, ax[2], title='calibrated image')




# %%
