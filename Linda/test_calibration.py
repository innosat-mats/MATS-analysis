#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from database_generation.experimental_utils import plot_CCDimage, calibrate_CCDitems
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
import datetime as DT
from mats_utils.plotting.plotCCD import orbit_plot, simple_plot, plot_image
import matplotlib.pyplot as plt
from mats_utils.plotting.plotCCD import all_channels_plot
from mats_utils.plotting.animate import generate_gif

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
    return image_calibrated
#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)



# start_time=DT.datetime(2023,2,13,9,5,0)
# stop_time=DT.datetime(2023,2,13,9,15,0)

start_time=DT.datetime(2023,2,12,1,0,0)
stop_time=DT.datetime(2023,2,12,1,14,0)

df=read_MATS_data(start_time, stop_time)



#%%
print(df.columns.tolist())
#CCDitems = dataframe_to_ccd_items(df)
CCDitems=df

#%%





CCDitems['NCBIN CCDColumns'] = CCDitems['NCBINCCDColumns'] 
CCDitems['GAIN Truncation'] = CCDitems['GAINTruncation'] 
CCDitems['NCBIN FPGAColumns'] = CCDitems['NCBINFPGAColumns']

mylist=[]
CCDitems['BC'] = [mylist for i in CCDitems.index]
CCDitems['GAIN Mode'] = CCDitems['GAINMode'] 


img_cal=calibrate(df.iloc[0], instrument)
#%%

CCDitems['image_calibrated']=CCDitems.apply(lambda CCDitem: calibrate(CCDitem,instrument ),axis=1)

#%%
ir1=CCDitems[CCDitems.channel=='IR1']

uv2=CCDitems[CCDitems.channel=='UV2']
#calibrate_CCDitems(CCDitems, Instrument(calibration_file))

fig, ax=plt.subplots(1)
plot_CCDimage(ir1.iloc[0]['image_calibrated'], fig, ax, clim=[0, 200])


# def plot_image(CCD, ax=None, fig=None, outpath=None,
#                nstd=2, cmap='inferno', ranges=None,
#                optimal_range=False, format='png', save=True,
#                fontsize=10):
# %%
#%%

all_channels_plot(uv2[:2], './output/test/', optimal_range=True)
generate_gif('./output/test/ALL','output/film_staticdots.gif' )