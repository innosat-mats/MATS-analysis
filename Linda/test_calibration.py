#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage, calibrate_CCDitems
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
import datetime as DT
from mats_utils.plotting.plotCCD import orbit_plot, simple_plot, plot_image
import matplotlib.pyplot as plt

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
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-L1-processing/scripts/calibration_data_linda.toml'




start_time=DT.datetime(2023,1,30,18,0,0)
stop_time=DT.datetime(2023,1,30,18,10,0)

df=read_MATS_data(start_time, stop_time)
print(df.columns.tolist())
#CCDitems = dataframe_to_ccd_items(df)
CCDitems=df

#%%

instrument=Instrument(calibration_file)



CCDitems['NCBIN CCDColumns'] = CCDitems['NCBINCCDColumns'] 
CCDitems['GAIN Truncation'] = CCDitems['GAINTruncation'] 
CCDitems['NCBIN FPGAColumns'] = CCDitems['NCBINFPGAColumns']
mylist=[]
CCDitems['BC'] = [mylist for i in CCDitems.index]
CCDitems['GAIN Mode'] = CCDitems['GAINMode'] 

CCDitems['image_calibrated']=CCDitems.apply(lambda CCDitem: calibrate(CCDitem,instrument ),axis=1)

#%%
ir1=CCDitems[CCDitems.channel=='IR1']

#calibrate_CCDitems(CCDitems, Instrument(calibration_file))

fig, ax=plt.subplots(1)
plot_CCDimage(ir1.iloc[0]['image_calibrated'], fig, ax)


# def plot_image(CCD, ax=None, fig=None, outpath=None,
#                nstd=2, cmap='inferno', ranges=None,
#                optimal_range=False, format='png', save=True,
#                fontsize=10):
# %%
