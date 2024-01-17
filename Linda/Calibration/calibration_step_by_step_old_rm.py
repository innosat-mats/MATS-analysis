#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from database_generation.experimental_utils import plot_CCDimage, calibrate_CCDitems
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from lindas_own_functions import rename_CCDitem_entries

def calibrate(CCDitem, instrument):
    (
    #     image_lsb, 
    #     image_se_corrected, 
    #     image_hot_pixel_corrected, 
    #     image_bias_sub, 
    #     image_desmeared, 
    #     image_dark_sub, 
    #     image_calib_nonflipped, 
    #     image_calib_flipped, 
    #     image_calibrated, 
    #     errors


    image_lsb, 
    image_bias_sub, 
    image_desmeared, 
    image_dark_sub, 
    image_calib_nonflipped, 
    image_calib_flipped, 
    image_calibrated, 
    errors

    ) = L1_calibrate(CCDitem, instrument)
    return   image_calibrated
#%%
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)

start_time=DT.datetime(2023,5,7,1,0,0)
stop_time=DT.datetime(2023,5,7,2,0,0)


filter_UV2={'CCDSEL': [6,6]} 
#df=read_MATS_data(start_t
# ime, stop_time, filter_UV2, level='1a')
df=read_MATS_data(start_time, stop_time, filter_UV2, level='1a',version='0.6')
#df=read_MATS_data(start_time, stop_time, level='1b')
print(len(df))


print(df.columns.tolist())
#CCDitems = dataframe_to_ccd_items(df)
CCDitems=df[:min(10,len(df))]

CCDitems=rename_CCDitem_entries(CCDitems)


#%%

img_cal=calibrate(CCDitems.iloc[0], instrument)
#%%

CCDitems['image_calibrated']=CCDitems.apply(lambda CCDitem: calibrate(CCDitem,instrument ),axis=1)




#%%

fig, ax=plt.subplots(2,1)
plot_CCDimage(CCDitems.iloc[0]['IMAGE'], fig, ax[0])
plot_CCDimage(CCDitems.iloc[0]['image_calibrated'], fig, ax[1])
#plot_CCDimage(ir1.iloc[0]['image_calibrated'], fig, ax[0], clim=[0,7])
#plot_CCDimage(ir2.iloc[0]['image_calibrated'], fig, ax[1],clim=[0,7])

# def plot_image(CCD, ax=None, fig=None, outpath=None,
#                nstd=2, cmap='inferno', ranges=None,
#                optimal_range=False, format='png', save=True,
#                fontsize=10):
# %%



#%%
"""
all_channels_plot(uv2[:2], './output/test/', optimal_range=True)
generate_gif('./output/test/ALL','output/film_staticdots.gif' )

"""
# %%
