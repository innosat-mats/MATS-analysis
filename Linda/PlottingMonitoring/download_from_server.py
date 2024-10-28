#%%
from mats_utils.rawdata.read_data import load_multi_parquet
import datetime as DT
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.plotting.plotCCD import plot_image

#%%
start = DT.datetime(2023, 2, 11, 12, 0, 0)
stop= DT.datetime(2023, 2, 11, 15, 0, 0)
#df=load_multi_parquet('/home/linda/MATS_data/CROPD_v0.6',start, stop)
filter_UV2={'CCDSEL': [6,6]} 

df=load_multi_parquet('../../../data/CROPD_v0.6',start, stop, filt=filter_UV2)


datafolder='/home/linda/data/tmp'
for index, CCD in df[:5].iterrows():
    #plot_CCDimage(CCD['ImageCalibrated'])
    plot_image(CCD, outpath=datafolder)
#''

# %%
