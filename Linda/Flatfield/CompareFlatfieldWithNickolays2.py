#%% Compares my flatfield to Nickolays

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

#%%

starttime=DT.datetime(2023,4,10,0,14,0)
endtime=DT.datetime(2023,4,10,0,15,0)

# # #%% Select on explicit time NLC
starttime = DT.datetime(2023, 2, 2, 19, 38, 0) 
endtime = DT.datetime(2023, 2, 2, 19, 50, 0)

#ir1=df[df.CCDSEL==1]
#ir2=df[df.CCDSEL==4]
#ir3=df[(df.CCDSEL==3)]
#ir4=df[(df.CCDSEL==2)]
#uv1=df[(df.CCDSEL==5)]
#uv2=df[(df.CCDSEL==6)]
#filter_IR1={'TPsza':[97,150],'CCDSEL': [1,1], 'NRBIN': [2, 2],'NCBINCCDColumns': [40, 40],'NCOL':[43,43], 'NROW':[187,187], 'NRSKIP':[109,109]} 
#ir1l1a=read_MATS_data(starttime, endtime,filter_IR1, level='1a',version='0.5')
df=read_MATS_data(starttime, endtime, level='1a',version='0.9')

print(len(df))
df=df[:20]
#%%

#Calibrate
if not 'instrument' in locals():

    calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml'
    instrument=Instrument(calibration_file)  
#calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
#instrument=Instrument(calibration_file)
rename_CCDitem_entries(df)
#%%

CCDlist= ['IR1', 'IR2','IR3','IR4','UV1','UV2']#,'NADIR']
fig, ax=plt.subplots(len(CCDlist),1)
for i, channel in enumerate(CCDlist):
    CCD =instrument.get_CCD(channel)
    aflatfield=CCD.flatfield() 
    plot_CCDimage(aflatfield, fig, ax[i], clim=[0.8, 1.2], title='fullframe '+channel)



# %%

for i, channel in enumerate(CCDlist):

    fig, ax=plt.subplots(2,1)
    CCD =instrument.get_CCD(channel)
    clim=[0.9, 1.1]
    CCDitem=df[df.channel==channel].iloc[0]
    CCDitem["CCDunit"] =instrument.get_CCD(CCDitem["channel"])
    unbinnedflatfield=CCD.flatfield() 
    
    image_flatf_fact,  error_flag_largeflatf = calculate_flatfield(CCDitem)


    plot_CCDimage(unbinnedflatfield, fig, ax[0], clim=clim, title='fullframe '+channel)
    plot_CCDimage(image_flatf_fact, fig, ax[1], clim=clim, title='binned '+channel)


    path='/Users/lindamegner/MATS/MATS-retrieval/data/flatfields_to_nickolay'
    image_flatf_fact_times_million=image_flatf_fact*1e6
    unbinnedflatfield_times_million=unbinnedflatfield*1e6
    image_flatf_fact_times_million.astype('uint32').tofile(path+'/mergedflatfield_binned_'+channel)# %%
    unbinnedflatfield_times_million.astype('uint32').tofile(path+'/mergedflatfield_unbinned_'+channel)# %%

    path='/Users/lindamegner/MATS/MATS-retrieval/data/flatfields_to_lukas'
# old stuff 

# flipped_channels=['IR1','IR3','UV1','UV2']
# for i, channel in enumerate(CCDlist):


#     fig, ax=plt.subplots(2,1)
#     CCD =instrument.get_CCD(channel)
#     clim=[0.9, 1.1]

#     #get fullframe flatfield
#     flatfield_full=CCD.flatfield('High') #First time I sent to Nickolay
#     #flatfield_morphed=flatfield.make_flatfield(channel, calibration_file, plotresult=True, plotallplots=False)

#     if channel in flipped_channels:
#          flatfield_full=np.fliplr(flatfield_full)
   

#     #get binned flatfield
#     CCDitem=df[df.channel==channel].iloc[0]
#     npixinbin=CCDitem['NCBINFPGAColumns']*CCDitem['NCBINCCDColumns']*CCDitem['NRBIN']
#     CCDitem["CCDunit"] =instrument.get_CCD(CCDitem["channel"])

#     #CCDitem["CCDunit"].flatfield('High') 
#     flatfield_binned =calculate_flatfield(CCDitem)
#     if channel in flipped_channels:
#         flatfield_binned=np.fliplr(flatfield_binned)
#     #flatfield_binned_normalised=flatfield_binned/flatfield_binned.mean()


#     plot_CCDimage(flatfield_full, fig, ax[0], clim=clim, title='fullframe '+channel)
#     plot_CCDimage(flatfield_binned/npixinbin, fig, ax[1], clim=clim, title='binned '+channel)


#     path='/Users/lindamegner/MATS/MATS-retrieval/data/flatfields_to_nickolay'
#     flatfield_binned.astype('uint32').tofile(path+'/mergedflatfield_binned_'+channel)

# %%
