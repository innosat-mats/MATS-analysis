
#%%

    # THIS SCRIPT MUST BE RUN FROM MATS-L1-processing in  src/database_generation/read_and_bin_flatfields.py
    # This script is what first what used to give flatfielss to Nickolay. 
    # It used to be run from MATS-L1-processing in  src/database_generation/read_and_bin_flatfields.py

    # makes flatfield using both a cold flatfield without baffle and a room temp flatfield with baffle.


from database_generation.flatfield import make_flatfield
import toml
from mats_l1_processing.instrument import Instrument
from database_generation.experimental_utils import plot_CCDimage
import matplotlib.pyplot as plt
from mats_l1_processing.L1_calibration_functions import bin_image_with_BC
import datetime as DT
from mats_utils.rawdata.read_data import read_MATS_data
import sys
sys.path.append('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from lindas_own_functions import rename_CCDitem_entries

#%%

calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
channels=['IR1','IR2','IR3','IR4','UV1','UV2']#,'NADIR' ]

sigmodes=['HSM','LSM']
starttime=DT.datetime(2023,4,10,0,14,0)
endtime=DT.datetime(2023,4,10,0,15,0)
df=read_MATS_data(starttime, endtime, level='1a',version='0.5')
rename_CCDitem_entries(df)
print(len(df))

#%%
#calibration_data = toml.load(calibration_file)
#instrument = Instrument(calibration_file)

CCDlist= ['IR1'],'IR2','IR3','IR4']#,'UV1','UV2','NADIR']
for i, channel in enumerate(CCDlist):

    flatfield_morphed, flatfield_w_baffle, flatfield_wo_baffle=make_flatfield(channel, sigmodes[0], calibration_file, plot=False, outputrawfields=True)



    fig, ax=plt.subplots(2,1)
    plot_CCDimage(flatfield_w_baffle, fig, ax[0],title='with baffle '+channel)
    plot_CCDimage(flatfield_wo_baffle, fig, ax[1],title='without baffle '+channel)

    # Bin the fields
    # Read in a CCDitem to read teh binning settings from
    CCDitem=df[df.channel==channel].iloc[0]


    bflatfield_w_baffle = bin_image_with_BC(CCDitem, flatfield_w_baffle)
    bflatfield_wo_baffle = bin_image_with_BC(CCDitem, flatfield_wo_baffle)
    fig, ax=plt.subplots(2,1)
    plot_CCDimage(bflatfield_w_baffle, fig, ax[0],title='binned with baffle '+channel)
    plot_CCDimage(bflatfield_wo_baffle, fig, ax[1],title='binned without baffle '+channel)

    # Save fields
    dirname='/Users/lindamegner/MATS/MATS-retrieval/data/irr/'
    bflatfield_wo_baffle.astype('uint32').tofile(dirname+'binned_flatf_wobaffle_'+channel)
    bflatfield_w_baffle.astype('uint32').tofile(dirname+'binned_flatf_wbaffle_'+channel)



# %%
