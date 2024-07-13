#%%
from mats_utils.rawdata.read_data import read_MATS_data, read_MATS_PM_data
import datetime as DT
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_l1_processing.read_in_functions import read_CCDitems 
import numpy as np
import pickle
import matplotlib.pyplot as plt
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.rawdata.calibration import calibrate_dataframe



# # #%% Select on explicit time NLC

start_time = DT.datetime(2023, 2, 2, 19, 38, 0)
stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)
#start_time = DT.datetime(2023, 1, 4, 19, 2, 10)
#stop_time = DT.datetime(2023, 1, 4, 19, 2, 40)

df = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)

n=2
uv1=df[df.channel=='UV1'][:n]
uv2=df[df.channel=='UV2'][:n]
ir1=df[df.channel=='IR1'][:n]
ir2=df[df.channel=='IR2'][:n]
ir3=df[df.channel=='IR3'][:n]
ir4=df[df.channel=='IR4'][:n]


# pickle.dump(CCDitemsuv1, open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'wb'))
# pickle.dump(CCDitemsuv2, open('testdata/CCD_items_in_orbit_NLCuv2.pkl', 'wb'))

# #%%
# #with open('testdata/CCD_items_in_orbit_UVIR.pkl', 'rb') as f:
# with open('testdata/CCD_items_in_orbit_NLCuv1.pkl', 'rb') as f:
#     CCDitems = pickle.load(f)
#%%
instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')

uv1cal=calibrate_dataframe(uv1, instrument, debug_outputs=True)
uv2cal=calibrate_dataframe(uv2, instrument, debug_outputs=True)



# %%


df_l1b = read_MATS_data(start_time,stop_time,version='0.6',level='1b',dev=False)
#%%
uv1_l1b=df_l1b[df_l1b.channel=='UV1'][:n]
uv2_l1b=df_l1b[df_l1b.channel=='UV2'][:n]



# %%
for i in range(n):
    fig, ax = plt.subplots(2,2)

    plot_CCDimage(uv1_l1b.iloc[i]["ImageCalibrated"], fig=fig, axis=ax[0,0], title=uv1_l1b.iloc[i]["channel"]+' AWScal '+str(uv1_l1b.iloc[i]["TMHeaderTime"])[0:19])
    ax[0,0].text(10, 60, 'max: '+str(np.max(uv1_l1b.iloc[i]["ImageCalibrated"])), color='white')
    ax[0,0].text(10, 80, 'min: '+str(np.min(uv1_l1b.iloc[i]["ImageCalibrated"])), color='white')
    plot_CCDimage(uv2_l1b.iloc[i]["ImageCalibrated"], fig=fig, axis=ax[0,1], title=uv2_l1b.iloc[i]["channel"]+' AWScal '+str(uv2_l1b.iloc[i]["TMHeaderTime"])[0:19])
    ax[0,1].text(10, 60, 'max: '+str(np.max(uv2_l1b.iloc[i]["ImageCalibrated"])), color='white')
    ax[0,1].text(10, 80, 'min: '+str(np.min(uv2_l1b.iloc[i]["ImageCalibrated"])), color='white')
    plot_CCDimage(uv1cal.iloc[i]["ImageCalibrated"], fig=fig, axis=ax[1,0],title='offline calib')
    ax[1,0].text(10, 60, 'max: '+str(np.max(uv1cal.iloc[i]["ImageCalibrated"])), color='white')
    ax[1,0].text(10, 80, 'min: '+str(np.min(uv1cal.iloc[i]["ImageCalibrated"])), color='white')
    plot_CCDimage(uv2cal.iloc[i]["ImageCalibrated"], fig=fig, axis=ax[1,1],title='offline calib')
    ax[1,1].text(10, 60, 'max: '+str(np.max(uv2cal.iloc[i]["ImageCalibrated"])), color='white')
    ax[1,1].text(10, 80, 'min: '+str(np.min(uv2cal.iloc[i]["ImageCalibrated"])), color='white')
    plt.tight_layout()
    


# %%
