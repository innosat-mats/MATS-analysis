#%%
from mats_utils.rawdata.read_data import read_MATS_data,read_ccd_data_in_interval
import datetime as DT
from mats_l1_processing.read_parquet_functions import dataframe_to_ccd_items
from mats_l1_processing.instrument import Instrument
import mats_utils.error_estimate.error_estimate as ee
import numpy as np
from matplotlib import pyplot as plt

#%% 

def get_flatfield_syserr(CCDitem, calibration_file):
    import toml
    #takes in CCDitem and calulates shot noise
    calibration_data = toml.load(calibration_file)
    channel = CCDitem["channel"]

   
    flatfield = np.load(
                calibration_data["flatfield"]["flatfieldfolder"]
                + "flatfield_"
                + channel
                + "_HSM.npy"
    )
    flatfield_syserr = np.load(
                        calibration_data["flatfield"]["flatfieldfolder"]
                + "flatfield_err_"
                + channel
                + "_HSM.npy"
    )

    return flatfield_syserr

calib_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml'
instrument = Instrument(calib_file)
    
#%% Select on explicit time
start_time = DT.datetime(2023, 3, 19, 20, 0)
stop_time = DT.datetime(2023, 3, 19, 21, 30)

#%%
#df_l1a_from_star= read_MATS_data(start_time,stop_time,level='1a',version='0.6',filter={'schedule_name': 'STAR'})
df_l1a = read_MATS_data(start_time,stop_time,level='1a',version='0.6')
df_l1a = df_l1a[df_l1a.reset_index().index % 10 == 0]
#%%
CCDitems = dataframe_to_ccd_items(df_l1a)


signal = []
shot_noise = []
readout_noise = []
digi_noise = []
jpeg_noise = []
flatfield_syserr = []


for i in range(len(CCDitems)):
    CCDitems[i]["CCDunit"] = instrument.get_CCD(CCDitems[i]["channel"])
    signal.append(CCDitems[i]["IMAGE"])
    shot_noise.append(ee.get_shot_noise(CCDitems[i]))
    readout_noise.append(ee.get_readout_noise(CCDitems[i])*np.ones(CCDitems[i]["IMAGE"].shape))
    digi_noise.append(ee.get_digitization_noise(CCDitems[i])*np.ones(CCDitems[i]["IMAGE"].shape))
    jpeg_noise.append(ee.get_compression_noise(CCDitems[i])*np.ones(CCDitems[i]["IMAGE"].shape))
    flatfield_syserr.append(get_flatfield_syserr(CCDitems[i],calib_file)*np.ones(CCDitems[i]["IMAGE"].shape))   

#%%
tot_noise = []
rel_noise = []
rel_shot = []
rel_ro = []
rel_digi = []
rel_jpeg = []
rel_flatfield = []
for i in range(len(CCDitems)):
    tot_noise.append(np.sqrt(shot_noise[i]**2 + readout_noise[i]**2 + digi_noise[i]**2 + jpeg_noise[i]**2))
    rel_noise.append(tot_noise[i]/signal[i])
    rel_shot.append(shot_noise[i]/signal[i])
    rel_ro.append(readout_noise[i]/signal[i])
    rel_digi.append(digi_noise[i]/signal[i])
    rel_jpeg.append(jpeg_noise[i]/signal[i])
    rel_flatfield.append(flatfield_syserr[i]/signal[i])

# %%
i = 0
fig = plt.figure()

plt.subplot(3, 3, 1)
plt.pcolor(signal[i])
plt.title('signal')
plt.colorbar()

plt.subplot(3, 3, 2)
plt.pcolor(shot_noise[i])
plt.title('shot noise')
plt.colorbar()

plt.subplot(3, 3, 3)
plt.pcolor(readout_noise[i])
plt.title('ro noise')
plt.colorbar()

plt.subplot(3, 3, 4)
plt.pcolor(digi_noise[i])
plt.title('digi noise')
plt.colorbar()

plt.subplot(3, 3, 5)
plt.pcolor(jpeg_noise[i])
plt.title('jpeg noise')
plt.colorbar()

plt.subplot(3, 3, 7)
plt.pcolor(tot_noise[i])
plt.title('tot noise')
plt.colorbar()

plt.subplot(3, 3, 8)
plt.pcolor(rel_noise[i]*100)
plt.title('rel noise')
plt.colorbar()

plt.subplot(3, 3, 9)
plt.pcolor(rel_flatfield[i]*100)
plt.title('rel flatfield')
plt.colorbar()


plt.show()
# %%
# %% take mean at each height
nrows=CCDitems[0]["IMAGE"].shape[0]
signal_z = np.zeros([nrows,len(CCDitems)])*np.nan
shot_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
readout_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
digi_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
jpeg_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
tot_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_noise_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_shot_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_ro_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_digi_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_jpeg_z = np.zeros([nrows,len(CCDitems)])*np.nan
rel_flatfield_z = np.zeros([nrows,len(CCDitems)])*np.nan


for i in range(len(CCDitems)):
    if signal[i].shape[0] == nrows:
        signal_z[:,i] = np.mean(signal[i],axis=1)
        shot_noise_z[:,i] = np.mean(shot_noise[i],axis=1)
        readout_noise_z[:,i] = np.mean(readout_noise[i],axis=1)
        digi_noise_z[:,i] = np.mean(digi_noise[i],axis=1)
        jpeg_noise_z[:,i] = np.mean(jpeg_noise[i],axis=1)

        tot_noise_z[:,i] = np.mean(tot_noise[i],axis=1)
        rel_noise_z[:,i] = np.mean(rel_noise[i],axis=1)

# %%
plt.pcolor(signal_z)
plt.title('Signal (LSB)')
plt.xlabel('time index')
plt.ylabel('vertical pixel')
plt.colorbar()

# %%
plt.pcolor(tot_noise_z)
plt.title('Noise (LSB)')
plt.xlabel('time index')
plt.ylabel('vertical pixel')
plt.colorbar()

# %%
plt.pcolor(rel_noise_z*100)
plt.title('Relative Noise (%)')
plt.xlabel('time index')
plt.ylabel('vertical pixel')
plt.colorbar()
# %%
day = 300

z = np.arange(0,nrows)
plt.plot(np.nanmean(shot_noise_z[:,:day],axis=1),z)
plt.plot(np.nanmean(readout_noise_z[:,:day],axis=1),z)
plt.plot(np.nanmean(digi_noise_z[:,:day],axis=1),z)
plt.plot(np.nanmean(jpeg_noise_z[:,:day],axis=1),z)
plt.plot(np.nanmean(tot_noise_z[:,:day],axis=1),z,'--')
plt.legend(['shot','ro','digi','jpeg','tot'])
plt.title('Mean noise (night)')
plt.xlabel('LSB')
plt.ylabel('vertical pixel')

# %%
z = np.arange(0,nrows)
plt.plot(np.nanmean(shot_noise_z[:,day:],axis=1),z)
plt.plot(np.nanmean(readout_noise_z[:,day:],axis=1),z)
plt.plot(np.nanmean(digi_noise_z[:,day:],axis=1),z)
plt.plot(np.nanmean(jpeg_noise_z[:,day:],axis=1),z)
plt.plot(np.nanmean(tot_noise_z[:,day:],axis=1),z,'--')
plt.legend(['shot','ro','digi','jpeg','tot'])
plt.title('Mean noise (day)')
plt.xlabel('LSB')
plt.ylabel('vertical pixel')

# %%
z = np.arange(0,nrows)
plt.plot(np.nanmean(rel_noise_z[:,:day]*100,axis=1),z)
plt.plot(np.nanmean(rel_noise_z[:,day:]*100,axis=1),z,'--')
plt.legend(['night','day'])
plt.title('Mean error')
plt.xlabel('relative error (%)')
plt.ylabel('vertical pixel')

# %%
