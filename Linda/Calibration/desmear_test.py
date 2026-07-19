#%%




from mats_l1b_tools.error import add_flags
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

#%%

# # #%% Select on explicit time NLC
start_time = DT.datetime(2023, 3, 14, 19, 38, 0) 
stop_time = DT.datetime(2023, 3, 14, 19, 40, 0)
#start_time = DT.datetime(2023, 2, 19, 18, 59, 39) 
#stop_time = DT.datetime(2023, 2, 20, 19, 3, 44)
#between 2023 03 22 09:15:00 UTC - 2023 03 22 09:25:00 UTC
#start_time = DT.datetime(2023, 3, 22, 9, 15, 0)
#stop_time = DT.datetime(2023, 3, 22, 9, 25, 0)



#start_time = DT.datetime(2023, 4, 22, 9, 15, 0)
#stop_time = DT.datetime(2023, 4, 22, 10, 45, 0)

#2023_4_16_1_11_0-2023_4_16_4_11_0 testdata from Lukas This has some dayglow, and some non-SAA nightglow (with 256 flag)
#start_time = DT.datetime(2023, 4, 16, 1, 11, 0)
#stop_time = DT.datetime(2023, 4, 16, 4, 11, 0)


# # #%% Select on explicit time nadir
#start_time = DT.datetime(2023, 4, 6, 0, 0, 0) 
#stop_time = DT.datetime(2023, 4, 6, 0, 30, 0)

#%%
pickle_filename = f"../output/df1a_{start_time.strftime('%Y%m%d_%H%M%S')}_to_{stop_time.strftime('%Y%m%d_%H%M%S')}.pkl"

#%%
# try to read_from_pickle, if it fails read from raw data and save to pickle
try:
    df_loc_l1a = pickle.load(open(pickle_filename, 'rb'))
    print('Level 1a read from pickle')
except:
    df_loc_l1a = read_MATS_data(start_time,stop_time,version='1.0',level='1a',dev=False)
    df_loc_l1a.to_pickle(pickle_filename)
    print('Level 1a read and saved to pickle')
    


#%%
#try reading the calibrated dataframe from pickle, if it fails, save it to pickle
try:
    df_loc = pickle.load(open(pickle_filename.replace('1a','1acal'), 'rb'))
    print('Calibrated dataframe read from pickle')
except:
    if not 'instrument' in locals():
        instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')
    df_loc=calibrate_dataframe(df_loc_l1a, instrument, debug_outputs=True)
    df_loc.to_pickle(pickle_filename.replace('1a','1acal'))
    print('Calibrated dataframe read and saved to pickle')
print('Level 1a calibrated')

#%%
#save calibrated dataframe to cvs
# calibrated_filename = pickle_filename.replace('1a','1acal').replace('.pkl','.csv')
# df_loc.to_csv(calibrated_filename)
# print('Calibrated dataframe saved to csv, filename: ', calibrated_filename)


#%%


#df_loc.to_pickle(pickle_filename.replace('1a','1acal'))

#print('Calibrated dataframe saved to pickle, filename: ', pickle_filename.replace('1a','1acal'))   
#%%

#load calibrated file
#calibratedfile = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/df1acal_20230219_185939_to_20230220_190344.pkl'

#df_loc = pickle.load(open(calibratedfile, 'rb'))

#select only fist 6000 rows of data 
#Select south atlantic anomaly
#df_loc = df_loc[(df_loc.satlon >-60) & (df_loc.satlon < -20) & (df_loc.satlat > -50) & (df_loc.satlat < -10)]

#%%
#Select nighttime data with SZA > 100



ir1_loc=df_loc[df_loc.channel=='IR1']
ir2_loc=df_loc[df_loc.channel=='IR2']
ir3_loc=df_loc[df_loc.channel=='IR3']
ir4_loc=df_loc[df_loc.channel=='IR4']
uv1_loc=df_loc[df_loc.channel=='UV1']
uv2_loc=df_loc[df_loc.channel=='UV2']
nadir_loc=df_loc[df_loc.channel=='NADIR']




#%%
#Plot and compare the two versions
channel = 'IR2'
if channel == 'IR2':
    df = ir2_loc
elif channel == 'IR1':
    df = ir1_loc
elif channel == 'IR3':
    df = ir3_loc
elif channel == 'IR4':
    df = ir4_loc
elif channel == 'UV1':
    df = uv1_loc
elif channel == 'UV2':
    df = uv2_loc
elif channel == 'NADIR':
    df = nadir_loc
else:
    raise ValueError('Channel not recognized')

#%%
#Flag
#    if 'CalibrationErrors' not in ds_slice.variables:
    #     raise KeyError('CalibrationErrors not available')
    # calibration_errors = ds_slice['CalibrationErrors'].astype(np.uint16)
    # # Create a new data variable for each flag
    # ds_slice['FlagBadColumns'] = (calibration_errors & 0b0000000000000001) > 0
    # ds_slice['FlagSingleEvent'] = (calibration_errors & 0b0000000000000010) > 0
    # ds_slice['FlagHotPixel'] = (calibration_errors & 0b0000000000000100) > 0
    # ds_slice['FlagNoHotPixel'] = (calibration_errors & 0b0000000000001000) > 0
    # ds_slice['FlagNegativeBias'] = (calibration_errors & 0b0000000000010000) > 0
    # ds_slice['FlagNonlinearCorrection'] = (calibration_errors & 0b0000000000100000) > 0
    # ds_slice['FlagSaturatedPixel'] = (calibration_errors & 0b0000000001000000) > 0
    # ds_slice['FlagDesmearNegative'] = (calibration_errors & 0b0000000010000000) > 0
    # ds_slice['FlagDesmearNoAtmosphere'] = (calibration_errors & 0b0000000100000000) > 0


df['FlagDesmearNoAtmosphere'] = df['CalibrationErrors'].apply(
    lambda calibration_errors: (np.asarray(calibration_errors, dtype=np.uint16) & 0b0000000100000000) > 0
)   

#Check if desmear no atmosphere flag is set anywhere in the image (then it should be set for the whole image)
df['FlagDesmearNoAtmosphereSetAnywhere'] = df['FlagDesmearNoAtmosphere'].apply(lambda flag_array: np.any(flag_array))

#%%
# # Find negative values in image_desmeared
# df_negatives = df[df['image_desmeared'].apply(lambda x: np.mean(x < 0))]
# print(f"Number of entries with negative values in image_desmeared: {len(df_negatives)}")


#%%

#select only entries with desmear no atmosphere flag set
df_desmear_no_atmosphere = df[(df['FlagDesmearNoAtmosphereSetAnywhere'] == True) ] #Also select only entries with SZA < 120, as we expect the desmear no atmosphere flag to be set for high SZA values where the atmosphere is very thin.
print(f"Number of entries with desmear no atmosphere flag set: {len(df_desmear_no_atmosphere)} out of {len(df)} total entries for channel {channel}")
fig, ax = plt.subplots(2,2, figsize=(10,10))

# panel 1: plot the SZA distribution for all entries
ax[0,0].hist(df.TPsza, bins=20)
ax[0,0].set_title(f'{channel} all {len(df)} entries')
ax[0,0].set_xlabel('Tangent point SZA (degrees)')
ax[0,0].set_ylabel('Count')
# panel 2: plot the SZA distribution for entries with desmear no atmosphere flag set
ax[0,1].hist(df_desmear_no_atmosphere.TPsza, bins=20)
ax[0,1].set_title(f'{channel} entries with desmear no atmosphere flag set')
ax[0,1].set_xlabel('Tangent point SZA (degrees)')
ax[0,1].set_ylabel('Count')

#Find what caused the desmear no atmosphere flag to be set for these entries, by comparing the image_desmeared to the image_linear 
# (which is the image prior to desmearing).
#'AnyIncrease' (np.any((image_desmear/image_linear) > 1.)) 
# '75PercentReduction' (np.mean(image_desmear)/np.mean(image_linear) < 0.25)
# 'Other reason' (neither of the above)

df_desmear_no_atmosphere['DesmearNoAtmosphereReason'] = df_desmear_no_atmosphere.apply(lambda row: 'AnyIncrease' if (np.any((row['image_desmeared']/row['image_linear']) > 1.)) else '75PercentReduction' if (np.mean(row['image_desmeared'])/np.mean(row['image_linear']) < 0.25) else 'Other reason', axis=1)
#plot histogram of the reasons for the desmear no atmosphere flag being set
#panel 3: plot the distribution of reasons for the desmear no atmosphere flag being set
ax[1,0].hist(df_desmear_no_atmosphere.DesmearNoAtmosphereReason, bins=20)
ax[1,0].set_title(f'{channel} reasons for desmear flag being set')
ax[1,0].set_xlabel('Reason')
ax[1,0].set_ylabel('Count')

plt.tight_layout()
plt.savefig('../output/desmear_no_atmosphere_analysis_'+channel+'.png')

#%%
# loop through the dataframe and plot the image_desmeared for entries with the desmear no atmosphere flag set, and compare to image_linear
inc=len(df_desmear_no_atmosphere)//8
nrimages=6
maxrows = nrimages*inc
endrow = min(maxrows, len(df_desmear_no_atmosphere))
for index, row in df_desmear_no_atmosphere.iloc[0:endrow:inc].iterrows():
   
    fig, axs = plt.subplots(3, 1, figsize=(10, 10))        
    plot_CCDimage(row['image_linear'], title=row.channel+ ' prior desmear '+str(row.TMHeaderTime), fig=fig, axis=axs[0])
    plot_CCDimage(row['image_desmeared'], title=row.channel+ ' desmeared '+str(row.TMHeaderTime), fig=fig, axis=axs[1])
    #difference = row['image_linear'] - row['image_desmeared']
    plot_CCDimage(row['image_linear'] - row['image_desmeared'], title=row.channel+ ' desmeared - linear '+str(row.TMHeaderTime), fig=fig, axis=axs[2])  
    plt.tight_layout()


#plot 10 examples of entries with the desmear no atmosphere flag set

fig, ax = plt.subplots(nrimages,1, figsize=(5,10))
i=0
for index, row in df_desmear_no_atmosphere[0:endrow:inc].iterrows():
    i=i+1
    if i>nrimages:
        break
    plot_CCDimage(row['image_linear'], title=row.channel+ ' prior desmear SZA'+ str(row.TPsza)+ ' '+str(row.TMHeaderTime), fig=fig, axis=ax[i-1])

    
plt.tight_layout()
plt.savefig('../output/desmear_no_atmosphere_examples_'+channel+'.png')



# %%
df['image_linear_median'] = df['image_linear'].apply(np.median)
df['image_desmeared_median'] = df['image_desmeared'].apply(np.median)
#plot linear and desmeared median values, as a function of time
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.plot(df.TMHeaderTime, df.image_linear_median, label='linear median')
ax.plot(df.TMHeaderTime, df.image_desmeared_median, label='desmeared median')
#plot only south atlantic anomaly
#df_loc = df_loc[(df_loc.satlon >-60) & (df_loc.satlon < -20) & (df_loc.satlat > -50) & (df_loc.satlat < -10)]
dfSSA=df[(df.satlon >-90) & (df.satlon < -20) & (df.satlat > -50) & (df.satlat < -0)]
ax.plot(dfSSA.TMHeaderTime, dfSSA.image_linear_median, label='linear median SSA', marker='o')
ax.plot(dfSSA.TMHeaderTime, dfSSA.image_desmeared_median, label='desmeared median SSA', marker='o')


ax.set_title(f'{channel} median values over time')
ax.set_xlabel('Time')
ax.set_ylabel('Median value')
ax.legend()
plt.tight_layout()

#plot satlat and satlon for the entries with desmear 
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.plot(df.TMHeaderTime, df.satlat, label='satlat')
ax.plot(df.TMHeaderTime, df.satlon, label='satlon')
ax.set_title(f'{channel} satellite latitude and longitude over time')
ax.set_xlabel('Time')
ax.set_ylabel('Degrees')
ax.legend()
plt.tight_layout()
# %%
