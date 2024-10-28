#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import datetime as DT

import os
import pandas as pd
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.selection_tools.select_at_random import select_random_images_all_channels
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.rawdata.calibration import calibrate_dataframe
from mats_l1_processing.instrument import Instrument
import glob

print(' Linda TODO: remove TBLNK bias for instance by running the calibration instead of removing it manually')
print('Linda TODO: cut on signal strength, see histograms at the bottom here')

    

#%%
def find_darkest_area_of_CCD(image, irowstart,irowend,icolstart,icolend):
# Function that takes the image, divides it into 7x7 areas and finds the darkest area
# The darkest area is defined as the area with the lowest mean value
# The function returns the mean value of the darkest area
    # Define the size of the areas

  
    if (icolend-icolstart)<10: #background channel
        npix_vertical = 2
    else:
        npix_vertical = 10
    npix_horizontal = 5
    # Define the number of areas in each direction
    n_ver = int(np.floor((icolend-icolstart) / 5))
    n_hor = int(np.floor((irowend-irowstart) / 5))
    if n_ver==0 or n_hor==0:
        Exception('Too few pixels in the image')

    # Create an array to store the mean values of the areas
    mean_values = np.zeros([n_ver, n_hor])
    # Loop over the areas
    for i in range(n_ver):
        for j in range(n_hor):
            # Calculate the mean value of the area
            mean_values[i, j] = np.mean(image[i * npix_vertical:(i + 1) * npix_vertical, j * npix_horizontal:(j + 1) * npix_horizontal])
    # Find the darkest area 
    darkest_area = np.unravel_index(np.argmin(mean_values), mean_values.shape)
    # Return the mean value of the darkest area
    return mean_values[darkest_area]




goodpointing=False
if goodpointing:
    crop='CROPD'
    #filter_channelcrop=make_crop_filter(channel, crop)
    starttime = DT.datetime(2023, 2, 10, 0, 0, 0)
    endtime = DT.datetime(2023, 5, 10, 0, 0, 0)
else:
    #Nonpointing
    #start_time = DT.datetime(2024, 1, 12, 9, 0, 0)
    #stop_time = DT.datetime(2024, 1, 12, 9, 30, 0)
    starttime = DT.datetime(2023, 10, 1, 1, 0, 0)
    endtime = DT.datetime(2024, 10, 1, 1, 0, 0)
    #starttime = DT.datetime(2024, 1, 12, 9, 0, 0)
    #endtime = DT.datetime(2024, 1, 12, 9, 30, 0)
    #dfnonpoint = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)
seed=42
nrimages=1000

readfrom='database'
if readfrom=='aws':
    dfchannelsdict=select_random_images_all_channels(starttime, endtime, nrimages, level='1a')
    os.makedirs('../output/testdata', exist_ok=True)
    if goodpointing:
        pickle.dump(dfchannelsdict, open('../output/testdata/df_random_allchannels_'+crop+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'wb'))
    else:
        pickle.dump(dfchannelsdict, open('../output/testdata/df_random_allchannels_badpointing_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'wb'))

elif readfrom=='pickledaws':
    # read in from pickle
    if goodpointing:
        nrimages=100
        dfchannelsdict = pickle.load(open('../output/testdata/df_random_allchannels_'+crop+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'rb'))
    else:
        dfchannelsdict = pickle.load(open('../output/testdata/df_random_allchannels_badpointing_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'rb'))

elif readfrom=='database':
    # Read in the data from downloaded files
    # Note that these files are not always all channels at the same time

    directory = '/Users/lindamegner/MATS/MATS-retrieval/data/every_few_days/'
    pkl_files = glob.glob(os.path.join(directory, "*.pkl"))
    dataframes = [pd.read_pickle(file) for file in pkl_files]
    combined_df = pd.concat(dataframes, ignore_index=True)
    len(combined_df)
    dfsel=combined_df[(combined_df['TPsza']>105) & (combined_df['afsTangentH_wgs84']>120)] 
    dfchannelsdict={}
    for channel in ['IR1','IR2','IR3','IR4','UV1','UV2']:
        dfchannelsdict[channel]=dfsel[dfsel['channel']==channel]    



#%%

if not 'instrument' in locals():
    instrument = Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_MATSinstrument.toml')
    #instrument= Instrument('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml')
#%%
channel='IR1'
df=dfchannelsdict[channel]
#plot heater temqeratu§
#re temperature, HTR8B against each other
fig, ax=plt.subplots(1)
ax.scatter(df['temperature'], df['HTR8B'])
plt.xlabel('Temperature of UV2 CCD housing [C]')
plt.ylabel('Temperature of UV1 CCD housing [C]')
plt.show()
# # %%
# starttime = DT.datetime(2023, 4, 7, 21, 48, 56)
# endtime=starttime+DT.timedelta(seconds=7)
# df= read_MATS_data(starttime, endtime, level='1b', pritfilesys=False)

# %%
# 
imagefield='IMAGE'
#imagefield='image_bias_sub'
#imagefield='image_desmeared'

if imagefield=='IMAGE':
    bias=293 #Bias is 293 for every readout
else:
    bias=0

channels=['IR1','IR2','IR3','IR4']#,'UV1','UV2']
#channels=['IR2']
# plot the topmean against the heater temperature
#fig, ax=plt.subplots(4,1, figsize=(3,13))

fig_gumbel, ax_gumbel=plt.subplots(1)
fig, ax=plt.subplots(1)
for inr, channel in enumerate(channels):
    df=dfchannelsdict[channel]
    df['temperature']=df['HTR8A']   

    #df=calibrate_dataframe(df, instrument, debug_outputs=True)
    #select TPsza>98, afsTangentH_wgs84°>100, (TPssa<90 or TPssa>270) TPSSa is the angle between the satellite TP pointing and the sun
    df=df[(df['TPsza']>98) & (df['afsTangentH_wgs84']>90)] #& ((df['TPssa']<90) | (df['TPssa']>270))]
    df=df[(df['TPssa']<90)]
    #df=df[(df['TPsza']>98) & (df['afsTangentH_wgs84']>100) & ((df['TPssa']<80) | (df['TPssa']>100))]

    if channel=='IR1':
        irowstart=175#230
        irowend=180#2a40
    elif channel=='IR2':
        irowstart=173
        irowend=176
        #taped
        #irowstart=216
        #irowend=217
    elif channel=='IR3':
        irowstart=55
        irowend=59
    elif channel=='IR4':
        irowstart=50
        irowend=55
    elif channel=='UV1' or channel=='UV2':
        Exception(' not implemented yet')

    #irowstart=50
    #irowend=60
    icolstart=int(df['NCOL'].iloc[0]/5)
    icolend=int(df['NCOL'].iloc[0]/5*4)
    df['topmean']= df[imagefield].apply(lambda x: np.mean(x[irowstart:irowend,icolstart:icolend]))


    #df['topmean'] =df[imagefield].apply(lambda x: find_darkest_area_of_CCD(x, 0, x.shape[0], 0, x.shape[0]))



    cutvalute=df['topmean'].mean()+3*df['topmean'].std()
    df=df[df['topmean']<cutvalute]
    # Remove data when df['temperature'] is NaN
    df = df.dropna(subset=['temperature'])
    df['nbins']=(df['NRBIN']*df['NCBINCCDColumns']*df['NCBINFPGAColumns'])
    df['topmean']=(df['topmean']-bias)/df['TEXPMS']*1000/df['nbins']

    df["CCDunit"] =instrument.get_CCD(channel)
    df['denominator']=df["CCDunit"].apply(lambda x: x.calib_denominator('High'))

    df['topmean_in_Gumbel']=df['topmean']/df['denominator']
    #ax[inr].scatter(df['temperature'], df['topmean'])




    ax.scatter(df['temperature'], df['topmean'], label=channel)
    
    
    # Fit a polynomial of degree 1 to the data

    coeffs = np.polyfit(df['temperature'], df['topmean'], 1)
    # Define the polynomial function
    poly_func = np.poly1d(coeffs)
    #plot the polynomial
    x = np.linspace(df['temperature'].min(), df['temperature'].max(), 100)
    y = poly_func(x)
    ax.plot(x, y, label=f'{poly_func[1]:.2f}x + {poly_func[0]:.2f}')


    ax_gumbel.scatter(df['temperature'], df['topmean_in_Gumbel'], label=channel)
    # Fit a polynomial of degree 1 to the data
    coeffs = np.polyfit(df['temperature'], df['topmean_in_Gumbel'], 1)
    # Define the polynomial function
    poly_func = np.poly1d(coeffs)
    #plot the polynomial
    x = np.linspace(df['temperature'].min(), df['temperature'].max(), 100)
    y = poly_func(x)
    ax_gumbel.plot(x, y, label=f'{poly_func[1]:.2f}x + {poly_func[0]:.2f}')



ax.set_xlabel('temperature [C]')
ax.set_ylabel('Topmean [Counts/s/unbinndedpixel]')
ax.legend()
ax_gumbel.set_xlabel('temperature [C]')
ax_gumbel.set_ylabel('Topmean [10**12 ph/s/m2/sr]')
ax_gumbel.legend()

  
#%%


df=df[(df['topmean']>0.5)]

for image in df[imagefield]:
    plot_CCDimage(image)

    #ax[inr].set_title(channel)

#%%
image=df[imagefield].iloc[0]
plot_CCDimage(image,title='image', clim=[0, 8000])
#
# %%

#nihgtglow
start_time = DT.datetime(2023, 2, 20, 19, 47, 0)
stop_time = DT.datetime(2023, 2, 20, 19, 58, 15)
start_time = DT.datetime(2024, 1, 12, 9, 0, 0)
stop_time = DT.datetime(2024, 1, 12, 9, 30, 0)
dfnonpoint = read_MATS_data(start_time,stop_time,version='0.7',level='1a',dev=False)
# %%
df=dfnonpoint[dfnonpoint['channel']=='IR1']
# %%
# plot the topmean against TPssa

fig, ax=plt.subplots(1)
ax.scatter(df['TPssa'], df['topmean'])
plt.xlabel('TPssa [deg]')
plt.ylabel('Topmean [Counts/s/unbinndedpixel]')
plt.show()
# %%
#plot historgram of mean of imagefield

channels=['IR1','IR2','IR3','IR4']#,'UV1','UV2']

for channel in channels:
    df=dfchannelsdict[channel]
    plt.hist(df[imagefield].apply(np.mean), bins=50)
    plt.xlabel('Mean of Image Field')
    plt.ylabel('Frequency')
    plt.title('Histogram of Mean of Image Field for '+channel)
    plt.show()
# %%
