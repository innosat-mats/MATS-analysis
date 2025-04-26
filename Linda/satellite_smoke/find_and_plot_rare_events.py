#%%
import requests
import numpy as np
import matplotlib.pyplot as plt
import pickle
import datetime as DT

import os
from mats_utils.plotting.animate import generate_gif 
from mats_utils.plotting.plotCCD import plot_image, orbit_plot
from sklearn.ensemble import IsolationForest
import pandas as pd
from database_generation.experimental_utils import plot_CCDimage
from mats_utils.rawdata.cropping import make_crop_filter
from mats_utils.selection_tools.select_at_random import select_random_images_all_channels



 #%%
# Function to download images from the MATS satellite
def download_image(url, save_path):
    response = requests.get(url)
    with open(save_path, 'wb') as f:
        f.write(response.content)








def make_model(df, field='ImageCalibrated', whatmodel='IsolationForest'):
    """
    Trains a model on the images in the given dataframe.

    Args:
        df (DataFrame): The dataframe containing the images.
        field (str): The column name of the images in the dataframe. Default is 'ImageCalibrated'.

    Returns:
        model (IsolationForest): The trained Isolation Forest model.
    """
    
    # Convert the images in df['ImageCalibrated'] to a 3D array
    images = np.stack(df[field].values)
    if whatmodel=='IsolationForest':
        # Reshape the images to 2D arrays
        reshaped_images = images.reshape(images.shape[0], -1)
        # Create an Isolation Forest model
        model = IsolationForest(contamination='auto')
        # Fit the model to the reshaped images
        model.fit(reshaped_images)
    else:
        print('No model was created, the input model was not recognized')
        model=None
    return model

def make_variance_model(df, field='ImageCalibrated', plot=False):
    """
    Calculates the standard deviation of variation across the field of view 
    at the bottom half of the images, for all images in the given dataframe. 
    A Gaussian function is fitted to the data, and the mean and standard deviation 
    are returned.


    Args:
        df (DataFrame): The dataframe containing the images.
        field (str): The column name of the images in the dataframe. Default is 'ImageCalibrated'.

    Returns:
        mu (float): The mean of the anomaly scores.
        std (float): The standard deviation of the anomaly scores.
    """
    from scipy.stats import norm
    images = np.stack(df[field].values)
    anomalyvalues=anomaly_calc(images)
    plt.hist(anomalyvalues, bins=50)
    # Fit a Gaussian function to the data
    mu, std = norm.fit(anomalyvalues)
    return mu, std

def anomaly_calc(images, maxaltpix=75):
    """
    Calculates the average standard deviation of signal strenght variation across the field of view 
    at the bottom half of the images up to a given altitude pixel. This 
    number is given ans an anolmaly score.
    Args:
        images (array): The images.
        maxaltpix (int): The maximum altitude pixel to consider. Default is 75.

    Returns:
        avstd_across (array or float): Anomaly scores for the images.
    """
    if len(images.shape)==2:
        avstd_across=np.mean(np.std(images[0:maxaltpix,:], axis=1), axis=0)
    elif len(images.shape)==3:
        avstd_across=np.mean(np.std(images[:,0:maxaltpix,:], axis=2), axis=1)# Average standard deviation across the field of view
    else:
        Exception('The input images have the wrong shape', images.shape)

    anomaly_score=avstd_across
    return anomaly_score

def rank_oddness_of_images(df, field, model=None):
    """
    Finds the most unlike images based on their anomaly scores.

    Args:
        df (DataFrame): The dataframe containing the images.
        field (str): The column name of the images in the dataframe.  The lower, the more abnormal.
    
    Returns:
        anomaly_scores (array): The anomaly scores for the images.
    """
    # Convert the images in df[field] to a 3D array
    images = np.stack(df[field].values)
    # Reshape the images to 2D arrays
    reshaped_images = images.reshape(images.shape[0], -1)

    # Predict the anomaly scores for each image. The lower, the more abnormal.
    if model is None:
    # make my own model, based on how much the field varies horizontally
        print('No model given, creating a new one')

    else: # use the model that was given as input
        anomaly_scores = model.decision_function(reshaped_images)
    
    # Create a new column 'anomaly_score' in df
    df.loc[:, 'anomaly_score'] = anomaly_scores

    return anomaly_scores

def rank_oddness_based_on_variance(df, field):
    import scipy.stats as stats

    for i in range(len(df)):
        image=df.iloc[i][field]
        anom_value=anomaly_calc(image)
        fraction_below = stats.norm.cdf(anom_value, loc=mu, scale=std)
        df.loc[i, 'anomaly_score'] = fraction_below
    anomaly_scores = df['anomaly_score'].values
    return anomaly_scores

def sort_images_by_anomaly_score(df, anomaly_scores):
    """  
    Function to sort the images based on their anomaly scores
    Args:
        df (DataFrame): The dataframe containing the images.
        anomaly_scores (array): The anomaly scores for the images.
    Returns:
        df_oddest_first (DataFrame): The dataframe containing the most unlike images.

    """
    sorted_indices = np.argsort(anomaly_scores)
    df_oddest_first = df.iloc[sorted_indices]
    return df_oddest_first



#%%
# Select random images to train the model on
if True:
    #channel='IR1'
    crop='CROPD'
    #filter_channelcrop=make_crop_filter(channel, crop)
    starttime = DT.datetime(2023, 2, 10, 0, 0, 0)
    endtime = DT.datetime(2023, 5, 10, 0, 0, 0)
    seed=42
    nrimages=4
    idiff=1
    dfchannelsdict=select_random_images_all_channels(starttime, endtime, nrimages,crop, seed, idifference=idiff)

#%%
pickle.dump(dfchannelsdict, open('testdata/df_random_allchannels_'+crop+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'wb'))

 #%%
# Load the data
#Select how many images to skip between the images where the difference is taken    
idiff=1
#Select which channel to serch for anomalies in 
channel='IR1'
nrimages=1600
seed=42
crop='CROPD'
with open('testdata/df_random_allchannels_'+crop+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'.pkl', 'rb') as f:

    dfchannelsdict= pickle.load(f)

df=dfchannelsdict[channel]



#%%
# Select which field to search for anomalies in and cut the dataset
field='ImageCalibratedDiff'+str(idiff)
#field='ImageCalibrated'
#df['MeanValue'] = df['ImageCalibrated'].apply(lambda x: np.mean(x))
#Remove data from South Atlantic Anomaly
df = df[~(((df['satlat'] <= 0) & (df['satlat'] >= -60) & ((df['satlon'] > 300) | (df['satlon'] < 30))))]

df_day = df[(df['TPsza'] <= 90)]
df_night = df[(df['TPsza'] >= 100)]




df_sel=df_day.copy()


#%%

# dfchannelsdictcut = {}
# for ichannel in channels:
#     dfchannelsdictcut[ichannel] = dfchannelsdict[ichannel].loc[df_sel.index]
useAI=False
if useAI:
    model=make_model(df_sel,field, whatmodel='IsolationForest')
    #use model to predict anomaly scores
    anomaly_scores=rank_oddness_of_images(df_sel,field, model=model)
    plt.hist(anomaly_scores, bins=50)
    anomalymethod='IsolationForest'
else:
    mu, std =make_variance_model(df_sel,field)
    anomaly_scores=rank_oddness_based_on_variance(df_sel,field)
    anomalymethod='Variance'


    #plt.hist(anomaly_scores, bins=50)
#%%
df_oddfirst=sort_images_by_anomaly_score(df_sel, anomaly_scores)

# for ichannel in channels:
#     dfchannelsdict_oddfirst[ichannel]=sort_images_by_anomaly_score(dfchannelsdictcut[ichannel], anomaly_scores)

#df_oddfirst['MeanValue'] = df_oddfirst['ImageCalibrated'].apply(lambda x: np.mean(x))


    
#%%
# Create the directory and plot images 
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

directory = os.path.join(data_folder, 'satellite_smoke'+anomalymethod+'_'+channel+'_idiff_'+str(idiff)+'_nimg_'+str(nrimages)+'_seed_'+str(seed)+'_'+crop)  # Create a directory to save the images in
os.makedirs(directory, exist_ok=True)



# Select corresponding rows in all the channels, including the one that was selected as odd
channels=['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']#, 'NADIR']
dfchannelsdict_oddfirst = {}        
for ichannel in channels:
    dfchannelsdict_oddfirst[ichannel]=dfchannelsdict[ichannel].loc[df_oddfirst.index] #reorder so that the image datat that has the oddest looking image in the selected channel becomes the first image in the dataframe


# Plot the most odd looking images in the selected channel, and the corresponding images in the other channels
for index in range(0, 5): #upper limit is how many of the oddest images you want to plot
    #plot_image(CCD, outpath=directory+'/')8790
    #CCD = df_oddfirst.iloc[index]
    fig, ax = plt.subplots(6, 2, figsize=(7, 15))  
    for i, ichannel in enumerate(channels):
        iCCD=dfchannelsdict_oddfirst[ichannel].iloc[index]
        if not np.isnan(iCCD.year): #If the data is not empty, all values are set to Nan if no image is found for that period
            
            plot_CCDimage(iCCD.ImageCalibrated,fig=fig, axis=ax[i,0], title=ichannel+'Time: '+str(iCCD['EXPDate']))
            ax[i,0].text(0.05, 0.9, f'Latitude: {iCCD.satlat}', transform=ax[i,0].transAxes, color='white')
            ax[i,0].text(0.05, 0.8, f'Longitude: {iCCD.satlon}', transform=ax[i,0].transAxes, color='white')
   
            
            plot_CCDimage(iCCD.ImageCalibratedDiff1,fig=fig, axis=ax[i,1], title='Differnce')













#%%

generate_gif(directory+'/', data_folder+'satellite_smoke_'+'orbit.gif')





# %%

values=['MeanValue','TPsza','TPssa']

for value in values:
    plt.figure()
    plt.plot(df_oddfirst.anomaly_score, df_oddfirst[value],'.')
    plt.xlabel('Anomaly score')
    plt.ylabel(value)

# plt.plot(df_oddfirst.anomaly_score, df_oddfirst.MeanValue,'.')
# plt.plot(df_oddfirst.anomaly_score, df_oddfirst.TPssa,'.')
# plt.plot(df_oddfirst.anomaly_score, df_oddfirst.TPsza,'.')
# %%

import matplotlib.pyplot as plt

times_IR1 = pd.to_datetime(dfchannelsdict['IR1']['EXPDate'])
times_IR2 = pd.to_datetime(dfchannelsdict['IR2']['EXPDate'])

plt.plot(times_IR1, label='IR1')
plt.plot(times_IR2[:len(times_IR1)], label='IR2')
plt.xlabel('Index')
plt.ylabel('Time')
plt.legend()
plt.show()

same_times = (times_IR1 == times_IR2)
print(f"Are the times the same? {same_times}")
# %%

filter_IR1={'TPsza':[97,150],'CCDSEL': [1,1], 'NRBIN': [2, 2],'NCBINCCDColumns': [40, 40],'NCOL':[43,43], 'NROW':[187,187], 'NRSKIP':[109,109]} 
crop='CROPD'
filter_channelcrop=make_crop_filter(channel, crop)
starttime = DT.datetime(2023, 2, 10, 0, 0, 0)
endtime = DT.datetime(2023, 5, 10, 0, 0, 0)