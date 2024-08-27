#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
import numpy as np

from database_generation.experimental_utils import plot_CCDimage


# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

readindata=False
if readindata:
    # general data with airglow and aurora
    start_time = DT.datetime(2023, 2, 11, 12, 0, 0)
    stop_time = DT.datetime(2023, 2, 11, 15, 0, 0)
    #and NLC
    #start_time = DT.datetime(2023, 2, 2, 19, 35, 0)
    #stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)

    # filter
    #filter={'CCDSEL': [1,4]}

    # read in measurements

    df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
    #df = read_MATS_data(start_time, stop_time,filter,level='1b',version='0.6')
    name='df_ir1ir2'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(data_folder+name+'.pkl')
else:
    airglowname='df_ir1ir220230211_120000_20230211_150000.pkl'
    nlcname='df_NLC20230202_193500_20230202_195000.pkl'
    df = pd.read_pickle(data_folder+airglowname)

#%%
#df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
channel='IR2'
df_channel = df[df['channel'] == channel].copy()

#%%'
diff=True
diff_from_earlier_mean=True
diff_from_first=True
diff_from_bg=False
if diff:
    # Create a field that is the difference between consecutive rows in df_channel['ImageCalibrated']
    df_channel['ImageCalibrated_diff'] = df_channel['ImageCalibrated'].diff()
if diff_from_first:
    #Create a field that is the difference between df_channel['ImageCalibrated'] and the first row in df_channel['ImageCalibrated']
    FirstImage=df_channel['ImageCalibrated'].iloc[0]
    df_channel['ImageCalibrated_diff_from_first'] = df_channel['ImageCalibrated'].apply(lambda x: x - FirstImage)
if diff_from_bg:
    #Create a field that is the difference between df_channel['ImageCalibrated'] and the first row in another dataframe
    start_time1 = DT.datetime(2023, 2, 2, 19, 35, 0)
    stop_time1 = DT.datetime(2023, 2, 2, 19, 37, 0)
    dfbg = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
    dfbg_channel = dfbg[dfbg['channel'] == channel].copy()
    bgImage=dfbg_channel.iloc[0].ImageCalibrated
    df_channel['ImageCalibrated_diff_from_bg'] = df_channel['ImageCalibrated'].apply(lambda x: x - bgImage)

if diff_from_earlier_mean:
    # Create a field that is the difference between df_channel['ImageCalibrated'] and the mean of the previous 10 rows
    def rolling_mean(images, window, skipnr=0):
        means = []
        for i in range(len(images)):
            if i < window:
                means.append(np.mean(images[:i+1], axis=0))
            elif (i  > window and i < skipnr+window):
                means.append(np.mean(images[window:i+1], axis=0))
            else:
                means.append(np.mean(images[i-skipnr-window+1:i+1], axis=0))
        return means

    rolling_means = rolling_mean(df_channel['ImageCalibrated'].tolist(), window=10, skipnr=10)
    df_channel['ImageCalibrated_diff_from_earlier_mean'] = [
        img - mean for img, mean in zip(df_channel['ImageCalibrated'], rolling_means)
    ]

# if diff_from_earlier_meanx:
#     #Create a field that is the difference between df_channel['ImageCalibrated'] and the mean of the previous 10 rows
#     df_channel['ImageCalibrated_diff_from_earlier_mean'] = df_channel['ImageCalibrated'].rolling(window=10).mean
#     df_channel['ImageCalibrated_diff_from_earlier_mean'] = df_channel['ImageCalibrated'] - df_channel['ImageCalibrated_diff_from_earlier_mean'] 



#%% testing orbit_plot changes
#simple_plot(df,data_folder)
imagedir=data_folder+'movie17/'+channel+'/'
orbit_plot(df_channel[300:350],imagedir,nbins=7, cmap='YlGnBu_r', plothistogram=False, field_of_choise='ImageCalibrated_diff_from_earlier_mean')
#cmap='YlGnBu_r'
# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2

if channel == 'IR1':
    generate_gif(imagedir+'CCDSEL1/', data_folder+'orbit_IR1.gif')
elif channel == 'IR2':
    generate_gif(imagedir+'CCDSEL4/', data_folder+'orbit_IR2.gif')
elif channel == 'UV2':
    generate_gif(imagedir+'CCDSEL6/', data_folder+'orbit_UV2.gif')
# %%
