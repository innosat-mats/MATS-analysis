#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean
import numpy as np

from database_generation.experimental_utils import plot_CCDimage

def add_field_with_subtracted_rolling_mean2(df, field, outfieldname,  window_before=10, window_after=20, skipbefore=0, skipafter=0):


    def rolling_mean(images, window_before=10, window_after=20, skipbefore=0, skipafter=0):
        means = []
        for i in range(len(images)):
            # Define the window range
            start = max(0, i - window_before)
            end = min(len(images), i + window_after + 1)

            combined_range = list(range(start, i - skipbefore)) + list(range(i + skipafter, end))
            
            # Calculate the rolling mean for the window
            rolling_mean = np.mean([images[j] for j in combined_range], axis=0)
            means.append(rolling_mean)
        return means

    rolling_means = rolling_mean(df[field].tolist(), window_before=window_before, window_after=window_after, skipbefore=skipbefore, skipafter=skipafter)
    df[outfieldname] = [
        img - mean for img, mean in zip(df_channel[field], rolling_means)
    ]
    


    # rolling_means = rolling_mean(df_channel['ImageCalibrated'].tolist(), window=90, skipnr=10)
    # df_channel['ImageCalibrated_diff_from_earlier_mean'] = [
    #     img - mean for img, mean in zip(df_channel['ImageCalibrated'], rolling_means)

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
rolling_mean=True
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

if rolling_mean:

    add_field_with_subtracted_rolling_mean(df_channel, 'ImageCalibrated', 'ImageCalibrated_minus_rolling_mean', window_before=20, window_after=20, skipbefore=5, skipafter=5)




#%% testing orbit_plot changes
#simple_plot(df,data_folder)
subfolder='movie21'
imagedir=data_folder+subfolder+'/'+channel+'/'
orbit_plot(df_channel[320:400],imagedir,nbins=7, cmap='YlGnBu_r', plothistogram=False, field_of_choise='ImageCalibrated_minus_rolling_mean', ranges=[-50,50])   
#cmap='YlGnBu_r'
# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2

if channel == 'IR1':
    generate_gif(imagedir+'CCDSEL1/', data_folder+'orbit_IR1_'+subfolder+'.gif')
elif channel == 'IR2':
    generate_gif(imagedir+'CCDSEL4/', data_folder+'orbit_IR2_'+subfolder+'.gif')
elif channel == 'UV2':
    generate_gif(imagedir+'CCDSEL6/', data_folder+'orbit_UV2_'+subfolder+'.gif')
# %%
