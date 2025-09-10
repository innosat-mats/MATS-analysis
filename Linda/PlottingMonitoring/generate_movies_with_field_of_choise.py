#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from mats_utils.imagetools.additional_fields import add_field_with_subtracted_rolling_mean, add_field_with_subtracted_rolling_mean2
import numpy as np

from database_generation.experimental_utils import plot_CCDimage


# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

readindata=True
level='1a' 
phenomenon='nightglow'
if readindata:
    if phenomenon=='airglow':
        # general data with airglow and aurora
        start_time = DT.datetime(2023, 2, 11, 12, 0, 0)
        stop_time = DT.datetime(2023, 2, 11, 15, 0, 0)
        filter={'CCDSEL': [1,4]}
        stringn='airglow'
    elif phenomenon=='dayglow':
        start_time = DT.datetime(2023, 2, 11, 12, 38, 0)
        stop_time = DT.datetime(2023, 2, 11, 13, 28, 0)
        filter={'CCDSEL': [1,4]}
        stringn='dayglow'
    elif phenomenon=='nightglow':
        start_time = DT.datetime(2023, 2, 11, 13, 48, 0)
        stop_time = DT.datetime(2023, 2, 11, 14, 1, 0)
        #start_time = DT.datetime(2023, 2, 11, 11, 58, 0)
        #stop_time = DT.datetime(2023, 2, 11, 12, 28, 0)
        filter={'CCDSEL': [1,4]}
        stringn='nightglow'       

    elif phenomenon=='nlc':

        start_time = DT.datetime(2023, 2, 2, 19, 35, 0)
        stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)

        # filter
        filter={'CCDSEL': [5,6]}
        stringn='nlc'

    # read in measurements

    #df = read_M
    df = read_MATS_data(start_time, stop_time,filter,level=level,version='1.0')
    name='df_'+stringn+'_'+level+'_'+start_time.strftime('%Y%m%d_%H%M%S')+'_'+stop_time.strftime('%Y%m%d_%H%M%S')
    df.to_pickle(data_folder+name+'.pkl')
else:
    airglowname='df_ir1ir220230211_120000_20230211_150000.pkl'
    airglowname='df_airglow_1b_20230211_120000_20230211_150000.pkl'

    nlcname='df_nlc_1a_20230202_193500_20230202_195000.pkl'
    nlcname='df_nlc_1b_20230202_193500_20230202_195000.pkl'

    dayglowname='df_nightglow_1b_20230211_134600_20230211_140100.pkl'

 
    if phenomenon == 'airglow':
        df = pd.read_pickle(data_folder+airglowname)
    elif phenomenon == 'dayglow':
        df = pd.read_pickle(data_folder+dayglowname)
    elif phenomenon == 'nightglow':
        df = pd.read_pickle(data_folder+dayglowname)
    elif phenomenon == 'nlc':
        df = pd.read_pickle(data_folder+nlcname)

#%%
#df = read_MATS_data(start_time, stop_time,level='1b',version='0.6')
if phenomenon in ['airglow', 'dayglow','nightglow']:
    channel='IR2'
else:   
    channel='UV2'
df_channel = df[df['channel'] == channel].copy()

#%%'
setting='IMAGE' # choose between 'diff', 'diff_from_first', 'diff_from_bg', 'rolling_mean' or 'ImageCalibrated' or if level 1a 'IMAGE'


if setting=='diff':
    # Create a field that is the difference between consecutive rows in df_channel['ImageCalibrated']
    df_channel['ImageCalibrated_diff'] = df_channel['ImageCalibrated'].diff()
    df_channel['ImageCalibrated_diff'].iloc[0]=df_channel['ImageCalibrated_diff'].iloc[1]
#if diff_from_first:
if setting=='diff_from_first':
    #Create a field that is the difference between df_channel['ImageCalibrated'] and the first row in df_channel['ImageCalibrated']
    FirstImage=df_channel['ImageCalibrated'].iloc[0]
    df_channel['ImageCalibrated_diff_from_first'] = df_channel['ImageCalibrated'].apply(lambda x: x - FirstImage)
if setting == 'diff_from_bg':
    #Create a field that is the difference between df_channel['ImageCalibrated'] and the first row in another dataframe
    start_time1 = DT.datetime(2023, 2, 2, 19, 35, 0)
    stop_time1 = DT.datetime(2023, 2, 2, 19, 37, 0)
    dfbg = read_MATS_data(start_time, stop_time,level=level,version='1.0')
    dfbg_channel = dfbg[dfbg['channel'] == channel].copy()
    bgImage=dfbg_channel.iloc[0].ImageCalibrated
    df_channel['ImageCalibrated_diff_from_bg'] = df_channel['ImageCalibrated'].apply(lambda x: x - bgImage)
if setting=='rolling_mean':
    add_field_with_subtracted_rolling_mean(df_channel, 'ImageCalibrated', 'ImageCalibrated_minus_rolling_mean', window_before=10, window_after=10, skipbefore=3, skipafter=3)

#%%
from matplotlib import pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(10, 10))
imgnr=100
if setting == 'ImageCalibrated' or setting == 'IMAGE':
    name=setting
else:
    name='ImageCalibrated_'+setting
if level == '1a':
    plot_CCDimage(df_channel['IMAGE'].iloc[imgnr], fig=fig, axis=ax[0], title='First image')
else:
    plot_CCDimage(df_channel['ImageCalibrated'].iloc[imgnr], fig=fig, axis=ax[0], title='First image')
plot_CCDimage(df_channel[name].iloc[imgnr], fig=fig, axis=ax[1], title='100th image')

#%% testing orbit_plot changes
#simple_plot(df,data_folder)

CCDs = df_channel[df_channel['CCDSEL'] == 4]
subfolder=setting+level+phenomenon
imagedir=data_folder+subfolder+'/'+channel+'/'
if phenomenon in ['airglow', 'dayglow', 'nightglow']:
    colorscheme='viridis'
else:
    colorscheme='Blues_r'
#orbit_plot(df_channel,imagedir,nbins=7, cmap='Blues_r', plothistogram=False, field_of_choise='ImageCalibrated_minus_rolling_mean', ranges=[-50,50]) 
if phenomenon == 'dayglow' and setting == 'ImageCalibrated':
    ranges=[0, 800]
    orbit_plot(df_channel,imagedir,nbins=7, cmap=colorscheme, plothistogram=False, printpositioninfo=False, ranges=ranges)
elif phenomenon == 'dayglow' and setting == 'diff':
    ranges=[-10, 10]
    orbit_plot(df_channel,imagedir,nbins=7, cmap=colorscheme, plothistogram=False, printpositioninfo=False, ranges=ranges, field_of_choise='ImageCalibrated_diff')
elif phenomenon == 'nightglow' and setting == 'ImageCalibrated':
    ranges=[0, 35]   
    orbit_plot(df_channel,imagedir,nbins=7, cmap=colorscheme, plothistogram=False, printpositioninfo=False, ranges=ranges)
elif phenomenon == 'nightglow' and setting == 'IMAGE':
    ranges=[300, 900]   
    orbit_plot(df_channel,imagedir,nbins=7, cmap=colorscheme, plothistogram=False, printpositioninfo=False, ranges=ranges)

else:
    orbit_plot(df_channel,imagedir,nbins=7, cmap=colorscheme, plothistogram=False, printpositioninfo=False)

print('Images are to be found in '+imagedir)
#cmap='YlGnBu_r'
# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2


if channel == 'IR1':
    generate_gif(imagedir+'CCDSEL1/', data_folder+'orbit_IR1_'+subfolder+'.gif')
elif channel == 'IR2':
    generate_gif(imagedir+'CCDSEL4/', data_folder+'orbit_IR2_'+subfolder+'.gif')
elif channel == 'UV2':
    generate_gif(imagedir+'CCDSEL6/', data_folder+'orbit_UV2_'+subfolder+'.gif')
print('Gif is to be found in '+data_folder+'orbit_'+channel+'_'+subfolder+'.gif')
# %%
