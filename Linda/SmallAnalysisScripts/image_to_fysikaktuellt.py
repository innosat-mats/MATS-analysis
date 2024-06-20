#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif

from database_generation.experimental_utils import plot_CCDimage

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

# times for start and stop
#start_time = DT.datetime(2023, 7, 18, 6, 0, 0)
#stop_time = DT.datetime(2023, 7, 18, 6, 1, 0)
start_time = DT.datetime(2023, 2, 2, 19, 38, 0)
stop_time = DT.datetime(2023, 2, 2, 19, 50, 0)
start_time1 = DT.datetime(2023, 2, 2, 19, 35, 0)
stop_time1 = DT.datetime(2023, 2, 2, 19, 37, 0)

# filter
filter={'CCDSEL': [5, 6]}

#%%
# read in measurements
df = read_MATS_data(start_time, stop_time,filter,level='1b')
df_UV1=df[df['CCDSEL']==5]
df_UV2=df[df['CCDSEL']==6]
#%%

start_time1 = DT.datetime(2023, 2, 2, 19, 35, 0)
stop_time1 = DT.datetime(2023, 2, 2, 19, 37, 0)


df_bg = read_MATS_data(start_time1, stop_time1,filter,level='1b')
df_bg_UV1=df_bg[df_bg['CCDSEL']==5]
df_bg_UV2=df_bg[df_bg['CCDSEL']==6]

#%%
#df.iloc
from matplotlib import pyplot as plt

preImage=df_bg_UV2.iloc[0].ImageCalibrated
for index, CCD in df_UV2[32:36].iterrows():
    print(index)
    if index == 1:
        #preImage=CCD.ImageCalibrated
        hej=0

    else:
        fig, ax=plt.subplots(1)
        diffimg=CCD.ImageCalibrated-preImage
        sp = ax.imshow(diffimg, cmap="Blues", origin="lower", interpolation="none")
        
        # set colormap to blues
        sp.set_cmap('Blues')
        # set size of figure
        fig.set_size_inches(9,3)
        ax.set_aspect("auto")
        ax.set_title('MATS image of Noctilucent Clouds')
        plt.ylabel( 'Vertical view [pixels]', fontsize=12)

        plt.xlabel( 'Horizontal view [pixels]', fontsize=12)
        fig.savefig(data_folder+'Fysikakturellt_'+str(index)+'.png')
        print('timestamp:',CCD.EXPDate)

        
        #preImage=CCD.ImageCalibrated

        #plot_CCDimage(CCD.ImageCalibrated, title='CCDSEL'+str(CCD.CCDSEL), borders=True, nrsig=3)

        # wr
        # remove colorbar




#simple_plot(df,data_folder)

#%% testing orbit_plot changes

orbit_plot(df,data_folder+'test2/',nbins=7)

# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2
generate_gif('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2/CCDSEL5/', data_folder+'orbit.gif')
# %%
