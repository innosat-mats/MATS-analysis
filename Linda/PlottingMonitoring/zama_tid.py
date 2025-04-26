#%%
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT
from mats_utils.plotting.plotCCD import simple_plot, plot_image, orbit_plot
from mats_utils.plotting.animate import generate_gif
from matplotlib import pyplot as plt
from database_generation.experimental_utils import plot_CCDimage



#Check available data with 
#aws s3 ls ops-payload-level1b-v0.6/2023/5/ --profile mats 

# data folder
data_folder = '/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/'

# times for start and stop
start_time = DT.datetime(2023, 3, 23, 14, 0, 0)
stop_time = DT.datetime(2023, 3, 23, 15, 0, 0)



# filter
filter={'CCDSEL': [5,6], 'TPlat':[-60, 60],  'channel':['IR1','NADIR']}
#'TPlon':[0, 40],
#%%
# read in measurements
df = read_MATS_data(start_time, stop_time,level='1b',version='0.9')


#%%
#df.iloc

for index, CCD in df06[:6].iterrows():
    fig, ax=plt.subplots(3,1)
    plot_CCDimage(df06.iloc[index].ImageCalibrated, axis=ax[0],fig=fig, title=df06.iloc[index].channel+' v0.6')
    unitfix=df05.iloc[index].TEXPMS/1000/10
    plot_CCDimage(df05.iloc[index].ImageCalibrated/unitfix, axis=ax[1],fig=fig, title=str(unitfix)+' * v0.5 '+df05.iloc[index].channel)
    diff=df06.iloc[index].ImageCalibrated-df05.iloc[index].ImageCalibrated/unitfix
    plot_CCDimage(diff, axis=ax[2],fig=fig, title='Difference')

    plt.savefig('../output/version05to06'+ df06.iloc[index].channel +'.png')
    
#%%
simple_plot(df06,data_folder)

#%% testing orbit_plot changes

orbit_plot(df06,data_folder+'reprocessingfeb/',nbins=7)

# %%
#/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/test2
for i in range(1,7):

    generate_gif('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/output/reprocessingfeb/CCDSEL'+str(i)+'/', data_folder+'orbit_feb_CCDSEL'+str(i)+'.gif')
# %%
