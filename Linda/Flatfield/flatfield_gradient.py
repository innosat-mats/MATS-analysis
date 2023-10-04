# Plots the horizontal average of several images, flipping the filpped channels
#%% Import modules
import os
#os.chdir('/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda')
from mats_utils.rawdata.read_data import read_MATS_data
from lindas_own_functions import create_imagecube
from mats_utils.plotting.plotCCD import all_channels_plot
from mats_utils.plotting.animate import generate_gif
#from lindas_own_functions import collapsandplot
import pandas as pd
import datetime as DT
from mats_utils.geolocation import satellite as satellite
import matplotlib.pyplot as plt
import numpy as np


def plot_CCDimage_hmean(axis, image, title="", clim=999):
    yax = range(0, image.shape[0])
    sp = axis.plot(image.mean(axis=1), yax)
    axis.set_title(title)
    return sp

def meanfit(image, ax=None):

    #Correct for slant due to uneven lighting of screen
    x=range(0,image.shape[1])

    meanmean=image.mean()
    mean_v=image.mean(0)/meanmean

    coef = np.polyfit(x, mean_v, 1) 

    if ax:
        poly1d_fn = np.poly1d(coef)
        ax.plot(x,mean_v)
        #col=axnr.lines[-1].get_color()
        ax.plot(x, poly1d_fn(x))
    return coef,mean_v
   

def sza(EXPDate):
    (satlat, satlon, satLT, nadir_sza, nadir_mza,
    TPlat, TPlon,TPLT, TPsza, TPssa) = satellite.get_position(EXPDate)
    return TPsza

def flipchannel(CCDitem):
    if CCDitem['channel'] in ['IR1','IR3','UV1','UV2']:
        return np.fliplr(CCDitem.IMAGE)

    else:
        return CCDitem.IMAGE
        


#%%
# Select Time

start_time=DT.datetime(2023,1,30,18,0,0)
stop_time=DT.datetime(2023,1,30,18,30,0)

df_long=read_MATS_data(start_time, stop_time)

print(df_long.columns.tolist())


# %% Prepare df

df=df_long#[8000:12100]
df['TPsza']=df.apply(lambda row: sza(row['EXPDate']),axis=1)
df['IMAGE']=df.apply(lambda row: flipchannel(row), axis=1)
#df=df[(df.TPsza>0)]
print('length df: ',len(df))



#%%


daynight='night'
if daynight=='night':
    dfdaynight=df[(df.TPsza>95)]
elif daynight=='day':
    dfdaynight=df[(df.TPsza<85)]

#%%

# all_channels_plot(dfdaynight[:50], './output/test/', optimal_range=True)
# generate_gif('./output/test/ALL','output/film_'+ daynight +'.gif' )

#%%
channels=['IR1','IR2','IR3','IR4']
sigmode='HSM'
fig_coef, ax_coef=plt.subplots(6,1, figsize=(10,14))
for ind, channel in enumerate(channels):

    dfchannel=dfdaynight[dfdaynight.channel==channel]
    if channel in ['IR1','IR2']:
        dfchannel=dfchannel[df.NROW==140] #only select the ones that have a certain nr of rows
    if channel in ['IR3','IR4']:
        dfchannel=dfchannel[df.NROW==47] 
    if len(dfchannel)>0:

        imagecube=create_imagecube(dfchannel, image_specification='IMAGE')
        meanimage=imagecube[0:int(imagecube.shape[0]/2+1),:,:].max(0) # take the max value of the lower half
        # if channel in ['IR3','IR4']: #background channels
        #     meanimage=meanimage[25:33,1:8]
        # else:
        #     meanimage=meanimage[150:200,10:40]



        axnr=ax_coef[ind]
        coef, meanv=meanfit(meanimage, ax=axnr)

        axnr.set_title(channel+' '+daynight)
        axnr.text(0, 1.0,'coef :'+str(coef) )
        axnr.text(0, .95,'nr of images :'+str(len(dfchannel)) )
    else:
        print('Warning empty df for channel', channel)

plt.tight_layout()

"""
    col=axnr.lines[-1].get_color()
    x=np.arange(0, 2047)
    poly1d_fn = np.poly1d(coef_wo)
    axnr.plot(poly1d_fn(x), 'b')
    print(coef_wo)
    poly1d_fn = np.poly1d(coef_w)
    axnr.plot(poly1d_fn(x), 'r')
    axnr.set_xlim([400, 1600])
    axnr.set_ylim([0.97, 1.03])

    axnr.text(1100, 1.02,'coef_wo :'+str(coef_wo) )
    axnr.text(1100, 0.98,'coef_w :'+str(coef_w) )
    print(coef_w)
    plt.tight_layout()

    plt.legend()

"""
# %%
