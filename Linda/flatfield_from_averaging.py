#%%
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l1_processing.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibration_functions import bin_image_with_BC
from mats_utils.rawdata.calibration import  calibrate_dataframe
from mats_utils.geolocation.altitude_correction import rows_to_altitudes
import datetime as DT
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from mats_utils.geolocation.coordinates import fast_heights ,heights
from lindas_own_functions import rename_CCDitem_entries, collapsandplot
from mats_utils.plotting.plotCCD import plot_image
from mats_utils.statistiscs.images_functions import create_imagecube
#from scipy.interpolate import CubicSpline
import math
import pickle
from scipy.integrate import quad
from scipy.interpolate import interp1d
#%%

def integrate_image_around_max(image, fixaltvec, deltay):
    """Returns integrated signal around max altitude"""
    # Initialize arrays to store the results
    max_altitude = np.zeros(image.shape[1])  # Maximum altitude for each column
    integrated_signal = np.zeros(image.shape[1])  # Integrated signal for each column

    # Iterate over each column
    for i in range(image.shape[1]):
        # Find the altitude range to fit the function to
        istart=image[:,i].argmax()-round(deltay/(fixaltvec[1]-fixaltvec[0]))
        if istart<0: istart=0
        istop=image[:,i].argmax()+round(deltay/(fixaltvec[1]-fixaltvec[0]))+1
        if istop>image.shape[0]: istop=image.shape[0]
        fitalt=fixaltvec[istart:istop]
        fitcol=image[istart:istop,i]


        # Perform spline interpolation
        interp_func = interp1d(fitalt, fitcol, kind='linear')

        # Generate interpolated values
        x_interp = np.linspace(fitalt[0], fitalt[-1], num=50)
        y_interp = interp_func(x_interp)
    
        # # Fit a polynomial of degree 5 to the column
        # coeffs = np.polyfit(fitalt, fitcol, deg=2)

        # # Define the polynomial function
        # poly_func = np.poly1d(coeffs)

        # # Find the altitude with the maximum signal
        # max_altitude[i] = -coeffs[1] / (2 * coeffs[0])


        max_altitude[i] = x_interp[y_interp.argmax()]


        # # Define the integration limits
        x = y_interp.argmax()  # Index of the maximum value
        # deltay
        x_lower = x - round(len(x_interp)/4)
        x_upper = x + round(len(x_interp)/4)


        # Perform the integration
        integrated_signal[i]= y_interp[x_lower:x_upper].mean()

        # plt.plot(y_interp[x_lower:x_upper],x_interp[x_lower:x_upper],'b')
        # plt.plot(fitcol,fitalt,'r') 

    return integrated_signal, max_altitude

#%%
starttime=DT.datetime(2023,4,10,0,0,0)
endtime=DT.datetime(2023,4,10,0,15,0)



#filter={'CCDSEL': [1, 4]} 
#filter={'CCDSEL': [1, 1], 'NRBIN': [1, 1], 'NCBINCCDColumns': [1, 1], 'NCOL':[2047,2048], 'NROW':[511,512]} 
#dfl1b = read_MATS_data(starttime, endtime,filter,level='1b',version='0.4')
filter_IR1={'CCDSEL': [1,1], 'NRBIN': [2, 2],'NCBINCCDColumns': [40, 40],'NCOL':[43,43], 'NROW':[187,187], 'NRSKIP':[109,109]} 
filter_IR2={'CCDSEL': [4,4], 'NRBIN': [2, 2],'NCBINCCDColumns': [40, 40],'NCOL':[43,43], 'NROW':[187,187], 'NRSKIP':[136,136]}
filter_IR3={'CCDSEL': [3,3], 'NRBIN': [6, 6],'NCBINCCDColumns': [200, 200],'NCOL':[8,8], 'NROW':[63,63], 'NRSKIP':[127,127]}
filter_IR4={'CCDSEL': [2,2], 'NRBIN': [6, 6],'NCBINCCDColumns': [200, 200],'NCOL':[8,8], 'NROW':[63,63], 'NRSKIP':[60,60]}

#dfl1a = read_MATS_data(starttime, endtime,filter, level='1a',version='0.5')

ir1l1a=read_MATS_data(starttime, endtime,filter_IR1, level='1a',version='0.5')
ir2l1a=read_MATS_data(starttime, endtime,filter_IR2, level='1a',version='0.5')
ir3l1a=read_MATS_data(starttime, endtime,filter_IR3, level='1a',version='0.5')
ir4l1a=read_MATS_data(starttime, endtime,filter_IR4, level='1a',version='0.5')

# ir1l1b=read_MATS_data(starttime, endtime,filter_IR1, level='1b',version='0.4')
# ir2l1b=read_MATS_data(starttime, endtime,filter_IR2, level='1b',version='0.4')
# ir3l1b=read_MATS_data(starttime, endtime,filter_IR3, level='1b',version='0.4')
# ir4l1b=read_MATS_data(starttime, endtime,filter_IR4, level='1b',version='0.4')




#%%

#df=copy.deepcopy(dfl1a)

calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'
instrument=Instrument(calibration_file)
#dfl1b=calibrate_dataframe(dfl1a, instrument)
ir1l1b=calibrate_dataframe(ir1l1a, instrument)
ir2l1b=calibrate_dataframe(ir2l1a, instrument)
ir3l1b=calibrate_dataframe(ir3l1a, instrument)
ir4l1b=calibrate_dataframe(ir4l1a, instrument)
#%%

#if len(ir1l1a)!=len(ir1l1b) or len(ir2l1a)!=len(ir2l1b) or len(ir3l1a)!=len(ir3l1b) or len(ir4l1a)!=len(ir4l1b):
if len(ir2l1a)!=len(ir2l1b) or len(ir3l1a)!=len(ir3l1b) or len(ir4l1a)!=len(ir4l1b):
#if len(ir2l1a)!=len(ir2l1b):
    raise ValueError('The length of the dataframes are not the same')




ir1l1b['IMAGE']=ir1l1a['IMAGE'].apply(np.fliplr) 
ir2l1b['IMAGE']=ir2l1a['IMAGE']
ir3l1b['IMAGE']=ir3l1a['IMAGE'].apply(np.fliplr)
ir4l1b['IMAGE']=ir4l1a['IMAGE']



#rename_CCDitem_entries(dfl1b)



# #%%
# df=copy.deepcopy(dfl1b)
# ir1=df[df.CCDSEL==1]
# ir2=df[df.CCDSEL==4]
# ir3=df[(df.CCDSEL==3)]
# ir4=df[(df.CCDSEL==2)]

#%%
# fig,ax=plt.subplots(2,1)
# CCDitem=ir1.iloc[0]
# plot_image(CCDitem, ax=ax[0], save=False)
# fixaltvec=np.arange(61000, 107000, 1000)
# rows_to_altitudes(CCDitem, fixaltvec)
# n=10
# ir1s=ir1l1b#[:n]
# ir2s=ir2l1b#[:n]
# ir3s=ir3l1b#[:n]
# ir4s=ir4l1b#[:n]



#imagefield='IMAGE'
fixaltvec=np.arange(61000, 107000, 1000)



ccdlist=[ir1l1b]#,ir2l1b,ir3l1b,ir4l1b]
#ccdlist=[ir1l1b]
meanimages={}
    

daynight='day_and_night'
for ccddf in ccdlist:
    if daynight=='daytime':
        ccddf=ccddf[ccddf['TPsza']<93.0]
    elif daynight=='nighttime':
        ccddf=ccddf[ccddf['TPsza']>99.0]       
    else:
        daynight='day_and_night' 

    ccddf.loc[:,'image_fixalt']= ccddf.apply(rows_to_altitudes, args=(fixaltvec, 'ImageCalibrated'), axis=1)
    ccddf.loc[:,'image_fixalt_nocal']= ccddf.apply(rows_to_altitudes, args=(fixaltvec, 'IMAGE'), axis=1)
    fig,ax=plt.subplots(5,1, figsize=(10,25))
    fig_nocal,ax_nocal=plt.subplots(5,1, figsize=(10,25))
    channel=ccddf.iloc[0]['channel']
    print(ccddf.iloc[0]['channel'])
    imagecube= create_imagecube(ccddf, image_specification='image_fixalt')
    imagecube_nocal= create_imagecube(ccddf, image_specification='image_fixalt_nocal')
    collapsandplot(imagecube,1, ax[0], signallabel='', title=channel+ ' calibrated '+daynight)
    collapsandplot(imagecube_nocal,1,ax_nocal[0], signallabel='',title=channel+' nocal '+daynight)

    with open('output/imagecube_'+channel, 'wb') as f:
        pickle.dump(imagecube, f)
    with open('output/imagecube_nocal_'+channel, 'wb') as f:
        pickle.dump(imagecube_nocal, f) 
    
    meanimage=imagecube.mean(axis=0)
    meanimage_nocal=imagecube_nocal.mean(axis=0)

    
    ax_nocal[1].plot(meanimage_nocal.argmax(axis=0), 'b', label='max')
    ax_nocal[1].set_title('row of max value')
    vertical_mean_nocal=meanimage_nocal.mean(axis=0)
    ax_nocal[2].plot(vertical_mean_nocal/vertical_mean_nocal.mean(), 'r', label='normalised mean')

    ax[1].plot(meanimage.argmax(axis=0), 'b', label='max')
    ax[1].set_title('row of max value')
    vertical_mean=meanimage.mean(axis=0)
    ax[2].plot(vertical_mean/vertical_mean.mean(), 'r', label='normalised mean')
   

    

    sp=plot_CCDimage(meanimage, fig, ax[3], title=channel, altvec=fixaltvec)
    sp_nocal=plot_CCDimage(meanimage_nocal, fig_nocal, ax_nocal[3], title=channel+' nocal', altvec=fixaltvec)
    #sp = ax[2].imshow(imagecube.mean(axis=0), cmap="magma", origin="lower", interpolation="none")
    ax[3].text(5,70000, 'starttime:'+str(starttime))
    ax[3].text(5,65000, 'endtime:'+str(endtime))
    ax_nocal[3].text(5,70000, 'starttime:'+str(starttime))
    ax_nocal[3].text(5,65000, 'endtime:'+str(endtime))



    if channel=='IR1' or channel=='IR2' or channel=='UV1' or channel=='UV2':
        downsizeratio=11
        arrshape=[meanimage.shape[0],meanimage.shape[1]//downsizeratio]
        meanimage_binned= np.mean(meanimage.reshape(-1,downsizeratio),axis=1).reshape(arrshape)
        meanimage_binned_nocal= np.mean(meanimage_nocal.reshape(-1,downsizeratio),axis=1).reshape(arrshape)
    else: #no need to bin
        meanimage_binned= meanimage
        meanimage_binned_nocal=meanimage_nocal

    for i in range(meanimage_binned.shape[1]):
        ax[4].plot(meanimage_binned[:,i], label='row '+str(i))
        #plt.ylim(400,550)
        #plt.xlim(19,28)
    ax[4].legend()
    ax[4].set_title('profiles calibrated')

    for i in range(meanimage_binned_nocal.shape[1]):    
        ax_nocal[4].plot(meanimage_binned_nocal[:,i], label='row '+str(i))
        #plt.ylim(400,550)
        #plt.xlim(19,28)
    ax_nocal[4].legend()
    ax_nocal[4].set_title('profiles not calibrated')
        





    deltay=10000
    integrated_signal, maxaltvec=integrate_image_around_max(meanimage, fixaltvec, deltay)    

    ax[2].plot(integrated_signal/integrated_signal.mean(), label='normalised mean signal over '+str(deltay)+' m at max altitude')
    ax[2].set_title('mean normalised signals calibrated')
    ax[2].legend()

    integrated_signal_nocal, maxalt_nocal=integrate_image_around_max(meanimage_nocal, fixaltvec, deltay)  
    ax_nocal[2].plot(integrated_signal_nocal/integrated_signal_nocal.mean(), label='normalised mean signal over '+str(deltay)+' m at max altitude')
    ax_nocal[2].set_title('mean normalised signals not calibrated')   
    ax_nocal[2].legend() 


    plt.tight_layout()




    fig.savefig('images/flatfield_av_cal_'+channel+daynight+'.png')
    fig_nocal.savefig('images/flatfield_av_nocal_'+channel+daynight+'.png')


    #Plot difference between cal and nocal
    figdiff,axdiff=plt.subplots(5,1, figsize=(10,20))

    meanimage_biassub=meanimage_nocal-291 #291 is an aproxamate value of TBLNK 
    sp=plot_CCDimage(meanimage/meanimage.mean(), figdiff, axdiff[0], title=channel+' cal', altvec=fixaltvec)
    sp=plot_CCDimage(meanimage_nocal/meanimage_nocal.mean(), figdiff, axdiff[1], title=channel+' nocal', altvec=fixaltvec)
    sp=plot_CCDimage(meanimage_biassub/meanimage_biassub.mean(), figdiff, axdiff[2], title=channel+' nocal but bias sub', altvec=fixaltvec) 
    sp=plot_CCDimage(meanimage/meanimage.mean()-meanimage_nocal/meanimage_nocal.mean(), figdiff, axdiff[3], title='cal-nocal '+channel, altvec=fixaltvec)
    sp=plot_CCDimage(meanimage/meanimage.mean()-meanimage_biassub/meanimage_biassub.mean(), figdiff, axdiff[4], title='cal-nocal(biassub) '+channel, altvec=fixaltvec) 
    axdiff[0].text(5,75000, daynight)
    axdiff[0].text(5,70000, 'starttime:'+str(starttime))
    axdiff[0].text(5,65000, 'endtime:'+str(endtime))
    
    
    fig.savefig('images/flatfield_av_cal_m_nocal'+channel+daynight+'.png')

    #meanimages[channel]=me

    # fig1,ax1=plt.subplots(3,1)
    # sp=plot_CCDimage(meanimage, fig1, ax1[0], title='cal '+channel)
    # sp=plot_CCDimage(meanimage_nocal, fig1, ax1[1], title='nocal '+channel)
    # sp=plot_CCDimage(meanimage/meanimage.mean()-meanimage_nocal/meanimage_nocal.mean(), fig1, ax1[2], title='diff '+channel)
    # plt.tight_layout()


#%%
# #%%
# namelist=["image_bias_sub","image_desmeared",
#     "image_dark_sub","image_calib_nonflipped","ImageCalibrated"]

# fig, ax=plt.subplots(5,1)
# for i, name in enumerate(namelist):
#     df.iloc[0][name].shape
#     plot_CCDimage(df.iloc[0][name], fig, ax[i], title=name)
   

# %%
