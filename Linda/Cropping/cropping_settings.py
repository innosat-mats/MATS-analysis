
#%% 


#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on Mon Sep 19 21:00:16 2022

@author: lindamegner

This script is used to spit out the cropping settings based on common field of view. 
Be aware of crappy coding though! :) 
# crop field 
This script is based on the old script test_alignment_launch_star.py

"""



from database_generation.experimental_utils import plotCCDitem
from mats_l1_processing.read_in_functions import read_CCDitems,read_CCDdata
from matplotlib import pyplot as plt
import numpy as np
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
from database_generation.experimental_utils import calibrate_CCDitems
import datetime as DT
from mats_utils.selection_tools.itemselect import select_on_time as seltime
from mats_utils.rawdata.time_tools import add_datetime as add_datetime
from mats_utils.geolocation import satellite as satellite
from mats_utils.imagetools.imagetools import shift_image
from mats_utils.plotting.plotCCD import orbit_plot, simple_plot
import pandas as pd
from mpldatacursor import datacursor
import copy

import warnings
#from mats_l1_processing.L1b_calibration_functions import shift_image
#%%
def plot_CCDimage(image, fig, axis, title="", clim=999, aspect="auto"):
    sp = axis.imshow(image, cmap="magma", origin="lower", interpolation="none")
    # sp=axis.pcolormesh(image, , cmap='viridis')
    if clim == 999:
        [col, row]=image.shape
        #Take the mean and std of the middle of the image, not boarders
        mean = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].mean()
        std = image[int(col/2-col*4/10):int(col/2+col*4/10), int(row/2-row*4/10):int(row/2+row*4/10)].std()
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    axis.set_aspect(aspect)
    return sp

def plot_CCDimage_transp(image, fig, axis, title="", clim=999, aspect="auto", alpha=1, 
    ncshift=0, nrshift=0, ncol=0, nrow=0, pos00=[0,0], colour='black'):
    sp = axis.imshow(image, cmap="viridis", origin="lower", interpolation="none", alpha=alpha)
    # sp=axis.pcolormesh(image, , cmap='viridis')
    if clim == 999:
        mean = image.mean()
        std = image.std()
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    #fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    rectangle = plt.Rectangle((ncshift, nrshift), ncol, nrow, facecolor='none', ec=colour)
    axis.add_patch(rectangle)
    #rectangle = plt.Rectangle((pos00[1], pos00[0]), 2047, 511, facecolor='none', ec=colour)
    #axis.add_patch(rectangle)
    axis.set_aspect(aspect)
    
    return sp

def nskip_to_nshift(ncskip, nrskip, ncol, pos00, channel):

    if channel=='IR1' or channel=='IR3' or channel=='UV1' or channel=='UV2':
        ncskiplr=2047-(ncol+ncskip)
    else:
        ncskiplr=ncskip
        
    ncshift=pos00[1]+ncskiplr
    nrshift=pos00[0]+nrskip

    return ncshift, nrshift

def nshift_to_nskip(ncshift, nrshift, ncol, pos00, channel):
    ncskiplr=ncshift-pos00[1]
    nrskip=nrshift-pos00[0]

    if channel=='IR1' or channel=='IR3' or channel=='UV1' or channel=='UV2':
        ncskip=2047-(ncol+ncskiplr)
    else:
        ncskip=ncskiplr

    return ncskip, nrskip    

class instrumentdata:
    def __init__(self):
        # class channelinfo:
        #     def __init__(self, x_pos, y_pos):
        #         self.x_pos=x_pos
        #         self.y_pos=y_pos

        
        # #Lindas settings as measured in the lab
        #self.x_pos = {'IR1':75,'IR2':-144,'IR3':-37,'IR4':0,'UV1':-88,'UV2':-156,'NADIR':0} 
        #self.y_pos = {'IR1': -47,'IR2':-76,'IR3':-66,'IR4':0,'UV1':-13,'UV2':-192,'NADIR':0} 
        #Donals settings from star measuremnets
        self.x_pos = {'IR1':83.47120648,'IR2':-141.29390513,'IR3':-46.07793043,'IR4':0,'UV1':-90.50990343,'UV2':-161.31804504,'NADIR':0} 
        self.y_pos = {'IR1': 48.92274474,'IR2':75.78736229,'IR3':67.06758131,'IR4':0,'UV1':5.37669736,'UV2':188.22050731,'NADIR':0}        
        for item in self.x_pos:
            self.x_pos[item]=round(self.x_pos[item])
        for item in self.y_pos:
            self.y_pos[item]=round(self.y_pos[item])


        # self.info_IR1=channelinfo(x_pos=-75,y_pos=47)
        # self.info_IR2=channelinfo(x_pos=144,y_pos=76)
        # self.info_IR3=channelinfo(x_pos=37,y_pos=66)
        # self.info_IR4=channelinfo(x_pos=0,y_pos=0)
        # self.info_UV1=channelinfo(x_pos=88,y_pos=13)
        # self.info_UV2=channelinfo(x_pos=156,y_pos=192) 



        #Stetting from analysis where baffle vinjetting start occuring
        #xmin = {'IR1': 20, 'IR2': 35,'IR3': 30,'IR4': 260,'UV1': 140,'UV2': 230,'NADIR': 0}
        #xmax = {'IR1': 1790, 'IR2': 2010,'IR3': 2010,'IR4': 2010,'UV1': 2010,'UV2': 2010,'NADIR': 0}
        xminread = {'IR1': 21, 'IR2': 26,'IR3': 26,'IR4': 216,'UV1': 140,'UV2': 271,'NADIR': 0}
        xmaxread = {'IR1': 1769, 'IR2': 2019,'IR3': 2019,'IR4': 2019,'UV1': 2019,'UV2': 2019,'NADIR': 0}


        flipped_channels = ['IR1','IR3','UV1','UV2']
        self.xmin=copy.deepcopy(xminread)
        self.xmax=copy.deepcopy(xmaxread)
        for ichannel in flipped_channels:
            self.xmin[ichannel]=2047-xmaxread[ichannel]
            self.xmax[ichannel]=2047-xminread[ichannel]

        #from lego analysis
        self.ymin = {'IR1': 180, 'IR2': 218,'IR3': 222,'IR4': 69,'UV1': 137,'UV2': 0,'NADIR': 0}
        self.ymax = {'IR1': 511, 'IR2': 511,'IR3': 511,'IR4': 434,'UV1': 511,'UV2': 511,'NADIR': 511}
    #def info_reset(channel, pos00):


#%%
def get_crop_positions(instr,pos00, channel):

    nrowmax=511
    ncolmax=2047
    ncol=ncolmax
    nrow=nrowmax

    ncskiplr=0
    nrskip=0
    
    if channel=='IR1':
        ncskiplr=-instr.x_pos['IR2']+instr.x_pos['IR1']
        ncol=ncolmax+(instr.x_pos['IR2']-instr.x_pos['IR1'])
        nrskip=instr.y_pos['IR1']-instr.y_pos['IR4']
        nrow=nrowmax-(instr.y_pos['IR2']-instr.y_pos['IR4'])
    elif channel=='IR2':
        ncskiplr=0
        ncol=ncolmax+(instr.x_pos['IR2']-instr.x_pos['IR1'])
        nrskip=instr.y_pos['IR2']-instr.y_pos['IR4']
        nrow=nrowmax-(instr.y_pos['IR2']-instr.y_pos['IR4'])   
    elif channel=='IR3':
        ncskiplr=-instr.x_pos['IR2']+instr.x_pos['IR3']
        ncol=ncolmax+(instr.x_pos['IR2']-instr.x_pos['IR1'])
        nrskip=instr.y_pos['IR3']-instr.y_pos['IR4']
        nrow=nrowmax-(instr.y_pos['IR2']-instr.y_pos['IR4'])
    elif channel=='IR4':
        ncskiplr=-instr.x_pos['IR2']+instr.x_pos['IR4']
        ncol=ncolmax+(instr.x_pos['IR2']-instr.x_pos['IR1'])
        nrskip=0
        nrow=nrowmax-(instr.y_pos['IR2']-instr.y_pos['IR4'])
    elif channel=='UV1':
        ncskiplr=-instr.x_pos['UV2']+instr.x_pos['UV1']
        ncol=ncolmax+(instr.x_pos['UV2']-instr.x_pos['UV1'])
        nrskip=0
        nrow=nrowmax-(instr.y_pos['UV2']-instr.y_pos['UV1'])
    elif channel=='UV2':
        ncskiplr=0
        ncol=ncolmax+(instr.x_pos['UV2']-instr.x_pos['UV1'])
        nrskip=instr.y_pos['UV2']-instr.y_pos['UV1']
        nrow=nrowmax-(instr.y_pos['UV2']-instr.y_pos['UV1'])
    
    CCDitem['ncskiplr_only_pos']=ncskiplr
    CCDitem['ncol_only_pos']=ncol
    CCDitem['nrskip_only_pos']=nrskip
    CCDitem['nrow_only_pos']=nrow

    cut_on_signal_study=True
    cut_on_lego_study=True
    if cut_on_signal_study: # cut in the x-direction
        if ncskiplr<instr.xmin[channel]:
            print('changing ncskiplr for '+channel+' :'+ str(ncskiplr))
            diff=instr.xmin[channel]-ncskiplr
            ncol=ncol-diff
            ncskiplr=instr.xmin[channel]
            print('new:  '+str(ncskiplr))
            
        if (ncskiplr+ncol)>instr.xmax[channel]:
            print('changing ncol for '+channel+' :'+ str(ncol))
            ncol=instr.xmax[channel]-ncskiplr
            print('new:  '+str(ncol))

    if cut_on_lego_study:    # cut in the y-direction
        if nrskip<instr.ymin[channel]:
            print('lego changing nrskip for '+channel+' :'+ str(ncskiplr))
            diff=instr.ymin[channel]-nrskip
            nrow=nrow-diff
            nrskip=instr.ymin[channel]
            print('new:  '+str(nrskip))
        if (nrskip+nrow)>instr.ymax[channel]:
            print('lego changing nrow for '+channel+' :'+ str(nrow))
            nrow=instr.ymax[channel]-nrskip
            print('new:  '+str(nrow))
            

# Add error margines:
    add_err_mag=False
    if add_err_mag:
        errmag=100
        if channel=='IR1' or channel=='IR3' or channel=='IR4':
            ncskiplr=ncskiplr-errmag
            ncol=ncol+errmag
        elif channel=='IR2':
            ncol=ncol+errmag



    ncshift=pos00[1]+ncskiplr #shift is in common fullframe ncskiplr is the cordinate of each ccd 
    #but always starting from the left of the image (ie needs to be changed for flipped channels)
    nrshift=pos00[0]+nrskip





    return ncshift, nrshift, ncskiplr,nrskip,ncol, nrow

def set_common_cropping(ncshift_max, nrshift_max, ncend_min, nrend_min, pos00, channel):
    ncol=ncend_min-ncshift_max
    nrow=nrend_min-nrshift_max
    ncshift=ncshift_max
    nrshift=nrshift_max
    ncskip, nrskip=nshift_to_nskip(ncshift, nrshift, ncol, pos00, channel)
    return ncshift, nrshift, ncskip,nrskip,ncol, nrow



def binned_to_unbinned(NCOL,NROW,NCBIN,NCBINFPGA, NRBIN):

    ncol=(NCOL+1)*NCBIN*2**NCBINFPGA
    nrow=NROW*NRBIN
    return ncol, nrow

def specify_settings(channel, cropversion):
 
    if cropversion=='CROPF':
        if channel=="UV1": 
            NRSKIP=5
            NRBIN=2
            NROW=162
            NCSKIP=201
            NCBIN=40
            NCOL=43-1
            NCBINFPGA=0

        if channel=="UV2":
            NRSKIP=188
            NRBIN=2
            NROW=162
            NCSKIP=271
            NCBIN=40
            NCOL=43-1
            NCBINFPGA=0

        if channel=="IR1": 
            NRSKIP=49
            NRBIN=2
            NROW=217
            NCSKIP=27
            NCBIN=40
            NCOL=43-1
            NCBINFPGA=0

        if channel=="IR2": 
            NRSKIP=76
            NRBIN=2
            NROW=217
            NCSKIP=75
            NCBIN=40
            NCOL=43-1
            NCBINFPGA=0

        if channel=="IR3": 
            NRSKIP=67
            NRBIN=6
            NROW=73
            NCSKIP=156
            NCBIN=215
            NCOL=8-1
            NCBINFPGA=0

        if channel=="IR4": 
            NRSKIP=0
            NRBIN=6
            NROW=73
            NCSKIP=216
            NCBIN=215
            NCOL=8-1
            NCBINFPGA=0



    elif cropversion=='CROPD': #det vi kört i januari till april 20223
        if channel=='UV1':
            NRSKIP=65
            NRBIN=2
            NROW=132
            NCSKIP=201
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='UV2':
            NRSKIP=247
            NRBIN=2
            NROW=132
            NCSKIP=271
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR1':
            NRSKIP=109
            NRBIN=2
            NROW=187
            NCSKIP=27
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR2': 
            NRSKIP=136
            NRBIN=2
            NROW=187
            NCSKIP=75
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR3':
            NRSKIP=127
            NRBIN=6
            NROW=63
            NCSKIP=156
            NCBIN=200
            NCOL=8
            NCBINFPGA=0

        elif channel=='IR4':
            NRSKIP=60
            NRBIN=6
            NROW=63 
            NCSKIP=216
            NCBIN=200
            NCOL=8
            NCBINFPGA=0


    elif cropversion=='CROP_TO_BOTTOM':
        if channel=='UV1':
            NRSKIP=0 #65
            NRBIN=2
            NROW=164 #132
            NCSKIP=201
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='UV2':
            NRSKIP=0 #247
            NRBIN=2
            NROW=256 #132
            NCSKIP=271
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR1':
            NRSKIP=0 #109
            NRBIN=2
            NROW=242 #187
            NCSKIP=27
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR2': 
            NRSKIP=0 #136
            NRBIN=2
            NROW=255 # 87
            NCSKIP=75
            NCBIN=40
            NCOL=43
            NCBINFPGA=0

        elif channel=='IR3':
            NRSKIP=0 #127
            NRBIN=6
            NROW=84 #63
            NCSKIP=156
            NCBIN=215
            NCOL=8
            NCBINFPGA=0

        elif channel=='IR4':
            NRSKIP=0 #60
            NRBIN=6
            NROW=73 #63 
            NCSKIP=216
            NCBIN=215
            NCOL=8
            NCBINFPGA=0


    ncol, nrow=binned_to_unbinned(NCOL,NROW,NCBIN,NCBINFPGA, NRBIN)
    pos00=CCDitem['set_pos00']
    ncskip=NCSKIP
    nrskip=NRSKIP
    ncshift, nrshift=nskip_to_nshift(ncskip, nrskip, ncol, pos00, channel)


    return ncshift, nrshift, ncskip,nrskip,ncol, nrow

#%%
########################################################################
#                        Program starts herehere                       #
########################################################################
star=7

instrument_analysis='/Users/lindamegner/MATS/retrieval/git/instrument_analysis/'
if star==999: 
    RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_November/'
else:
    RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacFiles_out_all/'
#RacOut='/Users/lindamegner/MATS/retrieval/FlightData/commissioning/RacOut_from25nov/'
#calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'
calibration_file='/Users/lindamegner/MATS/MATS-retrieval/MATS-analysis/Linda/calibration_data_linda.toml'


run_from_path=instrument_analysis+'Operational_analysis/commisioning/Lindas/'
image_path=run_from_path+'images/'


#%%
_,df = read_CCDdata(RacOut)
df_all=df.copy()


# %%

df=df_all


if star==4:
    date1 = DT.datetime(2022,11,29,8,30,00)
    date2 = DT.datetime(2022,11,29,10,30,00)
elif star==5:
    date1 = DT.datetime(2022,11,29,10,30,00)
    date2 = DT.datetime(2022,11,29,12,00,00)
elif star==1:
    date1 = DT.datetime(2022,11,29,12,00,00)
    date2 = DT.datetime(2022,11,29,13,00,00)
elif star==6:
    date1 = DT.datetime(2022,11,29,13,00,00)
    date2 = DT.datetime(2022,11,29,15,00,00)
elif star==7000:
    date1 = DT.datetime(2022,12,12,21,00,00)
    date2 = DT.datetime(2022,12,12,23,00,00)
elif star==7:
    date1 = DT.datetime(2022,12,5,19,50,00)
    date2 = DT.datetime(2022,12,5,20,30,00)
elif star==999: #full frame staring mode test LBFT
    date1 = DT.datetime(2022,11,22,14,19,00)
    date2 = DT.datetime(2022,11,23,14,51,00)

dfdatetime=add_datetime(df)
df = seltime(date1,date2,dfdatetime) #filtered dataframe between 11 and 12 UTC the 24th november 2022


CCDitems = read_CCDitems(RacOut,items=df.to_dict('records'))

# #%%
# CCDitemsdf = pd.DataFrame.from_dict(CCDitems)
# simple_plot(CCDitemsdf, image_path, nstd=2, cmap='inferno', custom_cbar=False,
#     ranges=[0, 1000], format='png')

#%%
image_specification='IMAGE'
for CCDitem in CCDitems:
    fig , ax= plt.subplots(1, 1, figsize=(8, 2))
    image=CCDitem[image_specification]
    sp=plot_CCDimage(image, fig, ax, title=CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000))

#%%
calibrate=True
if calibrate:
    calibrate_CCDitems(CCDitems, Instrument(calibration_file))
    for CCDitem in CCDitems:        
        totbin=CCDitem['NCBIN FPGAColumns'] *CCDitem['NCBIN CCDColumns']*CCDitem['NRBIN']
        CCDitem['image_calibrated']=CCDitem['image_calibrated']/totbin
        #Shift image, i.e. put image on common field of view
        image_common_fov, error_flags_flipnshift = shift_image(CCDitem, CCDitem['image_calibrated'])
        CCDitem['image_common_fov']=image_common_fov
        #pos00=np.argwhere(np.isfinite(CCDitem['image_common_fov']))[0]
        #CCDitem['pos00']=pos00
if calibrate:
    signallabel='Signal [10^10*ph/cm2/str/nm]'
    image_specification='image_calibrated'
else:
    signallabel='Counts'
    image_specification='IMAGE'

#%%



instr=instrumentdata()

fig1, ax1 = plt.subplots(6, 1)
fig2, ax2 = plt.subplots(1)
fig3, ax3 = plt.subplots(6, 1)
start=0
stop=5
fix_clim=[-10,100]
colours=['orange', 'blue', 'red', 'yellow','cyan','magenta', 'white']
for ind, CCDitem in enumerate(CCDitems):
    colour=colours[ind]
    plot_CCDimage(CCDitem['image_calibrated'], fig1, ax1[ind], CCDitem['channel'])    
    plot_CCDimage(CCDitem['image_common_fov'], fig3, ax3[ind], title=CCDitem['channel'], clim=fix_clim)  

    pos00=np.argwhere(np.isfinite(CCDitem['image_common_fov']))[0]
    ncshift, nrshift, ncskiplr,nrskip,ncol, nrow=get_crop_positions(instr,pos00, CCDitem['channel'])

    CCDitem['set_ncshift']=ncshift
    CCDitem['set_nrshift']=nrshift
    CCDitem['set_ncskiplr']=ncskiplr
    CCDitem['set_nrskip']=nrskip
    CCDitem['set_ncol']=ncol
    CCDitem['set_nrow']=nrow
    CCDitem['set_pos00']=pos00
    
    #Set the area to the common area
    if CCDitem['channel']=='IR1':
        ncend_min_ir=CCDitem['set_ncshift']+CCDitem['set_ncol']
    if CCDitem['channel']=='IR4':
        ncshift_max=CCDitem['set_ncshift']
    if CCDitem['channel']=='UV2':
        ncend_min=CCDitem['set_ncshift']+CCDitem['set_ncol']
        nrend_min_uv=CCDitem['set_nrshift']+CCDitem['set_nrow']
    #if CCDitem['channel']=='IR3':
        #nrshift_max=CCDitem['set_nrshift']
    if CCDitem['channel']=='IR4':
        nrend_min_ir=CCDitem['set_nrshift']+CCDitem['set_nrow']  
        nrshift_IR4_only_pos=pos00[0]+CCDitem['nrskip_only_pos']
        nrshift_max=nrshift_IR4_only_pos
 


    sp=plot_CCDimage_transp(CCDitem['image_common_fov'], fig2, ax2, 
    title='All channels', clim=fix_clim, alpha=0.3, ncshift=ncshift, nrshift=nrshift, 
        ncol=ncol, nrow=nrow, pos00=pos00, colour=colour)    
    
    ax2.text(pos00[1], pos00[0], str(CCDitem['channel']) , color=colour)
    rectangle = plt.Rectangle((ncshift, nrshift), ncol, nrow, facecolor='none', ec=colour)
    #ax2.add_patch(rectangle)
    ax3[ind].add_patch(rectangle)
    

    
    fig2.savefig('images/Alignment_crop_together_'+str(star)+'.png',  dpi=700) 
    fig3.savefig('images/Alignment_crop_all_'+str(star)+'.png',  dpi=700)   
  
#%%


fig4, ax4 = plt.subplots(1)
for ind, CCDitem in enumerate(CCDitems):
    # if CCDitem['channel']=='IR1' or CCDitem['channel']=='IR3' or CCDitem['channel']=='UV1' or CCDitem['channel']=='UV2':
    #     CCDitem['set_ncskip']=2047-(CCDitem['set_ncol']+CCDitem['set_ncskiplr'])
    # else:
    #     CCDitem['set_ncskip']=CCDitem['set_ncskiplr']
    # ncshift, nrshift= nskip_to_nshift(ncskip, nrskip, ncol, pos00, CCDitem['channel'])
    # print (CCDitem['channel']+' NCSKIP: '+str(CCDitem['set_ncskip'])
    #     +' NRSKIP: '+str(CCDitem['set_nrskip'])
    #     +' NCOL: '+str(CCDitem['set_ncol'])
    #     +' NROW: '+str(CCDitem['set_nrow']))



    # ncshift, nrshift= nskip_to_nshift(ncskip,nrskip, CCDitem['set_ncol'], CCDitem['set_pos00'], CCDitem['channel'])

    nrow=CCDitem['set_nrow']
    ncol=CCDitem['set_ncol']
    ncshift=CCDitem['set_ncshift']
    nrshift=CCDitem['set_nrshift']
    pos00=CCDitem['set_pos00']



    set_common_cropping_bool=False
    if set_common_cropping_bool:
        ##set_common_cropping(CCDitem[])
        if (CCDitem['channel']=='UV1' or CCDitem['channel']=='UV2'):
            nrend_min=nrend_min_uv
        else:
            nrend_min=nrend_min_ir
        ncshift, nrshift, ncskip,nrskip,ncol, nrow=set_common_cropping(ncshift_max, nrshift_max, ncend_min, nrend_min, pos00, CCDitem['channel'])
    else:
        # specify the used setting
        cropversion='CROPF'
        ncshift, nrshift, ncskip,nrskip,ncol, nrow= specify_settings(CCDitem['channel'],cropversion=cropversion)
       


    print ('Zero:'+CCDitem['channel']
        +' NCSKIP: '+str(ncskip)
        +' NRSKIP: '+str(nrskip)
        +' nr of col: '+str(ncol)
        +' NROW: '+str(nrow))

    print ('One:'+CCDitem['channel']
        +' NCSKIP: '+str(ncskip)
        +' NRSKIP: '+str(nrskip)
        +' nr oc col: '+str(ncol)
        +' NROW: '+str(nrow)
        +' NCSHIFT: '+str(ncshift)
        +' NRSHITF: '+str(nrshift)
        )





    colour=colours[ind]
    
    
    sp=plot_CCDimage_transp(CCDitem['image_common_fov'], fig4, ax4, 
    title=cropversion, clim=fix_clim, alpha=0.3, ncshift=ncshift, nrshift=nrshift, 
        ncol=ncol, nrow=nrow, pos00=pos00, colour=colour)    

    ax4.text(pos00[1], pos00[0], str(CCDitem['channel']) , color=colour)
    rectangle = plt.Rectangle((ncshift, nrshift), ncol, nrow, 
    facecolor='none', ec=colour)
    ax4.add_patch(rectangle)
    #ax3[ind].add_patch(rectangle)
  
fig4.savefig('images/cropping_alltogether.png',  dpi=700) 
  

#  rectangle = plt.Rectangle((10,10), 200, 200, fc='blue',ec="red")
  #  plt.gca().add_patch(rectangle) 

# %%
# Set current settings and to see how then look:

