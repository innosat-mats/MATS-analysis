#%% Import modules
#%matplotlib qt5
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, timezone
import warnings
import boto3
import pyarrow.dataset as ds
from pyarrow import fs
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from typing import (
    cast, Any, Dict, List, Optional, Sequence, SupportsFloat, Tuple, Union,
)
import numpy as np
import pyarrow as pa  # type: ignore
import pyarrow.dataset as ds  # type: ignore
import pyarrow.parquet as pq  # type: ignore
from pandas import DataFrame, Timestamp  # type: ignore
from PIL import Image
import time as time

pd.set_option('display.max_rows', 500)

#%%

# defining some parameters
channels = {1:'IR1',2:'IR4',3:'IR3',4:'IR2',5:'UV1',6:'UV2',7:'NADIR'}
CCDSEL = {'IR1':1,'IR2':4,'IR3':3,'IR4':2,'UV1':5,'UV2':6,'NADIR':7}
sampling_rates={'IR1':timedelta(seconds=6),
                    'IR2':timedelta(seconds=6),
                    'IR3':timedelta(seconds=6),
                    'IR4':timedelta(seconds=6),
                    'UV1':timedelta(seconds=6),
                    'UV2':timedelta(seconds=6),
                    'NADIR':timedelta(seconds=2)}

Port_dic = {'IR1' : ['CPRUA',0],
                'IR4' : ['CPRUA',1],
                'IR3' : ['CPRUA',2],
                'IR2' : ['CPRUA',3],
                'UV1' : ['CPRUB',0],
                'UV2' : ['CPRUB',1],
                'NADIR' : ['CPRUB',2]}
    
start_sza = 97 # sza for which nadir measurement starts
stop_sza = 95 # sza for which nadir measurement ends

start_TPlat = 48.5 # TPlat for which UV measurement starts
stop_TPlat = 42.5 # TPlat for which UV measurement starts



# %%

def timeline_stat(df,time_sampling,df_loc):    

    nb_expected_images = np.zeros((7,len(time_sampling)-1))
    nb_generated_images = np.zeros((7,len(time_sampling)-1))
    nb_expected_images_default = np.zeros((7,len(time_sampling)-1)) 

    for channel_ind in range(1,8):
    
        channel = channels[channel_ind]

        for i in range(len(time_sampling)-1):
            nb_expected_images_default[channel_ind-1,i] = (time_sampling[i+1]-time_sampling[i])/sampling_rates[channel]

        # compute the number of expected images in each time intervall
        if channel in ['IR1','IR2','IR3','IR4']:
            for i in range(len(time_sampling)-1):
                nb_expected_images[channel_ind-1,i] = (time_sampling[i+1]-time_sampling[i])/sampling_rates[channel]

        if channel == 'NADIR':
            # list having the start and end times of all the NADIR measurement windows
            nadir_measurements = [] 
            try :
                df_loc = df_loc.sort_values('EXPDate')
                orb_start = df_loc.iloc[0]['EXPDate']
                orb_end = df_loc.iloc[0]['EXPDate']
                for i in range(len(df_loc)-1):
                    if df_loc.iloc[i]['nadir_sza']<start_sza and df_loc.iloc[i+1]['nadir_sza']>start_sza:
                        orb_start = df_loc.iloc[i+1]['EXPDate']             
                    elif df_loc.iloc[i]['nadir_sza']>stop_sza and df_loc.iloc[i+1]['nadir_sza']<stop_sza:
                        orb_end = df_loc.iloc[i]['EXPDate']
                        nadir_measurements.append([orb_start,orb_end])
                if orb_end < orb_start:
                    orb_end = df_loc.iloc[-1]['EXPDate']
                    nadir_measurements.append([orb_start,orb_end])

                # compute number of expected nadir images
                for i in range(len(time_sampling)-1):            
                    nadir_duration = timedelta(seconds=0)
                    start = time_sampling[i]
                    end = time_sampling[i+1]
                    for j in range(len(nadir_measurements)):
                        start_win = nadir_measurements[j][0]
                        end_win = nadir_measurements[j][1]
                        if start_win<end and end_win>=start:
                            nadir_duration += min(end-start_win,end_win-start,end_win-start_win,end-start)
            
                    nb_expected_images[channel_ind-1,i] = nadir_duration/sampling_rates[channel]

            except :
                for i in range(len(time_sampling)-1):            
                    nb_expected_images[channel_ind-1,i] = 0
            
            


        if channel in ['UV1','UV2']:
            # list having the start and end times of all the NADIR measurement windows
            uv_measurements = [] 
            try :
                df_loc = df_loc.sort_values('EXPDate')
                orb_start = df_loc.iloc[0]['EXPDate']
                orb_end = df_loc.iloc[0]['EXPDate']
                for i in range(len(df_loc)-1):
                    if df_loc.iloc[i]['TPlat']<start_TPlat and df_loc.iloc[i+1]['TPlat']>start_TPlat:
                        orb_start = df_loc.iloc[i+1]['EXPDate']             
                    elif df_loc.iloc[i]['TPlat']>stop_TPlat and df_loc.iloc[i+1]['TPlat']<stop_TPlat:
                        orb_end = df_loc.iloc[i]['EXPDate']
                        uv_measurements.append([orb_start,orb_end])
                if orb_end < orb_start:
                    orb_end = df_loc.iloc[-1]['EXPDate']
                    uv_measurements.append([orb_start,orb_end])           
            
                # compute number of expected nadir images
                for i in range(len(time_sampling)-1):            
                    uv_duration = timedelta(seconds=0)
                    start = time_sampling[i]
                    end = time_sampling[i+1]
                    for j in range(len(uv_measurements)):
                        start_win = uv_measurements[j][0]
                        end_win = uv_measurements[j][1]
                        if start_win<end and end_win>=start:
                            uv_duration += min(end-start_win,end_win-start,end_win-start_win,end-start)
                
                    nb_expected_images[channel_ind-1,i] = uv_duration/sampling_rates[channel]
            
            except :
                for i in range(len(time_sampling)-1):            
                    nb_expected_images[channel_ind-1,i] = 0

        for i in range(len(time_sampling)-1):
            start = time_sampling[i]
            end = time_sampling[i+1]
            nb_generated_images[channel_ind-1,i] = len(df[(df['CCDSEL']==channel_ind) & (start<=df['EXPDate']) & (df['EXPDate']<end)])

    return(nb_expected_images_default,nb_expected_images,nb_generated_images)
        


#%%
    
def timeline_plot(data,time_sampling,title,line_labels,file=None,show_plot=False):

    lim_red = 0
    lim_orange = 0.5
    lim_blue = 0.8
    lim_green = 0.95
    lim_purple = 1.05
    width = 0.9


    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    # iteration over channels
    for line_ind in range(len(line_labels)):
        
        line_label = line_labels[line_ind]

        # compute the ratio nb of images/expected number of images

        ax.hlines(y=line_label,xmin=time_sampling[0],xmax=time_sampling[0],color='white') # some invisible line to have a working plot
        for i in range(len(time_sampling)-1):
            start = time_sampling[i]
            end = time_sampling[i+1]
            value = data[line_ind,i]
            color = 'white'
            #print(value)
            if type(value) == type(None):
                color = 'white'
            # if nb_expected_images[line_ind,i] == 0 and nb_generated_images[line_ind,i]>0:
            #     color = 'purple'
            # if nb_expected_images[line_ind,i] > 0:
            #     ratio = nb_generated_images[line_ind,i]/nb_expected_images[line_ind,i]                
            elif value == 0.0:
                color ='black'
            elif (lim_red<value) and (value<=lim_orange):
                color = 'red'
            elif (lim_orange<value) and (value<=lim_blue):
                color = 'orange'
            elif (lim_blue<value) and (value<=lim_green):
                color = 'tab:blue'
            elif (lim_green<value) and (value<=lim_purple):
                color = 'green'
            elif (lim_purple<value):
                color = 'purple'
                #print(ratio)    
            ax.add_patch(Rectangle((start,line_ind-width*0.5),end-start,width,color=color))

    # legend
    legend_elements = [Patch(facecolor='black',label=f"0 == ratio"),    
                    Patch(facecolor='red',label=f"{lim_red*100:.0f}% < ratio <= {lim_orange*100:.0f}%"),
                    Patch(facecolor='orange',label=f"{lim_orange*100:.0f}% < ratio <= {lim_blue*100:.0f}%"),
                    Patch(facecolor='tab:blue',label=f"{lim_blue*100:.0f}% < ratio <= {lim_green*100:.1f}%"),
                    Patch(facecolor='green',label=f"{lim_green*100:.0f}% < ratio <= {lim_purple*100:.1f}%"),
                    Patch(facecolor='purple',label=f"{lim_purple*100:.1f}% < ratio")]
    ax.set_xlabel("Date")
    ax.set_title(title)
    ax.legend(handles=legend_elements,loc='upper left')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)




#%%
def multi_timeline(dataframes,dataframe_labels,sampling_period=timedelta(seconds=120),output_folder=None,show_plot=False):

    df_loc = dataframes[0]

    # define time sampling
    start = min([min(df['TMHeaderTime']) for df in dataframes])
    end = max([max(df['TMHeaderTime']) for df in dataframes])
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)

    total_expected_images = np.zeros((len(dataframes),7,len(time_sampling)-1))
    total_generated_images = np.zeros((len(dataframes),7,len(time_sampling)-1))
    

    for i in range(len(dataframes)):
        df = dataframes[i]
        try :
            sza = df.iloc[0]['nadir_sza']
            df_loc = df
        except :
            print(f"There is no geolocation data in {dataframe_labels[i]}")

    for i in range(len(dataframes)):
        print(f"Checking data generation on {dataframe_labels[i]}")
        df = dataframes[i]
        title = f"nb of images/expected nb of images ({dataframe_labels[i]})"
        nb_expected_images_default,nb_expected_images,nb_generated_images = timeline_stat(df,time_sampling,df_loc)
        total_expected_images[i,:,:] = nb_expected_images
        total_generated_images[i,:,:] = nb_generated_images

        

        line_labels = {0:'IR1',1:'IR4',2:'IR3',3:'IR2',4:'UV1',5:'UV2',6:'NADIR'}
        data = np.zeros_like(nb_expected_images)
        data = np.where(nb_expected_images!=0,nb_generated_images/nb_expected_images,None)
        data = np.where((nb_expected_images==0)&(nb_generated_images>0),nb_generated_images/nb_expected_images_default,data)
        start = min(df['EXPDate'])
        end = max(df['EXPDate'])
        file_path = None
        if type(output_folder) != type(None):
            file_path = f"{output_folder}/image_generation_channels.png"
        timeline_plot(data,time_sampling,title,line_labels=line_labels,file=file_path,show_plot=show_plot)


    
        # for k in range(len(nb_expected_images)):
        #     print(nb_generated_images[k])

    nb_total_expected_images = np.sum(total_expected_images,axis=1)
    nb_total_generated_images = np.sum(total_generated_images,axis=1)
    total_data = np.zeros_like(nb_total_expected_images)
    total_data = np.where(nb_total_expected_images!=0,nb_total_generated_images/nb_total_expected_images,None)
    file_path = None
    if type(output_folder) != type(None):
        file_path = f"{output_folder}/image_generation_sum.png"
    timeline_plot(total_data,time_sampling,"nb of images/expected nb of images (all channels)",line_labels=dataframe_labels,file=file_path,show_plot=show_plot)


    # processing success

    print("Plotting processing success rate")
    nb_im_origin = np.zeros((len(dataframes),len(time_sampling)-1))
    nb_im_processed = np.zeros((len(dataframes),len(time_sampling)-1))
    processing_labels = []
   
    for i in range(len(dataframes)-1):
        processing_labels.append(f"{dataframe_labels[i]} --> {dataframe_labels[i+1]}")
        nb_im_origin[i,:] = np.sum(total_generated_images[i,:,:],axis=(0,1))
        nb_im_processed[i,:] = np.sum(total_generated_images[i+1,:,:],axis=(0,1))
    processing_labels.append(f"{dataframe_labels[0]} --> {dataframe_labels[-1]}")
    nb_im_origin[-1,:] = np.sum(total_generated_images[0,:,:],axis=(0,1))
    nb_im_processed[-1,:] = np.sum(total_generated_images[-1,:,:],axis=(0,1))
    #print(processing_labels)

    processing_data = np.zeros_like(nb_im_origin)
    processing_data = np.where(nb_im_origin!=0,nb_im_processed/nb_im_origin,None)
    #total_data = np.where((nb_im_origin==0)&(nb_im_processed>0),nb_im_processed/nb_expected_images_default,data)
    
    file_path = None
    if type(output_folder) != type(None):
        file_path = f"{output_folder}/image_processing.png"
    timeline_plot(processing_data,time_sampling,"processing success rate (nb processed images/nb of images)",line_labels=processing_labels,file = file_path,show_plot=show_plot)

    
#%%

def temperatureCRBD_plot(dataframe,title='',file=None,show_plot=False):

    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    for channel_ind in range(1,8):  
        channel = channels[channel_ind]
        df = dataframe[dataframe['CCDSEL']==channel_ind]
        temp_ADC = []

        try: 
            for i in range(len(df)):
                temp_ADC.append(df.iloc[i]["temperature_ADC"])
        except: # if it's level0 data            
            ADC_temp_in_mV = df["TEMP"] / 32768 * 2048
            tmp = 1.0 / 0.85 * ADC_temp_in_mV - 296
            for i in range(len(df)):
                temp_ADC.append(tmp.iloc[i])


        if len(temp_ADC)>2:
            # removing the first 2 measurements after the sensor is back on, as they are usually unrealable
            temp_ADC[0] = None
            temp_ADC[1] = None
            for i in range(1,len(temp_ADC)-1):
                if (df.iloc[i]['EXPDate']-df.iloc[i-1]['EXPDate']) > timedelta(seconds = 100):
                    temp_ADC[i] = None
                    temp_ADC[i+1] = None
            
            ax.plot(df['EXPDate'],temp_ADC,label=channel,marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Temperature in each CRB-D {title}')
    ax.set_ylabel('Temperature in C')
    ax.set_xlabel('EXPDate')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)
    

def temperatureHTR_plot(dataframe,title='',file=None,show_plot=False,sampling_period=timedelta(seconds=60)):

    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    temp_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)

    for HTR_name in ['HTR1A','HTR1B','HTR2A','HTR2B','HTR8A','HTR8B']:  
        temp_data = dataframe[HTR_name] 
        HTR_temp = np.zeros(len(temp_sampling)-1)    
        for i in range(len(temp_sampling)-1):
            start = temp_sampling[i]
            end = temp_sampling[i+1]
            HTR_temp[i] = np.nanmean(temp_data[(start<=dataframe['TMHeaderTime']) & (dataframe['TMHeaderTime']<=end)])
        
        ax.plot(temp_sampling[:-1],HTR_temp,label=HTR_name,marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Temperature in each heater sensor (averaged over {sampling_period.total_seconds():.0f} s) {title}')
    ax.set_ylabel('Temperature in C')
    ax.set_xlabel('TMHeaderTime')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)


  

def PWRV_plot(dataframe,title='',file=None,show_plot=False,sampling_period=timedelta(seconds=60)):

    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)

    for voltage in ['PWRP32V','PWRP16V','PWRM16V']:  
        volt_data = dataframe[voltage] 
        VOLT = np.zeros(len(time_sampling)-1)    
        for i in range(len(time_sampling)-1):
            start = time_sampling[i]
            end = time_sampling[i+1]
            VOLT[i] = np.nanmean(volt_data[(start<=dataframe['TMHeaderTime']) & (dataframe['TMHeaderTime']<=end)])
        
        ax.plot(time_sampling[:-1],VOLT,label=voltage,marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Voltage in each bus (averaged over {sampling_period.total_seconds():.0f} s) {title}')
    ax.set_ylabel('Voltage in V')
    ax.set_xlabel('TMHeaderTime')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)


def PWRC_plot(dataframe,title='',file=None,show_plot=False,sampling_period=timedelta(seconds=60)):

    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)
    for current in ['PWRP32C','PWRP16C','PWRM16C','PWRP3C3']:  
        curr_data = dataframe[current] 
        CURR = np.zeros(len(time_sampling)-1)    
        for i in range(len(time_sampling)-1):
            start = time_sampling[i]
            end = time_sampling[i+1]
            CURR[i] = np.nanmean(curr_data[(start<=dataframe['TMHeaderTime']) & (dataframe['TMHeaderTime']<=end)])
        
        ax.plot(time_sampling[:-1],CURR,label=current,marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Current in each bus (averaged over {sampling_period.total_seconds():.0f} s) {title}')
    ax.set_ylabel('Current in A')
    ax.set_xlabel('TMHeaderTime')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)



def PWRT_plot(dataframe,title='',file=None,show_plot=False,sampling_period=timedelta(seconds=60)):

    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)
    temp_data = dataframe['PWRT'] 
    TEMP = np.zeros(len(time_sampling)-1)    
    for i in range(len(time_sampling)-1):
        start = time_sampling[i]
        end = time_sampling[i+1]
        TEMP[i] = np.nanmean(temp_data[(start<=dataframe['TMHeaderTime']) & (dataframe['TMHeaderTime']<=end)])
        
    ax.plot(time_sampling[:-1],TEMP,label='PWRT',marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Power module temperature (averaged over {sampling_period.total_seconds():.0f} s) {title}')
    ax.set_ylabel('Temperature  in C')
    ax.set_xlabel('TMHeaderTime')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)


#%%

def CPRU_overvoltage_plot(dataframe,sampling_period=timedelta(seconds=600),file=None,show_plot=False):
    
    # define time sampling
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)
    
    # width of the line
    width = 0.9
        


    # plotting overvoltage   
    fig, ax = plt.subplots(figsize=(20,10),dpi=250)
    # iteration over channels
    for line_ind in range(7):
        ccdsel = line_ind + 1
        channel = channels[ccdsel]        
        portnum = Port_dic[channel][1]
        cprunum = Port_dic[channel][0]
        df = dataframe[dataframe['SID']==cprunum]
        
        # compute the ratio nb of images/expected number of images

        ax.hlines(y=channel,xmin=time_sampling[0],xmax=time_sampling[0],color='white') # some invisible line to have a working plot
        for i in range(len(time_sampling)-1):
            start = time_sampling[i]
            end = time_sampling[i+1]
            overvoltages = df[(start<=df['TMHeaderTime']) & (df['TMHeaderTime']<end)][f"Overvoltage{portnum}"]
            color = 'white'
            #print(value)
            if len(overvoltages)>0:
                if any(overvoltages):
                    color = 'red'
                else: color='green'
            ax.add_patch(Rectangle((start,line_ind-width*0.5),end-start,width,color=color))

    # legend
    legend_elements = [Patch(facecolor='white',label="no data"),    
                    Patch(facecolor='red',label="Overvoltage"),
                    Patch(facecolor='green',label="No overvoltage")]
    ax.set_xlabel("Date")
    ax.set_title("Overvoltage")
    ax.legend(handles=legend_elements,loc='upper left')
    if type(file) != type(None):
        fig.savefig(file)
    if show_plot:
        plt.show(block=False)


def CPRUV_plot(dataframe,output_folder=None,show_plot=False,sampling_period=timedelta(seconds=600)):    

    for voltage in ['VGATE','VSUBS','VRD','VOD']:
        print(f"Plotting {voltage} values (accross all channels)")
        if type(output_folder) != type(None):
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            file = f"{output_folder}/{voltage}.png"
        # plotting overvoltage   
        fig, ax = plt.subplots(figsize=(20,10),dpi=250)
        # iteration over channels
        for line_ind in range(7):
            ccdsel = line_ind + 1
            channel = channels[ccdsel]        
            portnum = Port_dic[channel][1]
            cprunum = Port_dic[channel][0]
            df = dataframe[dataframe['SID']==cprunum]
            ax.plot(df['TMHeaderTime'],df[f"{voltage}{portnum}"],label=channel,linestyle='',marker='o')        

        # legend
        ax.set_xlabel("TMHeaderTime")
        ax.set_ylabel(f"{voltage} (V)")
        ax.set_title(f"{voltage}")
        ax.legend()
        if type(file) != type(None):
            fig.savefig(file)
        if show_plot:
            plt.show(block=False)


    for channel in ['IR1','IR2','IR3','IR4','UV1','UV2','NADIR']:
        print(f"Plotting {channel} summary")
        if type(output_folder) != type(None):
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            file = f"{output_folder}/{channel}_summary.png"

        ccdsel = CCDSEL[channel]
        portnum = Port_dic[channel][1]
        cprunum = Port_dic[channel][0]
        df = dataframe[dataframe['SID']==cprunum]
        #fig.tight_layout()
        
        fig, (ax_volt,ax_bool) = plt.subplots(2, gridspec_kw={'height_ratios': [4, 1]},figsize=(20,10),dpi=250)
        fig.suptitle(f"{channel} summary")
        # iteration over voltages
        for voltage in ['VGATE','VSUBS','VRD','VOD']:            
            ax_volt.plot(df['TMHeaderTime'],df[f"{voltage}{portnum}"],label=voltage,linestyle='',marker='o')        

        # legend
        ax_volt.set_xlabel("TMHeaderTime")
        ax_volt.set_ylabel(f"Voltage (V)")
        ax_volt.legend()


        # define time sampling
        start = min(df['TMHeaderTime'])
        end = max(df['TMHeaderTime'])
        time_sampling = pd.date_range(start=start,
                    end=end,
                    periods=(end-start).total_seconds()/sampling_period.total_seconds() + 1,tz=timezone.utc)
        
        #print(time_sampling)
        # width of the line
        width = 0.9

        for key in ['Overvoltage','Power']:
            #ax_bool.hlines(y=key,xmin=time_sampling[0],xmax=time_sampling[0],color='white') # some invisible line to have a working plot
            for i in range(len(time_sampling)-1):
                start = time_sampling[i]
                end = time_sampling[i+1]
                data = df[(start<=df['TMHeaderTime']) & (df['TMHeaderTime']<end)][f"{key}{portnum}"]
                #print(data)
                color = 'white'
                #print(value)
                if len(data)>0:
                    if any(data):
                        color = 'green'
                    else: color='red'
                #print(f"--> {color}")
                #ax_bool.add_patch(Rectangle((start,line_ind-width*0.5),end-start,width,color=color))
                ax_bool.hlines(y=key,xmin=start,xmax=end,color=color)

        # legend
        legend_elements = [Patch(facecolor='white',label="no data"),    
                        Patch(facecolor='red',label="False"),
                        Patch(facecolor='green',label="True")]
        ax_bool.set_xlabel("Date")
        ax_bool.legend(handles=legend_elements,loc='upper left')
        ax_bool.set_ylim([-0.5,1.5])


        if type(file) != type(None):
            fig.savefig(file)
        if show_plot:
            plt.show(block=False)



  

# %%


def read_MATS_payload_data(start_date,end_date,data_type='HTR',filter=None,version='0.3'):
    """Reads the payload data between the specified times. 

    Args:
        start (datetime):           Read payload data from this time (inclusive).
        stop (datetime):            Read payload data up to this time (inclusive).
        data_type (str):            key describing the different types of data :
                                    CCD, CPRU, HTR, PWR, STAT, TCV, PM
        filter (Optional[dict]):    Extra filters of the form:
                                    `{fieldname1: [min, max], ...}`
                                    (Default: None)

    Returns:
        DataFrame:      The payload data.
    """

    session = boto3.session.Session(profile_name="mats")
    credentials = session.get_credentials()
    filesystem = f'ops-payload-level0-v{version}'
    file = f"{filesystem}/{data_type}"

    s3 = fs.S3FileSystem(
        secret_key=credentials.secret_key,
        access_key=credentials.access_key,
        region=session.region_name,
        session_token=credentials.token)
    
    if start_date.tzinfo == None:
        start_date = start_date.replace(tzinfo=timezone.utc)
    if end_date.tzinfo == None:
        end_date = end_date.replace(tzinfo=timezone.utc)

    dataset = ds.dataset(
        file,
        filesystem=s3,
        )
    filterlist = (
        (ds.field("TMHeaderTime") >= pd.Timestamp(start_date))
        & (ds.field("TMHeaderTime") <= pd.Timestamp(end_date))
    )
    if filter != None:
        for variable in filter.keys():
            filterlist &= (
                (ds.field(variable) >= filter[variable][0])
                & (ds.field(variable) <= filter[variable][1])
            )

    table = dataset.to_table(filter=filterlist)
    dataframe = table.to_pandas()
    dataframe.reset_index(inplace=True)
    dataframe.set_index('TMHeaderTime',inplace=True)
    dataframe.sort_index(inplace=True)
    dataframe.reset_index(inplace=True)

    return dataframe


def read_MATS_data_custom(start_date,end_date,filter=None,version='0.4',level='1a',columns=None):
    '''
    WARNING !!! No calibration !!!
    '''
    session = boto3.session.Session(profile_name="mats")
    credentials = session.get_credentials()

    s3 = fs.S3FileSystem(
        secret_key=credentials.secret_key,
        access_key=credentials.access_key,
        region=session.region_name,
        connect_timeout=10,
        session_token=credentials.token)

    if start_date.tzinfo == None:
        start_date = start_date.replace(tzinfo=timezone.utc)
    if end_date.tzinfo == None:
        end_date = end_date.replace(tzinfo=timezone.utc)

    if level == '1b' and version == "0.4":
        filesystem = f"ops-payload-level{level}-v{version}" + "/ops-payload-level1a-v0.5"
    elif level == '1b' and version == "0.3":
        filesystem = f"ops-payload-level{level}-v{version}" + "/ops-payload-level1a-v0.4"
    elif level == '0':
        filesystem = f"ops-payload-level{level}-v{version}" + "/CCD"
    else:
        filesystem = f"ops-payload-level{level}-v{version}"
    
    ccd_data = read_ccd_data_in_interval_custom(start_date, end_date, filesystem, s3,filter=filter,columns=columns)

    return (ccd_data)


def read_ccd_data_in_interval_custom(
    start: datetime,
    stop: datetime,
    path: str,
    filesystem: Optional[pa.fs.FileSystem] = None,
    filter: Optional[Dict[str, Sequence[float]]] = None,
    metadata: bool = False,
    columns: Optional[Sequence[str]] = None,
) -> Union[DataFrame, Tuple[DataFrame, pq.FileMetaData]]:
    """Reads the CCD data and metadata from the specified path or S3 bucket
    between the specified times. Optionally read file metadata.

    Args:
        start (datetime):           Read CCD data from this time (inclusive).
        stop (datetime):            Read CCD data up to this time (inclusive).
        path (str):                 Path to dataset. May be a directory or a
                                    bucket, depending on filesystem.
        filesystem (FileSystem):    Optional. File system to read. If not
                                    specified will assume that path points to
                                    an ordinary directory disk. (Default: None)
        filter (Optional[dict]):    Extra filters of the form:
                                    `{fieldname1: [min, max], ...}`
                                    (Default: None)
        metadata (bool):            If True, return Parquet file metadata along
                                    with data frame. (Default: False)

    Returns:
        DataFrame:      The CCD data.
        FileMetaData:   File metadata (optional).
    """

    if start.tzinfo is None:
        start.replace(tzinfo=timezone.utc)
    if stop.tzinfo is None:
        stop.replace(tzinfo=timezone.utc)

    partitioning = ds.partitioning(
        schema=pa.schema(
            [
                ("year", pa.int16()),
                ("month", pa.int8()),
                ("day", pa.int8()),
                ("hour", pa.int8()),
            ]
        ),
    )

    dataset = ds.dataset(
        path,
        filesystem=filesystem,
        partitioning=partitioning,
    )

    start_with_margin =  start - timedelta(hours=1)
    stop_with_margin = stop + timedelta(hours=1)

    partition_filter = (
        ds.field("year") * 1000000
        + ds.field("month") * 10000
        + ds.field("day") * 100
        + ds.field("hour")
        >= start_with_margin.year * 1000000
        + start_with_margin.month * 10000
        + start_with_margin.day * 100
        + start_with_margin.hour
    ) & (
        ds.field("year") * 1000000
        + ds.field("month") * 10000
        + ds.field("day") * 100
        + ds.field("hour")
        <= stop_with_margin.year * 1000000
        + stop_with_margin.month * 10000
        + stop_with_margin.day * 100
        + stop_with_margin.hour
    )


    filterlist = (
        (ds.field("EXPDate") >= Timestamp(start))
        & (ds.field("EXPDate") <= Timestamp(stop))
    )
    if filter != None:
        for variable in filter.keys():
            filterlist &= (
                (ds.field(variable) >= filter[variable][0])
                & (ds.field(variable) <= filter[variable][1])
            )

    if columns != None and 'TMHeaderTime' not in columns: # the column TMHeaderTime has to be always selected as it is set as index
        columns.append('TMHeaderTime')

    table = dataset.to_table(filter=partition_filter & filterlist,columns=columns)

    dataframe = table.to_pandas()
    dataframe.reset_index(inplace=True)
    dataframe.set_index('TMHeaderTime',inplace=True)
    dataframe.sort_index(inplace=True)
    dataframe.reset_index(inplace=True)

    if metadata:
        return dataframe, table.schema.metadata
    return dataframe

# %%


