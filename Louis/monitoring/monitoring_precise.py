#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT
from mats_utils.plotting.plotCCD import *
from mats_utils.statistiscs.images_functions import create_imagecube
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, timezone
import warnings

import boto3
import re
import pyarrow.parquet as pq  # type: ignore
import pyarrow.dataset as ds
from pyarrow import fs
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

pd.set_option('display.max_rows', 500)

#%%

# defining some parameters
channels = {1:'IR1',2:'IR4',3:'IR3',4:'IR2',5:'UV1',6:'UV2',7:'NADIR'}
sampling_rates={'IR1':timedelta(seconds=6),
                    'IR2':timedelta(seconds=6),
                    'IR3':timedelta(seconds=6),
                    'IR4':timedelta(seconds=6),
                    'UV1':timedelta(seconds=6),
                    'UV2':timedelta(seconds=6),
                    'NADIR':timedelta(seconds=2)}
    
start_sza = 97 # sza for which nadir measurement starts
stop_sza = 95 # sza for which nadir measurement ends

start_TPlat = 48.5 # TPlat for which UV measurement starts
stop_TPlat = 42.5 # TPlat for which UV measurement starts

# folders to store monitoring data
monitoring_folder = "/home/louis/MATS/MATS-Data/Monitoring"


#%%

start_time = datetime(2023, 5, 1, 0, 0, 0)
stop_time = datetime(2023, 5, 1, 0, 10, 0)


# folders to store figures
data_folder = f"{monitoring_folder}{start_time.strftime('%Y:%m:%d')}_{stop_time.strftime('%Y:%m:%d')}"
if not os.path.exists(data_folder):
        os.mkdir(data_folder)


sampling = 'custom'

custom_period = timedelta(minutes=2)

if sampling == 'day':
    start=start_time.replace(hour=0, minute=0, second=0, microsecond=0)
    end=stop_time.replace(hour=0, minute=0, second=0, microsecond=0)
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).days + 1)

elif sampling == 'orbit':
    orbit_filter = {'CCDSEL': [1,7],'satlat':[-1.0,+1.0]}
    df = read_MATS_data(start_time, stop_time,filter=orbit_filter,level='1a',version='0.5')
    df = df[~np.isnan(df['satlat'])].sort_values('EXPDate')
    time_sampling = [start_time.replace(tzinfo=timezone.utc)]
    for i in range(len(df)-1):
        if df.iloc[i]['satlat']<0.0 and df.iloc[i+1]['satlat']>0.0:
            time_sampling.append(df.iloc[i]['EXPDate'])  
    time_sampling.append(stop_time.replace(tzinfo=timezone.utc))    

if sampling == 'custom':
    start=start_time
    end=stop_time
    time_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/custom_period.total_seconds() + 1,tz=timezone.utc)



#%%   



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
        

    
    

# nb_expected_images,nb_generated_images = timeline_stat(df0,time_sampling,df0)


#%%
    
def timeline_plot(data,time_sampling,title,line_labels,file=None):

    lim_red = 0
    lim_orange = 0.5
    lim_blue = 0.8
    lim_green = 0.95
    lim_purple = 1.05
    width = 0.9

    
    fig, ax = plt.subplots(figsize=(20,10))
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
    plt.show(block=False)




# line_labels = {0:'IR1',1:'IR4',2:'IR3',3:'IR2',4:'UV1',5:'UV2',6:'NADIR'}

# timeline_plot(nb_expected_images,nb_generated_images,time_sampling,'test',line_labels=line_labels)
    

#%%
def multi_timeline(dataframes,dataframe_labels,time_sampling):

    df_loc = dataframes[0]

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
        file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_nb_images.png"
        timeline_plot(data,time_sampling,title,line_labels=line_labels,file=file_path)


    
        # for k in range(len(nb_expected_images)):
        #     print(nb_generated_images[k])

    nb_total_expected_images = np.sum(total_expected_images,axis=1)
    nb_total_generated_images = np.sum(total_generated_images,axis=1)
    total_data = np.zeros_like(nb_total_expected_images)
    total_data = np.where(nb_total_expected_images!=0,nb_total_generated_images/nb_total_expected_images,None)
    
    file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_nb_image_sum.png"
    timeline_plot(total_data,time_sampling,"nb of images/expected nb of images (all channels)",line_labels=dataframe_labels,file=file_path)


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

    file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_processing.png"
    timeline_plot(processing_data,time_sampling,"processing success rate (nb processed images/nb of images)",line_labels=processing_labels,file = file_path)

    
#%%

def temperatureCRBD_plot(dataframe,title='',file=None):

    fig, ax = plt.subplots(figsize=(20,10))
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
    ax.set_title(f'Temperature in each CRB-D ({title})')
    ax.set_ylabel('Temperature in C')
    ax.set_xlabel('EXPDate')
    if type(file) != type(None):
        fig.savefig(file)
    plt.show(block=False)

def temperatureHTR_plot(dataframe,title='',file=None):

    fig, ax = plt.subplots(figsize=(20,10))
    T_period = timedelta(seconds=60)
    start = min(dataframe['TMHeaderTime'])
    end = max(dataframe['TMHeaderTime'])
    temp_sampling = pd.date_range(start=start,
                  end=end,
                  periods=(end-start).total_seconds()/T_period.total_seconds() + 1,tz=timezone.utc)

    for HTR_name in ['HTR1A','HTR1B','HTR2A','HTR2B','HTR8A','HTR8B']:  
        temp_data = dataframe[HTR_name] 
        HTR_temp = np.zeros(len(temp_sampling)-1)    
        for i in range(len(temp_sampling)-1):
            start = temp_sampling[i]
            end = temp_sampling[i+1]
            HTR_temp[i] = np.nanmean(temp_data[(start<=dataframe['TMHeaderTime']) & (dataframe['TMHeaderTime']<=end)])
        
        ax.plot(temp_sampling[:-1],HTR_temp,label=HTR_name,marker='o',linestyle='')
    
    ax.legend()
    ax.set_title(f'Temperature in each heater sensor (averaged over {T_period.total_seconds():.0f} s) ({title})')
    ax.set_ylabel('Temperature in C')
    ax.set_xlabel('EXPDate')
    if type(file) != type(None):
        fig.savefig(file)
    plt.show(block=False)


def temp_plot(dataframe,dataframe_label):
    try:
        print(f"Plotting CRB-D temperatures for {dataframe_label}")
        start = min(dataframe['EXPDate'])
        end = max(dataframe['EXPDate'])
        file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_CRBD_temp.png"
        temperatureCRBD_plot(dataframe,title=dataframe_label,file=file_path)
    except:
        print(f"Unable to plot CRB-D temperatures for {dataframe_label}")

    try:
        print(f"Plotting HTR temperatures for {dataframe_label}")
        start = min(dataframe['EXPDate'])
        end = max(dataframe['EXPDate'])
        file_path = f"{data_folder}/{start.strftime('%Y:%m:%d')}_{end.strftime('%Y:%m:%d')}_HTR_temp.png"
        temperatureHTR_plot(dataframe,title=dataframe_label,file=file_path)
    except:
        print(f"Unable to plot HTR temperatures for {dataframe_label}")


    
    

    

# multi_timeline([df0,df1a,df1b],['level 0 v0.3','level 1a v0.5','level 1b v0.4'],time_sampling)

#%%




def multi_level_data(start_time,stop_time,sampling='custom',custom_period = timedelta(minutes=2)):

    # defining the timesampling
    if sampling == 'day':
        start=start_time.replace(hour=0, minute=0, second=0, microsecond=0)
        end=stop_time.replace(hour=0, minute=0, second=0, microsecond=0)
        time_sampling = pd.date_range(start=start,
                    end=end,
                    periods=(end-start).days + 1)

    elif sampling == 'orbit':
        orbit_filter = {'CCDSEL': [1,7],'satlat':[-1.0,+1.0]}
        df = read_MATS_data(start_time, stop_time,pltfilter=orbit_filter,level='1a',version='0.5')
        df = df[~np.isnan(df['satlat'])].sort_values('EXPDate')
        time_sampling = [start_time.replace(tzinfo=timezone.utc)]
        for i in range(len(df)-1):
            if df.iloc[i]['satlat']<0.0 and df.iloc[i+1]['satlat']>0.0:
                time_sampling.append(df.iloc[i]['EXPDate'])  
        time_sampling.append(stop_time.replace(tzinfo=timezone.utc))    

    if sampling == 'custom':
        start=start_time
        end=stop_time
        time_sampling = pd.date_range(start=start,
                    end=end,
                    periods=(end-start).total_seconds()/custom_period.total_seconds() + 1,tz=timezone.utc)
        

    
    dataframes = []
    dataframe_labels = []
    try :
        print("Importing level 0 data")
        df0 = read_MATS_data(start_time, stop_time,level='0',version='0.3')
        dataframes.append(df0)
        dataframe_labels.append('l0 v0.3')
    except :
        print('No level 0 data')

    try :
        print("Importing level 1a data")
        df1a = read_MATS_data(start_time, stop_time,level='1a',version='0.5')
        dataframes.append(df1a)
        dataframe_labels.append('l1a v0.5')
    except :
        print('No level 1a data')

    try :
        print("Importing level 1b data")
        df1b = read_MATS_data(start_time, stop_time,level='1b',version='0.4')
        dataframes.append(df1b)
        dataframe_labels.append('l1b v0.4')
    except :
        print('No level 1b data')
    
    
    if len(dataframes)>0:
        multi_timeline(dataframes,dataframe_labels,time_sampling)  



    return dataframes,dataframe_labels,time_sampling





#%%



dataframes,dataframe_labels,time_sampling = multi_level_data(start_time,stop_time)


#%%


print('Importing L1b data')
df1b = read_MATS_data(start_time, stop_time,level='1b',version='0.4')
print('Plotting temperature data')
temp_plot(df1b,'l1b v0.4')


        
# multi_timeline(dataframes,dataframe_labels,time_sampling)





# %%

# check no geolocation ?
# check temp



plt.show(block=True)




# %%

def HTR_temp(start_time,stop_time,file=None):
    session = boto3.session.Session(profile_name="mats")
    credentials = session.get_credentials()

    s3 = fs.S3FileSystem(
        secret_key=credentials.secret_key,
        access_key=credentials.access_key,
        region=session.region_name,
        session_token=credentials.token)


    dataset = ds.dataset(
        "ops-payload-level0-v0.3/HTR",
        filesystem=s3,
        )

    df = dataset.to_table().to_pandas()
    df = df[(df['TMHeaderTime']>start_time.replace(tzinfo=timezone.utc))&(df['TMHeaderTime']<stop_time.replace(tzinfo=timezone.utc))]
    temperatureHTR_plot(df,file=file)


temp_start = datetime(2023, 2, 19, 0, 0, 0)
temp_stop = datetime(2023, 2, 26, 0, 0, 0)
data_folder = f"{monitoring_folder}{temp_start.strftime('%Y:%m:%d')}_{temp_stop.strftime('%Y:%m:%d')}"
if not os.path.exists(data_folder):
        os.mkdir(data_folder)
HTR_temp(temp_start,temp_stop,file=f"{data_folder}/test.png")

# %%
