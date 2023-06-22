#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data#,read_MATS_payload_data
import os
import pandas as pd
import numpy as np
from datetime import datetime, timedelta, timezone
import argparse
import matplotlib.pyplot as plt
from monitoring_functions import multi_timeline,temperatureCRBD_plot,temperatureHTR_plot,read_MATS_payload_data,PWRT_plot,PWRC_plot,PWRV_plot,CPRU_overvoltage_plot,CPRUV_plot,read_MATS_data_custom,schedule_plot,mean_image



pd.set_option('display.max_rows', 500)

#%%

# defining some parameters

# folders to store monitoring data
def_monitoring_folder = "/home/louis/MATS/MATS-Data/Monitoring"

plt.ioff()




#%%
# parsing arguments 

parser = argparse.ArgumentParser(description='arguments for weekly and daily routine monitoring scripts')

parser.add_argument('--outdir', type=str, default=def_monitoring_folder,
                    help='output directory')
parser.add_argument('--start_time', type=str, default='',
                    help='start of the studied time intervall')
parser.add_argument('--stop_time', type=str, default='',
                    help='end of the studied time intervall')
parser.add_argument('--show_plots', type=str, default='False',
                    help='if matplotlib plots are shown')
parser.add_argument('--sampling_period', type=int, default=120,
                    help='time sampling period in seconds')
parser.add_argument('--mode', type=str, default='daily',
                    help='Monitoring mode : daily, weekly, temp, all, power')


args = parser.parse_args()

start_time = args.start_time
stop_time = args.stop_time
output_folder = args.outdir
show_plot = args.show_plots=='True'
sampling_period = timedelta(seconds=args.sampling_period)
mode = args.mode


if mode == 'daily':
    data_processing = True # plots data generation for each channel (compares it to the expected nb of images) and the processing succes rate between levels
    mean_im = True # plots the mean image for each channel
    histo = False # plots the mean image for each channel aswell as the histogram of pixel values
    CRBD_temp = False # plots the temperature in each CRB-D
    HTR = True # plots the temperature for each Heater and raises some warnings if the temperature isn't nominal
    PWRT = False # plots the temperature in the power module
    PWRC = False # plots the currents from the power module
    PWRV = False # plots the voltage values in the power module
    overvoltage = True # plots a summary of overvoltage events for each channel
    CPRU = False # plots several voltage data coming from the CPRU
    mode_schedule = True # plots the scheduled instrument mode 

elif mode == 'temp':
    data_processing = False
    mean_im = False
    histo = False
    CRBD_temp = True
    HTR = True
    PWRT = True
    PWRC = False
    PWRV = False
    overvoltage = False
    CPRU = False
    mode_schedule = False

elif mode == 'all':
    data_processing = True
    mean_im = True
    histo = True
    CRBD_temp = True
    HTR = True
    PWRT = True
    PWRC = True
    PWRV = True
    overvoltage = True
    CPRU = True
    mode_schedule = True

elif mode == 'power':
    data_processing = False
    mean_im = False
    histo = False
    CRBD_temp = False
    HTR = False
    PWRT = False
    PWRC = True
    PWRV = True
    overvoltage = True
    CPRU = True
    mode_schedule = False    


l0_version = '0.3'
l1a_version = '0.6'
l1b_version = '0.5'


if start_time != '' and stop_time != '':
        start_time = datetime.strptime(start_time,'%Y:%m:%d_%H:%M:%S')
        stop_time = datetime.strptime(stop_time,'%Y:%m:%d_%H:%M:%S')
else: # no time range given : take last day
    dt = datetime.now()
    stop_time = datetime(dt.year, dt.month, dt.day, 0, 0, 0)
    start_time = stop_time - timedelta(days=1)


print('===========================================')
print(f"Monitoring from {start_time} to {stop_time}")

# folders to store figures
if not os.path.exists(output_folder):
        os.mkdir(output_folder)

print(f"Output directory : {output_folder}")

#%%


dataframes = []
dataframe_labels = []

columns_l1a = ['TMHeaderTime', 'EXPDate', 'CCDSEL', 'TEXPMS', 
            'satlat', 'satlon', 'TPlat', 'TPlon', 'nadir_sza','schedule_start_date','schedule_end_date','schedule_id','schedule_name']

columns_l1b = columns_l1a

columns_l0 = ['TMHeaderTime', 'CCDSEL',
            'EXPDate', 'TEXPMS', 'TEMP']

temp_range = {'HTR1A':[10,25],
              'HTR1B':[10,25],
              'HTR2A':[5,20],
              'HTR2B':[5,20],
              'HTR8A':[-25,-5],
              'HTR8B':[-25,-5]}

max_variability = 5 # in K


#%%
if data_processing:
    try :
        print("Importing level 1b data")
        df1b = read_MATS_data_custom(start_time, stop_time,level='1b',version=l1b_version,columns=columns_l1b)
        # df1b = df1b.drop('ImageCalibrated', axis=1)
        if len(df1b)>0:
            dataframes.append(df1b)
            dataframe_labels.append('l1b v0.4')
        else :
            print('No level 1b data')
    except :
        print('No level 1b data')

    try :
        print("Importing level 1a data")
        df1a = read_MATS_data_custom(start_time, stop_time,level='1a',version=l1a_version,columns=columns_l1a)
        # df1a = df1a.drop(columns=['IMAGE','ImageData','id'], axis=1)
        if len(df1a)>0:
            dataframes.append(df1a)
            dataframe_labels.append('l1a v0.5')
        else :
            print('No level 1a data')
    except :
        print('No level 1a data')

    try :
        print("Importing level 0 data")
        df0 = read_MATS_data_custom(start_time, stop_time,level='0',version=l0_version,columns=columns_l0)
        # df0 = df0.drop('ImageData', axis=1)
        if len(df0)>0:
            dataframes.append(df0)
            dataframe_labels.append('l0 v0.3')
        else :
            print('No level 0 data')
    except :
        print('No level 0 data')    
        
    if len(dataframes)>0:
        multi_timeline(dataframes,dataframe_labels,output_folder=output_folder,show_plot=show_plot,sampling_period=sampling_period)
        if not show_plot:
            plt.close('all')


if mean_im:
    try:
        print("Importing level 0 data")
        df0 = read_MATS_data(start_time, stop_time,level='0/CCD',version='0.3')
        print(f"Plotting the mean image for each channel")
        file_path = f"{output_folder}/mean_im.png"
        mean_image(df0,file=file_path,show_plot=show_plot)
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot mean images from l0 v0.3")


if CRBD_temp:
    try:
        print("Importing level 0 data")
        df0 = read_MATS_data_custom(start_time, stop_time,level='0',version=l0_version,columns=columns_l0)
        print(f"Plotting CRB-D temperatures")
        file_path = f"{output_folder}/CRBD_temp.png"
        temperatureCRBD_plot(df0,title='',file=file_path,show_plot=show_plot)
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot CRB-D temperatures from l0 v0.3")


if HTR:
    try:
        print(f"Importing HTR temperature data")
        HTR_df = read_MATS_payload_data(start_time,stop_time,data_type='HTR')
        file_path = f"{output_folder}/HTR_temp.png"
        print(f"Plotting HTR temperatures")
        temperatureHTR_plot(HTR_df,file=file_path,show_plot=show_plot,sampling_period=sampling_period)
        for heater in ['HTR1A','HTR1B','HTR2A','HTR2B','HTR8A','HTR8B']:
            nominal_range = temp_range[heater]
            if max(HTR_df[heater]) > nominal_range[1] :
                print(f"{heater} temperature anormally high")
            if min(HTR_df[heater]) < nominal_range[0] :
                print(f"{heater} temperature anormally low")
            if max(HTR_df[heater]) - min(HTR_df[heater]) > max_variability:
                print(f"High variability on {heater} temperature between {start_time} and {stop_time}")
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot HTR temperatures")

if PWRV or PWRC or PWRT:
    try:
        print(f"Importing PWR data")
        PWR_df = read_MATS_payload_data(start_time,stop_time,data_type='PWR')
        if PWRV:
            print(f"Plotting PWR voltages")
            PWRV_plot(PWR_df,file=f"{output_folder}/PWR_voltage.png",show_plot=show_plot,sampling_period=sampling_period)
        if PWRT:
            print(f"Plotting PWR temperature")
            PWRT_plot(PWR_df,file=f"{output_folder}/PWR_temp.png",show_plot=show_plot,sampling_period=sampling_period)
        if PWRC:
            print(f"Plotting PWR currents")
            PWRC_plot(PWR_df,file=f"{output_folder}/PWR_current.png",show_plot=show_plot,sampling_period=sampling_period)
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot PWR temperatures")

if CPRU or overvoltage:
    try:
        print(f"Importing CPRU data")
        CPRU_df = read_MATS_payload_data(start_time,stop_time,data_type='CPRU')    
        if CPRU:
            print(f"Plotting CPRU data")
            CPRUV_plot(CPRU_df,output_folder=output_folder,show_plot=show_plot,sampling_period=sampling_period)
        if overvoltage:
            print(f"Plotting Overvoltage summary")
            CPRU_overvoltage_plot(CPRU_df,file=f"{output_folder}/CPRU_overvoltage.png",show_plot=show_plot,sampling_period=sampling_period)
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot CPRU data")


if mode_schedule: 
    try:
        print(f"Importing df1a v0.6 data")
        df1a = read_MATS_data(start_time,stop_time,version=l1a_version,level='1a')
        df1a = df1a.drop(columns=['IMAGE','ImageData','id'], axis=1)
        print(f"Plotting schedule data")
        schedule_plot(df1a,file=f"{output_folder}/schedule.png",show_plot=show_plot,column='schedule_name')
        if not show_plot:
            plt.close('all')
    except:
        print(f"Unable to plot schedule data")
        

       

if show_plot:
    plt.show(block=True)
else:
    plt.close('all')
     


# %%
