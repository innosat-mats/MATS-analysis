#%% Import modules
#%matplotlib qt5
from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import os
import numpy as np
from datetime import datetime, timedelta
import argparse
from monitoring_functions import timeline_stat,timeline_plot,multi_timeline,temperatureCRBD_plot,temperatureHTR_plot,read_MATS_payload_data,PWRT_plot,PWRC_plot,PWRV_plot

pd.set_option('display.max_rows', 500)

#%%

# defining some parameters

# folders to store monitoring data
def_monitoring_folder = "/home/louis/MATS/MATS-Data/Monitoring"

sampling = 'custom'

custom_period = timedelta(minutes=2)


#%%
# parsing arguments

parser = argparse.ArgumentParser(description='arguments for weekly and daily routine monitoring scripts')

parser.add_argument('--outdir', type=str, default=def_monitoring_folder,
                    help='output directory')
parser.add_argument('--start_time', type=str, default='',
                    help='start of the studied time intervall')
parser.add_argument('--stop_time', type=str, default='',
                    help='end of the studied time intervall')
parser.add_argument('--show_plots', type=bool, default=False,
                    help='if matplotlib plots are shown')

args = parser.parse_args()

start_time = args.start_time
stop_time = args.stop_time
monitoring_folder = args.outdir
show_plot = args.show_plots

if start_time != '' and stop_time != '':
        week_start = datetime.strptime(start_time,'%Y:%m:%d_%H:%M:%S')
        week_end = datetime.strptime(stop_time,'%Y:%m:%d_%H:%M:%S')
else: # no time range given : take last week
    dt = datetime.now() 
    day = datetime(dt.year, dt.month, dt.day, 0, 0, 0)
    week_start = day - timedelta(days=day.weekday()+7)
    week_end = week_start + timedelta(days=7)





#%%


print('===========================================')
print(f"Monitoring from {week_start} to {week_end}")


data_folder = f"{monitoring_folder}/weekly_monitoring_{week_start.strftime('%Y_%m_%d')}_{week_end.strftime('%Y_%m_%d')}"
if not os.path.exists(data_folder):
        os.mkdir(data_folder)

print(f"Output directory : {data_folder}")

#%%

try:
    print(f"Importing HTR temperature data")
    HTR_df = read_MATS_payload_data(week_start,week_end,data_type='HTR')
    file_path = f"{data_folder}/HTR_temp.png"
    print(f"Plotting HTR temperatures")
    temperatureHTR_plot(HTR_df,file=file_path,show_plot=show_plot)
except:
    print(f"Unable to plot HTR temperatures")


try:
    print(f"Importing PWR data")
    PWR_df = read_MATS_payload_data(week_start,week_end,data_type='PWR')
    print(f"Plotting PWR voltages")
    PWRV_plot(PWR_df,file=f"{data_folder}/PWR_voltage.png",show_plot=show_plot)
    print(f"Plotting PWR temperature")
    PWRT_plot(PWR_df,file=f"{data_folder}/PWR_temp.png",show_plot=show_plot)
    print(f"Plotting PWR currents")
    PWRC_plot(PWR_df,file=f"{data_folder}/PWR_current.png",show_plot=show_plot)
except:
    print(f"Unable to plot PWR temperatures")
# %%
