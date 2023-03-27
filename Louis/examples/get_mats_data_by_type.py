from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
from matplotlib import pyplot as plt
import datetime as DT

# load schedule
schedule = pd.read_csv('/home/louis/MATS/timeline_schedule.csv')

start_time = DT.datetime(2022, 11, 22, 8, 15, 0)
stop_time = DT.datetime(2022, 11, 22, 9, 25, 0)

# filter schedule
star_schedule = schedule[schedule.name=='STAR'] 

star_schedule = star_schedule.reset_index(drop=True) #reset indexing of dataframe
star_schedule.start_date = pd.to_datetime(star_schedule.start_date)
star_schedule.end_date = pd.to_datetime(star_schedule.end_date)

start_time2 = star_schedule.start_date[0]
print(star_schedule.start_date[0])
print(start_time)

# add all star measurements
df = read_MATS_data(start_time,stop_time,level='1a', version=0.5, filter=None)
for i in range(1,len(star_schedule)):
    df = df.append(read_MATS_data(star_schedule.start_date[i],star_schedule.end_date[i]))


plt.imshow(df.IMAGE[5],clim=[300,370],origin='lower')

