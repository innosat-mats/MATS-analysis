#%%
from mats_utils.rawdata.timeline_tools import plot_schedule,load_schedule
import pandas as pd
import datetime as DT
import os

#%%
directory = '/Users/lindamegner/MATS/MATS-retrieval/data/TimelineSchedules/'


# Get a list of all files in the directory
file_list = os.listdir(directory)

dflong=pd.DataFrame()
# Loop through each file and read it
for file_name in file_list:
    file_path = os.path.join(directory, file_name)
    if os.path.isfile(file_path):

        df = load_schedule(file_path)
        print(f"Reading file: {file_path}")
        dflong = pd.concat([dflong,df],ignore_index=True)


#df = load_schedule('/Users/lindamegner/MATS/MATS-retrieval/data/Schedules/20230309_timeline_schedule.csv')

# %%
#plot_schedule(dflong)
plot_schedule(dflong,start_date = DT.datetime(2023,2,1),end_date = DT.datetime(2023,3,1))



# %%
