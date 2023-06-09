#%%
from mats_utils.rawdata.timeline_tools import plot_schedule,load_schedule
import pandas as pd
import datetime as DT

#%%
df = load_schedule('/home/olemar/Projects/Universitetet/MATS/MATS-utility-functions/data/20221221_timeline_schedule.csv')
# %%
plot_schedule(df,start_date = DT.datetime(2023,3,1),end_date = DT.datetime(2023,4,1))
# %%
