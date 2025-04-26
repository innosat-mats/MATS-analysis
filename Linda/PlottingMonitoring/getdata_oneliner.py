#%%
from mats_utils.rawdata.read_data import read_MATS_data
import datetime as DT

#filter={'CCDSEL': [5,6]}
df1= read_MATS_data(DT.datetime(2023, 2, 9, 4, 0, 0),  DT.datetime(2023, 2, 9, 20, 0, 0),level='1b',version='0.9')
print('df1',df1.schedule_name.unique())
df2= read_MATS_data(DT.datetime(2023, 2, 9, 17, 0, 0),  DT.datetime(2023, 2, 9, 20, 0, 0),level='1b',version='0.9')
print('df2',df2.schedule_name.unique())


# %%
# find first instance with schedule_name = 'CROPD'
df[df.schedule_name == 'CROPD'].index[0]
df[df.schedule_name == 'CROPD'].TMHeaderTime[0]
