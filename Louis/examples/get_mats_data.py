from mats_utils.rawdata.read_data import read_MATS_data
import pandas as pd
import datetime as DT

# times for start and stop
start_time = DT.datetime(2023, 2, 12, 6, 0, 0)
stop_time = DT.datetime(2023, 2, 12, 6, 2, 0)

# filter
# filter={'CCDSEL': [7,7]}

df = read_MATS_data(start_time,stop_time)

df

df.columns

df.IMAGE[0]

from matplotlib import pyplot as plt
plt.pcolor(df.IMAGE[0])

image = df.IMAGE[df.CCDSEL==1].iloc[0]
plt.pcolor(image)
plt.show()

