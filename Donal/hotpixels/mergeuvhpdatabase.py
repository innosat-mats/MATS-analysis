#%%
import numpy as np
import sqlite3 as sqlite
from glob import glob

import os
import pickle
import datetime as DT
#%%
homecat = os.environ.get('HOME')
#%%
dbmain = sqlite.connect(homecat + '/Downloads/hpms/hpms.db')
curmain = dbmain.cursor()
insertstr = "insert or replace into hotpixelmaps values (?, ?, ?, ? )"
dbuv1 = sqlite.connect(homecat + '/Downloads/hpms/ir3.db')
curuv1 = dbuv1.cursor()
#%%
# Delete rows in dbmain with channel "UV1" or "UV2"
curmain.execute("DELETE FROM hotpixelmaps WHERE channel IN ('IR3', 'IR4')")
dbmain.commit()

# Fetch rows from dbuv1
curuv1.execute("SELECT * FROM hotpixelmaps WHERE channel IN ('IR3', 'IR4')")
rows = curuv1.fetchall()

# Insert rows into dbmain with formatted datetime
for row in rows:
    formatted_row = list(row)
    formatted_row[1] = DT.datetime.strptime(
        formatted_row[1], '%Y-%m-%d %H:%M:%S.%f%z'
    ).strftime('%Y-%m-%d %H:%M:%S')
    curmain.execute(insertstr, formatted_row)

dbmain.commit()
dbuv1.close()
dbmain.close()
# %%
# Reopen the database connection
dbmain = sqlite.connect(homecat + '/Downloads/hpms/hpms.db')
curmain = dbmain.cursor()

# Fetch the earliest date and corresponding HPM for each channel
curmain.execute("SELECT channel, MIN(datetime) FROM hotpixelmaps GROUP BY channel")
earliest_dates = curmain.fetchall()
print (earliest_dates)
#%%
# Insert a new row for each channel with the date '2023-02-08 00:00' and the same HPM as the earliest date
date_str = '2023-02-08 00:00'
date = DT.datetime.strptime(date_str, '%Y-%m-%d %H:%M')
datekey = date.year * 100000000 + date.month * 1000000 + date.day * 10000 + date.hour * 100 + date.minute
for channel, earliest_date in earliest_dates:
    curmain.execute("SELECT hpm FROM hotpixelmaps WHERE channel = ? AND datetime = ?", (channel, earliest_date))
    hpm = curmain.fetchone()[0]
    new_row = (datekey,date, channel,  hpm)
    curmain.execute(insertstr, new_row)

dbmain.commit()
dbmain.close()

# %%
