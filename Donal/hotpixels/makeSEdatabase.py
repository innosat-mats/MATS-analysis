#%%
import numpy as np
import sqlite3 as sqlite
from glob import glob
import os
import pickle
import datetime as DT
#%%
homecat = os.environ.get('HOME')
files= glob(homecat + '/Downloads/OneDrive_2_10-10-2023/*')
#%%
os.remove(homecat + '/Downloads/hpms/SE.db')
db = sqlite.connect(homecat + '/Downloads/hpms/SE.db')
cur = db.cursor()
cur.execute("create table if not exists SingleEvents ('datetime' NUMERIC ,'channel' TEXT,'BildNumber' INTEGER,'X' INTEGER, 'Y' INTEGER,'Value' INTEGER,'Sigma' INTEGER)")#,PRIMARY KEY (datetime,channel))")
insertstr = "insert or replace into SingleEvents values (?, ?, ?, ?, ?, ?, ?)"

# %%
nfiles=len(files)
for ifile,f in enumerate(files):
    filename= f.split("/")[-1]
    ch=filename[3:6]
    print("{:s} {:d}/{:d} ".format(filename,ifile,nfiles),end="\r")
    data=np.loadtxt(f)  
    
    for i in range(data.shape[0]):   
        #print(i,filename,data[i,:])  
        date=DT.datetime.strptime( filename[7:-4]+f"{data[i,1]:06.0f}", "%Y%m%d%H%M%S")  
        #print(i,filename,ch,date)  
        cur.execute(insertstr,(date,ch,data[i,0],data[i,2],data[i,3],data[i,4],data[i,5]))
db.commit()
                     
# %%
