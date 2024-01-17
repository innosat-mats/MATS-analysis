#%%
import numpy as np
import sqlite3 as sqlite
from glob import glob
import hpfunctions
import datetime as DT
import matplotlib.pylab as plt
from mats_utils.rawdata.read_data import read_MATS_data
from mats_utils.geolocation.coordinates import satpos
import pickle
#%%
startdate=DT.datetime(2023,4,22,10,0)
enddate=DT.datetime(2023,4,22,11,10)
startdate=DT.datetime(2023,3,31,21,0)
stopdate=DT.datetime(2023,3,31,22,35)
startdate=DT.datetime(2023,5,5,1,0)
stopdate=DT.datetime(2023,5,5,2,0)
# dftop=read_MATS_data(starttime,stoptime,level="1b",version="0.4")
df=read_MATS_data(startdate,end_date=stopdate,version=0.6)
#%%
print('Number of images = ',df.shape[0])
clim=999
plt.close('all')
ccdnames=('IR1','IR4','IR3','IR2','UV1','UV2','NADIR')
flip=(True,False,True,False,True,True,False)
ir1=df[df.CCDSEL==1].reset_index(drop=True)
ir2=df[df.CCDSEL==4].reset_index(drop=True)
ir3=df[(df.CCDSEL==3)].reset_index(drop=True)
ir4=df[(df.CCDSEL==2)].reset_index(drop=True)
uv1=df[(df.CCDSEL==5)].reset_index(drop=True)
uv2=df[(df.CCDSEL==6)].reset_index(drop=True)
for ch in [ir1,ir2,ir3,ir4]:
    print (ccdnames[ch.CCDSEL.iloc[0]-1],  ch.shape[0])
uv2.shape
# %%
i=399

ch=ir3;clims=[0,60]
ch=ir2;clims=[0,2000]
#ch=ir4;clims=[0,60]
#ch=ir3;clims=[0,10]
#ch=uv2;clims=[0,800]

fig,axis=plt.subplots(1,1,figsize=[8,2])
#image=np.stack(ch.ImageCalibrated.iloc[i])
image=ch.IMAGE.iloc[i]
sp=plt.imshow(image, cmap="viridis", origin="lower", interpolation="none")
axis.axis("auto")
TPlat,TPlon,satalt=satpos(ch.iloc[i])
plt.title("{:4s} Lat = {:8.3f} Lon = {:8.3f} {:s}".format(ccdnames[ch.iloc[i].CCDSEL - 1 ], TPlat, TPlon ,ch.EXPDate.iloc[i].isoformat()))

plt.clim(clims)
plt.colorbar()
plt.show()
# %%
mapdate,HP=hpfunctions.gethpm(startdate,'IR2',)
print(mapdate)
# %%
fig,axis=plt.subplots(1,1,figsize=[8,2])
sp=plt.imshow(image-HP, cmap="viridis", origin="lower", interpolation="none")
axis.axis("auto")
plt.title("{:4s} Lat = {:8.3f} Lon = {:8.3f} {:s}".format(ccdnames[ch.iloc[i].CCDSEL - 1 ], TPlat, TPlon ,ch.EXPDate.iloc[i].isoformat()))
plt.clim(clims)
plt.colorbar()
plt.show()
# %%
fig,axis=plt.subplots(1,1,figsize=[8,2])
sp=plt.imshow(HP, cmap="magma", origin="lower", interpolation="none")
axis.axis("auto")
plt.title("{:4s} Lat = {:8.3f} Lon = {:8.3f} {:s}".format('HotPixels', TPlat, TPlon ,ch.EXPDate.iloc[i].isoformat()))

#plt.clim(clims)
plt.colorbar()
plt.show()
# %%
#ir2f=ir2[ir2.NROW==187].reset_index(drop=True)
for i in range(len(ir1)):
    image=ir1.IMAGE.iloc[i]
    ir1.IMAGE.iloc[i]=image - HP
# %%
with open('ir1mars31.pickle', 'wb') as handle:
     pickle.dump(ir1, handle)
# %%
