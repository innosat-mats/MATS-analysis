#%%
import xarray as xr

import datetime as DT
from matplotlib import pyplot as plt


path = 'https://bolin.su.se/data/s3/data/mats-level-1b-limb-cropd-1.0/mats-level-1b-limb-cropd-IR1.zarr'

# Get xarray Dataset using zarr engine
ds = xr.open_zarr(path)
data = ds.sel(time=slice('2023-03-01T00:00:00', '2023-04-01T10:00:00'))


#channel='IR2'
#dsfile='/Users/lindamegner/MATS/testdatadownload/'+channel+'_zarr.nc'
#ds = xr.open_dataset(dsfile)
#%%
# select satlat 67.9 deg N plus minus 1.8 deg (100 km is about 1.8 deg in latitude) 
# select satlon 21.1 deg E plus minus 2.4 deg (200 km at 67.9 deg N is about 2.4 deg in longitude)
satlat=67.9
satlon=21.1
satlatdelta=0.5
satlondelta=2.4
dssel = data.where((data.satlat > satlat - satlatdelta) & (data.satlat < satlat + satlatdelta)
                    & (data.satlon > satlon - satlondelta) & (data.satlon < satlon + satlondelta), drop=True)   


#%%
ds_ascen=dssel #.where(dssel['satlat']<dssel['TPlat'])#ascending node
ds_ascen.satlon.plot(linestyle='', marker='o')
plt.title('Times over Esrange area (67.9N, 21.1E) +- ('+str(satlatdelta)+' deg, '+str(satlondelta)+' deg)')

print(ds_ascen.sel(time=slice('2023-03-07T00:00:00', '2023-03-11T10:00:00')).time.values)

# %%
