#%%
import xarray as xr

import datetime as DT
from matplotlib import pyplot as plt
#%%
path = 'https://bolin.su.se/data/s3/data/mats-level-1b-limb-cropd-1.0/mats-level-1b-limb-cropd-IR1.zarr'

# Get xarray Dataset using zarr engine
ds = xr.open_zarr(path)
data = ds.sel(time=slice('2023-03-01T00:00:00', '2023-05-30T10:00:00'))


#channel='IR2'
#dsfile='/Users/lindamegner/MATS/testdatadownload/'+channel+'_zarr.nc'
#ds = xr.open_dataset(dsfile)
#%%
# select TPlat 67.9 deg N plus minus 1.8 deg (100 km is about 1.8 deg in latitude) 
# select TPlon 21.1 deg E plus minus 2.4 deg (200 km at 67.9 deg N is about 2.4 deg in longitude)
TPlat=-53.8
TPlon=-67.8
TPlatdelta=0.5
TPlondelta=1.
dssel = data.where((data.TPlat > TPlat - TPlatdelta) & (data.TPlat < TPlat + TPlatdelta)
                    & (data.TPlon > TPlon - TPlondelta) & (data.TPlon < TPlon + TPlondelta), drop=True)


#%%

dssel.TPlon.plot(linestyle='', marker='o')
plt.title('Times over Rio Grande area ('+str(TPlat)+'N, '+str(TPlon)+'E) +- ('+str(TPlatdelta)+' deg, '+str(TPlondelta)+' deg)')

#%%
dssel.TPlat.plot(linestyle='', marker='o')
plt.title('Times over Rio Grande area ('+str(TPlat)+'N, '+str(TPlon)+'E) +- ('+str(TPlatdelta)+' deg, '+str(TPlondelta)+' deg)')    



# %%
#plot the selected images

for i in range(len(dssel.time)):
    fig, ax = plt.subplots()
    img=ax.imshow(dssel.ImageCalibrated.isel(time=i), cmap='viridis', aspect=1/4)
    #print TPlat, TPlon, solar zenith angle TPsza  and time in title
    ax.set_title('TPlat: '+str(dssel.TPlat.values[i])+' deg N, TPlon: '+str(dssel.TPlon.values[i])+' deg E, time: '+str(dssel.time.values[i]))
    ax.invert_yaxis()  # Invert y-axis to match image orientation
    plt.show()
   
# %%
