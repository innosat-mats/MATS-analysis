
# %%
import xarray as xr
from database_generation.experimental_utils import plot_CCDimage
import datetime as DT
from matplotlib import pyplot as plt
# Open the NetCDF file as an Xarray dataset

#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(6, 1, figsize=(10, 15))
fig2, ax2 = plt.subplots(6, 1, figsize=(10, 15))
for i, channel in enumerate(['UV2', 'UV1', 'IR1', 'IR2']):

    dsfile='/Users/lindamegner/MATS/testdatadownload/'+channel+'_zarr.nc'
    ds = xr.open_dataset(dsfile)

    print('This is channel:',channel)
    print(ds)

    start_time = DT.datetime(2023, 3, 8, 9, 0, 0)
    stop_time = DT.datetime(2023, 3, 8, 9, 15, 0)
    image_mean = ds['ImageCalibrated'].sel(time=slice(start_time, stop_time)).mean(dim='time')
    #label axis

    image_mean.plot(ax=ax[i])
    ax[i].set_xlabel(image_mean.dims[1])
    ax[i].set_ylabel(image_mean.dims[0])
    ax[i].set_title(channel)

    attrs = ds['ImageCalibrated'].attrs
    sp=plot_CCDimage(image_mean, fig=fig2, axis=ax2[i], title=channel+attrs.get("long_name", "Calibrated Image"))


#%%
channel='IR2'
dsfile='/Users/lindamegner/MATS/testdatadownload/'+channel+'_zarr.nc'
ds = xr.open_dataset(dsfile)
start_time = DT.datetime(2023, 3, 8, 9, 17, 0)
stop_time = DT.datetime(2023, 3, 8, 9, 18, 0)
images = ds['ImageCalibrated'].sel(time=slice(start_time, stop_time))
for i in images[:]:
    print(i)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    #ax.imshow(i, cmap='inferno')

    i.plot(ax=ax)
    #plot_CCDimage(images.mean(dim='time'), title=channel+ds['ImageCalibrated'].attrs.get("long_name", "Calibrated Image"))
#%%
# compare to aws data
from mats_utils.rawdata.read_data import read_MATS_data
filter={'channel': 'UV1'}
df = read_MATS_data(start_time, stop_time,filter, level='1a',version='1.0')
#%%
for i in range (0,8):
    plot_CCDimage(df['IMAGE'][i], title=df['EXPDate'][i])


#%%
first_image =  ds['ImageCalibrated'].isel(time=0)
# Plot the first image using the attributes for labels and title
plt.figure(figsize=(10, 6))
first_image.plot()
plt.title(f"{ds['ImageCalibrated'].attrs.get('long_name', 'Calibrated Image')} at {ds['time'].values[0]}")
plt.xlabel("Image Column")
plt.ylabel("Image Row")
plt.show()



#%%
channel='IR2'
dsfile='/Users/lindamegner/MATS/testdatadownload/'+channel+'_zarr.nc'

ds = xr.open_dataset(dsfile)
# Assuming 'dataset' is your xarray Dataset
image_calibrated = ds["ImageCalibrated"]
# Select the first and last images
first_image = image_calibrated.isel(time=0)
last_image = image_calibrated.isel(time=-1)

# Plot the images
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

first_image.plot(ax=ax1)
ax1.set_title("First Image")
ax1.set_xlabel(image_calibrated.dims[2])
ax1.set_ylabel(image_calibrated.dims[1])

last_image.plot(ax=ax2)
ax2.set_title("Last Image")
ax2.set_xlabel(image_calibrated.dims[2])
ax2.set_ylabel(image_calibrated.dims[1])

plt.tight_layout()
plt.show()

# %%
#lukas script
import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import sys


channel = sys.argv[2]


with nc.Dataset(f"{sys.argv[1]}", 'r') as nf:
    images = nf["ImageCalibrated"][:]
    #channel = nf["channel"]

for img in [int(x) for x in sys.argv[3:]]:
    plt.figure(figsize=(12, 8))
    assert not (images[img, :, :].mask).all(), "Whole image is masked!"
    vrange = [np.percentile(images[img, :, :].data, p) for p in [1, 99]]
    plt.imshow(images[img, :, :], cmap="inferno", vmin=vrange[0], vmax=vrange[1], aspect='auto', origin='lower')
    plt.ylabel("Row number")
    plt.xlabel("Column number")
    plt.title(channel)
    plt.colorbar()
    plt.savefig(f"{channel}_{img}.png", dpi=400)
# %%
