
# %%
import xarray as xr
from database_generation.experimental_utils import plot_CCDimage
# Open the NetCDF file as an Xarray dataset
UV2='/Users/lindamegner/MATS/MATS-retrieval/data/UV2_zarr.nc'
UV2_day='/Users/lindamegner/MATS/testdatadownload/UV2_zarr_day.nc'
ds = xr.open_dataset(UV2_day)

# Access the 'temperature' variable
temperature = ds['temperature']
print(temperature)

# Access the 'ImageCalibrated' variable
image_calibrated = ds['ImageCalibrated']
#%%
# find images where 
# plot the mean value of the images as a fuction of time
image_mean = image_calibrated.mean(dim=('time'))
#set colorscale limits to -500 to 500
image_mean.plot(vmin=-100, vmax=100)


#plot_CCDimage(image_calibrated[3, :, :])


# %%
import xarray as xr
# Open the Zarr file as an Xarray dataset

xr.open_zarr('zip::https://bolin.su.se/data/s31/data/mats-level-1b-limb-cropd-1.0-root/mats-level-1b-limb-cropd-UV2.zarr.zip')
# %%
import xarray as xr

# Assign web address
path = 'zip::https://bolin.su.se/data/s31/data/mats-level-1b-limb-cropd-1.0-root/mats-level-1b-limb-cropd-UV2.zarr.zip'

# Get xarray Dataset using zarr engine
data = xr.open_zarr(path)

# Print data structure
print(data)

# Retrieve first time step
time = data.time[0].values

# Retrieve first calibrated image
image = data.ImageCalibrated[:,:,0].values

file='/Users/lindamegner/MATS/testdatadownload/mats-level-1b-limb-cropd-UV2.zarr.zip'
data2 = xr.open_zarr(file)
print(data2)
# Retrieve first time step
time = data2.time[0].values

# Retrieve first calibrated image
image = data2.ImageCalibrated[:,:,0].values
# %%
