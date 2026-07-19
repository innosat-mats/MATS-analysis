
#%%

import xarray as xr
from database_generation.experimental_utils import plot_CCDimage
import datetime as DT
from matplotlib import pyplot as plt

#%%
dsfile='/Users/lindamegner/MATS/MATS-retrieval/data/MATS_L1b_IR1_20230210T180607_20230210T190607.nc'
# Open the NetCDF file as an Xarray dataset
ds = xr.open_dataset(dsfile)


image_calibrated = ds["ImageCalibrated"]
# Select the first and last images
first_image = image_calibrated.isel(time=0)
last_image = image_calibrated.isel(time=-1)

# Plot the images
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

first_image.plot(ax=ax1)
first_time = str(ds["time"].isel(time=0).values)
first_tplat = float(ds["TPlat"].isel(time=0).values)
first_tplon = float(ds["TPlon"].isel(time=0).values)
ax1.set_title(f"First Image | Time: {first_time} | TPlat: {first_tplat:.4f} | TPlon: {first_tplon:.4f}")
ax1.set_xlabel(image_calibrated.dims[2])
ax1.set_ylabel(image_calibrated.dims[1])

last_image.plot(ax=ax2)
last_time = str(ds["time"].isel(time=-1).values)
last_tplat = float(ds["TPlat"].isel(time=-1).values)
last_tplon = float(ds["TPlon"].isel(time=-1).values)
ax2.set_title(f"Last Image | Time: {last_time} | TPlat: {last_tplat:.4f} | TPlon: {last_tplon:.4f}")
ax2.set_xlabel(image_calibrated.dims[2])
ax2.set_ylabel(image_calibrated.dims[1])

plt.tight_layout()
plt.show()
# %%
