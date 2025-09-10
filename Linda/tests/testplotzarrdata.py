#%%
import xarray as xr
# Open the Zarr file as an Xarray dataset

#Read from Bolin
#data=xr.open_zarr('zip::https://bolin.su.se/data/s31/data/mats-level-1b-limb-cropd-1.0-root/mats-level-1b-limb-cropd-UV2.zarr.zip')
#Read from downloaded s3 bucket
data=xr.open_zarr('s3://test-release-v1.0.1/mats-level-1b-limb-cropd-UV2-rev-2.zarr/')
#Read from local file
#dsfile='/Users/lindamegner/MATS/testdatadownload/mats-level-1b-limb-cropd-UV2-rev-2.zarr'
data=xr.open_zarr(dsfile)

# Print data structure
print(data)
# %%
meantpheight=data["TPheight"].mean(dim="time")
print(meantpheight)
meaniamge=data["ImageCalibrated"].mean(dim="time")

# %%
