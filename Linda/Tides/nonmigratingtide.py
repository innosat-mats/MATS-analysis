#%%
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
# Open the Zarr file as an Xarray dataset

#Read from Bolin
#data=xr.open_zarr('zip::https://bolin.su.se/data/s31/data/mats-level-1b-limb-cropd-1.0-root/mats-level-1b-limb-cropd-UV2.zarr.zip')
#Read from downloaded s3 bucket

#data=xr.open_zarr('s3://test-release-v1.0.1/mats-level-1b-limb-cropd-UV2.zarr/')
#Read from local file
#dsfile='/Users/lindamegner/MATS/testdatadownload/mats-level-1b-limb-cropd-UV2-rev-2.zarr'
#data=xr.open_zarr(dsfile)


channels = ['UV2', 'UV1', 'IR1', 'IR2', 'IR3', 'IR4']
channel='IR2'
dataall=xr.open_zarr('s3://test-release-v1.0.1/mats-level-1b-limb-cropd-'+channel+'.zarr/',storage_options={"profile":"mats"})

dataall.attrs['long_name'] = channel
# Print data structure
print(dataall)

#%%

data = dataall.sel(time=slice("2023-02-11", "2023-02-22"))



data['TPlon360'] = data.TPlon + 180
data['DW3phase'] = xr.where(data.TPlon360 <= 120, data.TPlon360,
                            xr.where(data.TPlon360 <= 240, data.TPlon360 - 120, data.TPlon360 - 240))



# %%
data['mean_signal'] = data.ImageCalibrated.mean(dim='im_col')


#%%
datafull = data.copy()
#%%
dataeq = data.where((data.TPlat >= -5) & (data.TPlat <= 5))
datamidlatn = data.where((data.TPlat >= 30) & (data.TPlat <= 60))
datamidlats = data.where((data.TPlat >= -60) & (data.TPlat <= -30))
datanorth = data.where(data.TPlat >= 60)
datasouth = data.where(data.TPlat <= -60)




#%%
data=datamidlatn
# Define the bin edges for im_row and DW3phase
im_row_bins = np.arange(data.im_row.min(), data.im_row.max() + 1, 1)
DW3phase_bins = np.arange(np.floor(data.DW3phase.min()), np.ceil(data.DW3phase.max()) + 1, 1)

# Digitize the im_row and DW3phase data into bins
im_row_binned = xr.DataArray(np.digitize(data.im_row, im_row_bins) - 1, dims='im_row')
DW3phase_binned = xr.DataArray(np.digitize(data.DW3phase, DW3phase_bins) - 1, dims='time')

# Add the binned data as data variables
data = data.assign(im_row_binned=im_row_binned, DW3phase_binned=DW3phase_binned)

# Ensure the binned variables are accessible
data = data.set_coords(['im_row_binned', 'DW3phase_binned'])

# Group by the binned coordinates and calculate the mean
mean_signal_binned = data.mean_signal.groupby(['im_row_binned', 'DW3phase_binned']).mean()

print(mean_signal_binned)
#%%



fig, ax = plt.subplots(figsize=(10, 6))
# plot and set color limits
mean_signal_binned.plot(ax=ax)
ax.set_xlabel('DW3phase')
ax.set_ylabel('im_row')
ax.set_title('Mean Signal Strength by DW3phase and im_row')

plt.show()

fig2, ax2 = plt.subplots(figsize=(10, 6))
# get the mean value between im_row 90 and 100
mean_at_max=mean_signal_binned.sel(im_row_binned=slice(90, 100)).mean(dim='im_row_binned')
mean_at_max.plot(ax=ax2)





#%%
# Define the bin edges for im_row and DW3phase
im_row_bins = np.arange(data.im_row.min(), data.im_row.max() + 1, 1)
DW3phase_bins = np.arange(np.floor(data.DW3phase.min()), np.ceil(data.DW3phase.max()) + 1, 1)

# Digitize the im_row and DW3phase data into bins
data['im_row_binned'] = xr.DataArray(np.digitize(data.im_row, im_row_bins) - 1, dims='im_row')
data['DW3phase_binned'] = xr.DataArray(np.digitize(data.DW3phase, DW3phase_bins) - 1, dims='time')

# Group by the binned coordinates and calculate the mean
#mean_signal_binned = data.mean_signal.groupby(['im_row_binned', 'DW3phase_binned']).mean()



# Add the binned data to the dataset
data = data.assign_coords(im_row_binned=('im_row_binned', im_row_bins[data.im_row_binned]))
data = data.assign_coords(DW3phase_binned=('DW3phase_binned', DW3phase_bins[data.DW3phase_binned]))

# Group by the binned coordinates and calculate the mean
mean_signal_binned = data.mean_signal.groupby(['im_row_binned', 'DW3phase_binned']).mean()

print(mean_signal_binned)
#%%
# Plot the mean signal strength as a function of im_row and DW3phase
plt.figure(figsize=(10, 6))
plt.pcolor(data.DW3phase, data.im_row, mean_signal.T, cmap='viridis')
plt.colorbar(label='Mean Signal Strength')
plt.xlabel('DW3phase')
plt.ylabel('im_row')
plt.title('Mean Signal Strength by DW3phase and im_row')
plt.show()

#%%
# Create a 2D histogram of mean_signal with im_row on the y-axis and DW3phase on the x-axis
mean_signal_hist = np.histogram2d(data.DW3phase, data.im_row, bins=(100, 100), weights=mean_sig
# %%
import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data
x = np.random.rand(1000) * 120  # DW3phase values
y = np.random.rand(1000) * 180  # im_row values

# Create a 2D histogram
hist, xedges, yedges = np.histogram2d(x, y, bins=[120, 180])

# Plot the 2D histogram using pcolormesh
plt.figure(figsize=(10, 6))
plt.pcolormesh(xedges, yedges, hist.T, cmap='viridis')
plt.colorbar(label='Count')
plt.xlabel('DW3phase')
plt.ylabel('im_row')
plt.title('2D Histogram of DW3phase and im_row')
plt.show()
# %%
