
#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
#read in gloassac xarray file
filename="/Users/lindamegner/MATS/MATS-retrieval/data/GloSSAC-V1.1.nc"

#%%
#ds=xr.open_dataset(filename)
data = xr.open_dataset(filename, decode_times=False)

# xarray.DataArray'DerivedProducts_N1'altitude1: 70latitude: 32years: 456
nd=data.DerivedProducts_N1
# mean over years
nd_mean=nd.mean(dim='years')

# print out latitudes available
print(data.lat.values)
#[-77.5 -72.5 -67.5 -62.5 -57.5 -52.5 -47.5 -42.5 -37.5 -32.5 -27.5 -22.5
# -17.5 -12.5  -7.5  -2.5   2.5   7.5  12.5  17.5  22.5  27.5  32.5  37.5
#  42.5  47.5  52.5  57.5  62.5  67.5  72.5  77.5]

#%%
# Select latitude 70 deg N
nd_mean_latselect=nd_mean.sel(latitude=29)

#%%
# plot nd vs altitude (altitude in km on y-axis)
plt.figure()
plt.plot(nd_mean_latselect.values, nd_mean_latselect.altitude1.values)
plt.xlabel('Aerosol Number Density (cm^-3)')
plt.ylabel('Altitude (km)')
plt.title('GloSSAC Aerosol Number Density at 62 deg N (mean over years)')
plt.xscale('log')
plt.grid(True)
plt.gca().invert_yaxis()
plt.show()

# mean over 
# %%
