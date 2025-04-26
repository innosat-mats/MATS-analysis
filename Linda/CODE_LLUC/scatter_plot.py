import pickle
import numpy as np

def cargar_datos(nombre_archivo):

    """
    Function that loads data from a file
    Args:
    - nombre archivo: file path
    Returns:
    - data array
    """

    datos = np.load(nombre_archivo)['datos']
    return datos

heights = cargar_datos('/Users/lluccampinssastre/Documents/MATS/rectangle_good_data/proportions_3_days_240430.npz')
latitudes = cargar_datos('/Users/lluccampinssastre/Documents/MATS/rectangle_good_data/latitudes_3_days_240430.npz')
longitudes = cargar_datos('/Users/lluccampinssastre/Documents/MATS/rectangle_good_data/longitudes_3_days_240430.npz')


# Approximate latitudes and longitudes to the nearest integer.
rounded_latitudes = np.round(latitudes).astype(int)
rounded_longitudes = np.round(longitudes).astype(int)

heights_dict = {}

# Group heights by latitude and longitude.
for lat, lon, height in zip(rounded_latitudes, rounded_longitudes, heights):
    key = (lat, lon)
    if key not in heights_dict:
        heights_dict[key] = []
    heights_dict[key].append(height)

# Calculate the mean of heights for each latitude and longitude.
mean_heights = []
for lat, lon in zip(rounded_latitudes, rounded_longitudes):
    key = (lat, lon)
    height_group = heights_dict[key]
    mean_height = np.mean(height_group)
    mean_heights.append(mean_height)


import matplotlib.pyplot as plt
import cartopy.crs as ccrs


fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

sc = ax.scatter(rounded_longitudes, rounded_latitudes, c=mean_heights, cmap='viridis',  vmin=0, vmax=1, linewidths = 0.01, transform=ccrs.PlateCarree())
#sc = ax.scatter(rounded_longitudes, rounded_latitudes, c=mean_heights, cmap='viridis', linewidths = 0.01, transform=ccrs.PlateCarree())  #without max value

cbar = plt.colorbar(sc, ax=ax, orientation='vertical', label='Proportion good data', ticks = np.arange(0,1,0.1))

ax.coastlines()
ax.gridlines()

plt.show()

