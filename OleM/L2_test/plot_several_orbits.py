#%%
import numpy as np
import pandas as pd
from mats_l2_processing.forward_model import calc_jacobian, cart2sph
import pickle
from datetime import datetime
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l2_processing.inverse_model import do_inversion
from mats_l2_processing.grids import center_grid,localgrid_to_lat_lon_alt_3D
from matplotlib import pyplot as plt

filename = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/0.pkl"

with open(filename, "rb") as file:
    result = pickle.load(file)
#%%
[x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local] = result

# %%
altitude_grid = center_grid(altitude_grid_edges)
alongtrack_grid = center_grid(alongtrack_grid_edges)
acrosstrack_grid = center_grid(acrosstrack_grid_edges)
x_hat_reshaped = np.array(x_hat).reshape(len(altitude_grid),len(acrosstrack_grid),len(alongtrack_grid))

[alt,lat,lon,r] = localgrid_to_lat_lon_alt_3D(altitude_grid,alongtrack_grid,acrosstrack_grid,ecef_to_local)
lat = np.rad2deg(lat)
lon = np.rad2deg(lon)

# %%
[alt,lat,lon,r] = localgrid_to_lat_lon_alt_3D(altitude_grid,alongtrack_grid,acrosstrack_grid,ecef_to_local)
lat = np.rad2deg(lat)
lon = np.rad2deg(lon)

#%%
import numpy as np
from scipy.interpolate import griddata

# Define your 3D data points and corresponding values
data_points = np.array([alt.flatten(),lon.flatten(),lat.flatten()]).T  # Shape: (N, 3)
values = x_hat       # Shape: (N,)

# Define the regular grid
alt_grid = np.arange(alt.min(), alt.max(), 1e3)
lon_grid = np.arange(lon.min(), lon.max(), 0.5)
lat_grid = np.arange(lat.min(), lat.max(), 0.5)
x_grid, y_grid, z_grid = np.meshgrid(alt_grid, lon_grid, lat_grid, indexing='ij')
interpgrid = np.array([x_grid.flatten(),y_grid.flatten(),z_grid.flatten()]).T
# Interpolate using linear interpolation
interpolated_values = griddata(data_points, values, interpgrid, method='linear')
# %%
interpolated_values = interpolated_values.reshape(len(alt_grid),len(lon_grid),len(lat_grid))
# %%
