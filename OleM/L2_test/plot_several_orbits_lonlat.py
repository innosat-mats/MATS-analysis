#%%
import numpy as np
import pandas as pd
from mats_l2_processing.forward_model import calc_jacobian
import pickle
from datetime import datetime
from mats_utils.rawdata.read_data import read_MATS_data
from mats_l2_processing.inverse_model import do_inversion
from mats_l2_processing.grids import center_grid,localgrid_to_lat_lon_alt_3D, sph2cart, cart2sph
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import xarray as xr
folder = "/home/olemar/Projects/Universitetet/MATS/MATS-analysis/L2_test/"
filenames = ["0","50","100","150","200","250","300","350","400","450","500","550","600","650","700","750","800"]
#filenames = ["500","550","600","650","700","750","800"]

all_data = []
all_alt_grid = []
all_lon_grid = []
all_lat_grid = []

#%% Define the regular grid
alt_grid = np.arange(70e3, 110e3, 2e3)
#lon_grid = np.arange(-90, 0, 0.02)
#lat_grid = np.arange(-90, 0, 0.1)
along_grid = np.arange(-0.5, 0.5, 0.01)
across_grid = np.arange(-0.05, 0.05, 0.005)

all_interpolated = np.zeros([len(alt_grid),len(across_grid),len(along_grid),len(filenames)])
#all_interpolated = np.zeros([len(alt_grid),len(lon_grid),len(lat_grid),len(filenames)])

all_lat = np.zeros([len(alt_grid),len(across_grid),len(along_grid),len(filenames)])
all_lon = np.zeros([len(alt_grid),len(across_grid),len(along_grid),len(filenames)])

#%%
for i, filename in enumerate(filenames):
    with open(folder + filename + ".pkl", "rb") as file:
        result = pickle.load(file)
    [x_hat, y, ks, altitude_grid_edges, alongtrack_grid_edges,acrosstrack_grid_edges, ecef_to_local] = result

    radius_grid = center_grid(altitude_grid_edges)
    acrosstrack_grid = center_grid(acrosstrack_grid_edges)
    alongtrack_grid = center_grid(alongtrack_grid_edges)
    x_hat_reshaped = np.array(x_hat).reshape(len(radius_grid),len(acrosstrack_grid),len(alongtrack_grid))

    crop_r = np.array([6.4,6.5])*1e6
    crop_across = np.array([-0.03,0.03])
    crop_along = np.array([-0.5,0.5])

    condition = lambda x: (x > crop_r[0]) and (x < crop_r[1])
    r_index = [index for index, value in enumerate(radius_grid) if condition(value)]

    condition = lambda x: (x > crop_across[0]) and (x < crop_across[1])
    acrosstrack_index = [index for index, value in enumerate(acrosstrack_grid) if condition(value)]

    condition = lambda x: (x > crop_along[0]) and (x < crop_along[1])
    alongtrack_index = [index for index, value in enumerate(alongtrack_grid) if condition(value)]

    radius_grid = radius_grid[r_index]
    acrosstrack_grid = acrosstrack_grid[acrosstrack_index]
    alongtrack_grid = alongtrack_grid[alongtrack_index]
    x_hat_reshaped = x_hat_reshaped[r_index,:,:][:,acrosstrack_index,:][:,:,alongtrack_index]

    [alt,lon,lat,r] = localgrid_to_lat_lon_alt_3D(radius_grid,acrosstrack_grid,alongtrack_grid,ecef_to_local)
    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon)

    # Define your 3D data points and corresponding values
    _, acrosstrack_grid_expanded, alongtrack_grid_expanded = np.meshgrid(radius_grid, acrosstrack_grid, alongtrack_grid, indexing='ij')

    data_points = np.array([alt.flatten(),acrosstrack_grid_expanded.flatten(),alongtrack_grid_expanded.flatten()]).T  # Shape: (N, 3)
    #data_points = np.array([alt.flatten(),lon.flatten(),lat.flatten()]).T  # Shape: (N, 3)
    values = x_hat_reshaped.flatten()       # Shape: (N,)

    x_grid, y_grid, z_grid = np.meshgrid(alt_grid, across_grid,along_grid, indexing='ij')
    # x_grid, y_grid, z_grid = np.meshgrid(alt_grid,lon_grid,lat_grid, indexing='ij')
    interpgrid = np.array([x_grid.flatten(),y_grid.flatten(),z_grid.flatten()]).T
    # Interpolate using linear interpolation
    interpolated_values = griddata(data_points, values, interpgrid, method='linear')
    interpolated_lat = griddata(data_points, lat.flatten(), interpgrid, method='linear')
    interpolated_lon = griddata(data_points, lon.flatten(), interpgrid, method='linear')

    interpolated_values = interpolated_values.reshape(len(alt_grid),len(across_grid),len(along_grid))
    interpolated_lon = interpolated_lon.reshape(len(alt_grid),len(across_grid),len(along_grid))
    interpolated_lat = interpolated_lat.reshape(len(alt_grid),len(across_grid),len(along_grid))

    #interpolated_values = interpolated_values.reshape(len(alt_grid),len(lon_grid),len(lat_grid))

    all_lon[:,:,:,i] = interpolated_lon
    all_lat[:,:,:,i] = interpolated_lat
    all_interpolated[:,:,:,i] = interpolated_values
# %%
for i in range(all_interpolated.shape[3]):
    plt.pcolor(along_grid-0.4*i,alt_grid,all_interpolated[:,10,:,i],vmin=0,vmax=3e12)

plt.show()
# %%
ds = xr.Dataset(
    data_vars=dict(
        VER=(["alt", "across", "along","batch"], all_interpolated),
    ),
    coords=dict(
        lon=(["alt", "across", "along","batch"], all_lon),
        lat=(["alt", "across", "along","batch"], all_lat),
        alt=alt_grid,
        across=across_grid,
        along=along_grid,
        batch=np.arange(len(filenames)),
    ),)
# %%
lat_plot = ds.lat.sel(alt=90e3)
VER_plot = ds.VER.sel(alt=90e3)
lon_plot = ds.lon.sel(alt=90e3)

plt.scatter(lon_plot.values.flatten(),lat_plot.values.flatten(),c=VER_plot.values.flatten(),s=1)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('First tomo orbit from MATS')
# %%
